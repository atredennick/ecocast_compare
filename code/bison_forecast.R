##  R script to fit a population growth model for YNP,
##  forecast 10 new years, and partition the forecast variance.
##
##  Based on Dietze et al. (2017)
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: October 19, 2016
##


rm(list=ls(all.names = TRUE))


####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(ggthemes)  # Pleasing themese for ggplot2
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
library(ggmcmc)    # MCMC-to-dataframe functions
library(stringr)
library(cowplot)
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze; modified by A. Tredennick)



####
####  LOAD DATA ----------------------------------------------------------------
####
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1)
weather_dat <- read.csv("../data/PRISM_ppt_tmin_tmean_tmax_tdmean_vpdmin_vpdmax_provisional_4km_197001_201711_44.8090_-110.5728.csv", skip = 10)
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")

##  Reformat weather data from PRISM
weather_dat <- weather_dat %>%
  dplyr::select(-tdmean..degrees.F., -vpdmin..hPa., -vpdmax..hPa.) %>%
  dplyr::rename(date = Date,
                ppt_in = ppt..inches.,
                tmin_F = tmin..degrees.F.,
                tmean_F = tmean..degrees.F.,
                tmax_F = tmax..degrees.F.) %>%
  separate(date, into = c("year", "month"), sep = "-")

precip_dat <- weather_dat %>%
  dplyr::select(year, month, ppt_in) %>%
  filter(month %in% c("01")) %>%
  mutate(year = as.integer(year))

##  Reformat bison data and combine with weather data
bison_dat <- bison_raw %>% 
  dplyr::select(-ends_with("source")) %>%     # drop the source column
  mutate(set = ifelse(year < 2011, "training", "validation")) %>% # make new column for data splits
  left_join(snow_ynp, by="year") %>% # merge in SNOTEL data
  left_join(precip_dat, by="year")



####
####  JAGS State-Space Model ---------------------------------------------------
####
r_mu_prior <- log(1.11) # lambda = 1.11 in Hobbs et al. 2015
r_sd_prior <- sd(log(rnorm(100000,1.11,0.024))) # sd_lambda = 0.024 in Hobbs et al. 2015

my_model <- "  
model{

  #### Variance Priors
  sigma_proc ~ dunif(0,5)
  tau_proc   <- 1/sigma_proc^2
  eta        ~ dunif(0, 50)
  
  #### Fixed Effects Priors
  r  ~ dnorm(0.1, 1/0.02^2)  # intrinsic growth rate, informed prior
  b  ~ dnorm(0,1/2^2)I(-2,2) # strength of density dependence, bounded
  b1 ~ dnorm(0,0.0001)       # effect of Jan. precip in year t
  
  #### Initial Conditions
  z[1]    ~ dnorm(Nobs[1], tau_obs[1]) # varies around observed abundance at t = 1
  zlog[1] <- log(z[1]) # set first zlog
  
  #### Process Model
  for(t in 2:npreds){
    # Calculate log integration of extractions
      e[t] = log( abs( 1 - (E[t] / z[t-1]) ) ) 

      # Gompertz growth, on log scale
      mu[t]   <- zlog[t-1] + e[t] + r + b*(zlog[t-1] + e[t]) + b1*x[t]
      zlog[t] ~ dnorm(mu[t], tau_proc)
      z[t]    <- exp(zlog[t]) # back transform to arithmetic scale
  }
  
  #### Data Model
  for(j in 2:n){
    p[j]     <- eta/(eta + z[j]) # calculate NB centrality parameter
    Nobs[j]  ~ dnegbin(p[j], eta) # NB likelihood
  }

}"



####
####  Fit Bison Forecasting Model ----------------------------------------------
####
##  For years without observation error, set to max observed standard deviation
na_sds <- which(is.na(bison_dat$count.sd)==T)
bison_dat[na_sds,"count.sd"] <- max(bison_dat$count.sd, na.rm=T)

##  Split into training and validation sets
training_dat   <- filter(bison_dat, set == "training")
validation_dat <- filter(bison_dat, set == "validation")

##  Set up SWE knowns (2011-2017), relative to scaling of observations
ppt_mean     <- mean(training_dat$ppt_in)
ppt_sd       <- sd(training_dat$ppt_in)
forecast_ppt <- precip_dat %>%
  filter(year %in% validation_dat$year) %>%
  pull(ppt_in)
scl_fut_ppt  <- (forecast_ppt - ppt_mean) / ppt_sd

##  Set initial values for unkown parameters
inits <- list(
  list(sigma_proc = 0.01,
       r = 0.05,
       b = -0.001,
       b1 = -0.5),
  list(sigma_proc = 0.3,
       r = 0.4,
       b = -0.1,
       b1 = -0.01),
  list(sigma_proc = 0.1,
       r = 0.7,
       b = -0.00001,
       b1 = -0.2)
)



####
####  FIT AND FORECAST WITH KNOWN JAN PPT --------------------------------------
####
##  Prepare data list
extractions <- training_dat$wint.removal
extractions[is.na(extractions) == TRUE] <- 0

mydat <- list(Nobs    = round(training_dat$count.mean), # mean counts
              E       = c(extractions,rep(0,7)),
              n       = nrow(training_dat), # number of observations
              tau_obs = 1/training_dat$count.sd^2, # transform s.d. to precision
              x       = c(as.numeric(scale(training_dat$ppt_in)),scl_fut_ppt), # snow depth, plus forecast years
              npreds  = nrow(training_dat)+nrow(validation_dat)) # number of total predictions (obs + forecast)

##  Random variables to collect
out_variables <- c("r", "b", "b1", "sigma_proc", "z")

##  Send to JAGS
mc3     <- jags.model(file = textConnection(my_model), 
                      data = mydat, 
                      n.chains = length(inits), 
                      n.adapt = 5000, 
                      inits = inits) 
update(mc3, n.iter = 10000) 
mc3.out <- coda.samples(model=mc3, 
                        variable.names=out_variables, 
                        n.iter=10000) 




## Split output
out          <- list(params=NULL, predict=NULL)
mfit         <- as.matrix(mc3.out,chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
predictions        <- rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = bison_dat$year,
                                 set = bison_dat$set,
                                 observation = bison_dat$count.mean,
                                 upper_observation = bison_dat$count.mean+bison_dat$count.sd,
                                 lower_observation = bison_dat$count.mean-bison_dat$count.sd,
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)



####
####  PLOT DATA AND POSTERIOR PREDICTIONS --------------------------------------
####
pred_color <- "#CF4C26"
obs_color  <- "#278DAF"
pred_color <- "black"
obs_color  <- "black"
calibration_plot <- ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction),
              fill=pred_color, 
              color=NA, 
              alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color, size = 0.2)+
  geom_errorbar(data = filter(prediction_df, year <2011), aes(x = year, ymin=lower_observation, ymax=upper_observation), 
                width=0.5, 
                color=obs_color, 
                size=0.2)+
  geom_point(data = filter(prediction_df, year <2011), aes(x = year, y=observation), color=obs_color, size=0.5)+
  geom_vline(aes(xintercept=2010), linetype=2,color="grey55")+
  # geom_col(data = bison_dat, aes(x = year, y = wint.removal), color = "grey55", fill = "grey55", width = 0.3)+
  scale_y_continuous(breaks = seq(0,15000,2500))+
  # scale_x_continuous(breaks = seq(1970,2015,5))+
  ylab("Number of bison")+
  xlab("Year")+
  theme_few()



####
####  SET UP GCM PROJECTION MATRIX ---------------------------------------------
####
# Set up column names for GCM projection file
col_names <- c("year",
               "month",
               "day",
               as.character(as.data.frame(read.table("../data/CMIP_YNP/bcca5/COLS_pr.txt"))[,1])
)

# Read in GCM projections and format as matrix
gcm_precip <- read_csv("../data/CMIP_YNP/bcca5/pr.csv", col_names = col_names) %>%
  gather(key = model, value = ppt, -year, -month, -day) %>%
  separate(model, into = c("model_name", "rep", "scenario"), sep = "[.]") %>%
  group_by(year, month, model_name, scenario) %>%
  summarise(total_ppt_mm = sum(ppt),
            total_ppt_in = total_ppt_mm*0.0393701) %>%
  ungroup() %>%
  dplyr::filter(month == 1) %>%
  dplyr::filter(year %in% c(2011,2012,2013,2014,2015,2016,2017)) %>%
  dplyr::select(model_name, scenario, year, month, total_ppt_mm, total_ppt_in) %>%
  dplyr::arrange(model_name, scenario, year, month) %>%
  dplyr::mutate(stdzd_precip = (total_ppt_in-ppt_mean) / ppt_sd) %>%
  dplyr::mutate(model_rcp = paste(model_name, scenario,"::")) %>%
  dplyr::select(model_rcp, year, stdzd_precip) %>%
  spread(model_rcp, stdzd_precip)



####
####  PARTITION FORECAST UNCERTAINTY -------------------------------------------
####
##  Function for the ecological process (Gompertz population growth)
iterate_process <- function(Nnow, xnow, r, b, b1, sd_proc) { 
  xnow[xnow>5] <- 5
  mu <- log(Nnow) + r + b*log(Nnow) + b1*xnow # determinstic process; log scale
  zlog <- rnorm(length(mu), mu, sd_proc) # stochastic process; log scale
  N <- exp(zlog) # back transform to arithmetic scale
}


##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values and no process error.
forecast_steps <- 7
num_iters      <- 50000
x              <- sample(predictions[,nrow(training_dat)], num_iters, replace = TRUE)
param_summary  <- summary(fitted_model$params)$quantile
r              <- param_summary[3,3]
b              <- param_summary[1,3]
b1             <- param_summary[2,3]
sd_proc        <- param_summary[4,3]
z              <- scl_fut_ppt
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, xnow = z[t], r, b, b1, sd_proc = 0)
  forecasts[,t] <- x
}
varI <- apply(forecasts,2,var)


##  Initial conditions and parameter uncertainty
x              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
r              <- params[sample_params,"r"]
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
sd_proc        <- param_summary[4,3]
z              <- scl_fut_ppt
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, xnow = z[t], r, b, b1, sd_proc = 0)
  forecasts[,t] <- x
}
varIP <- apply(forecasts,2,var)


##  Initial conditions, parameter, and driver uncertainty
x              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
r              <- params[sample_params,"r"]
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
sd_proc        <- params[sample_params,"sigma_proc"]
zsamps         <- sample(x = ncol(gcm_precip[2:ncol(gcm_precip)]), size = num_iters, replace = TRUE)
z              <- as.matrix(gcm_precip[2:ncol(gcm_precip)])
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, xnow = as.numeric(z[t,zsamps]), r, b, b1, sd_proc = 0)
  forecasts[,t] <- x
}
varIPD <- apply(forecasts,2,var)

##  Initial conditions, parameter, driver, and process uncertainty
x              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
r              <- params[sample_params,"r"]
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
sd_proc        <- param_summary[4,3]
z              <- as.matrix(gcm_precip[2:ncol(gcm_precip)])
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, xnow = as.numeric(z[t,zsamps]), r, b, b1, sd_proc = sd_proc)
  forecasts[,t] <- x
}
varIPDE <- apply(forecasts,2,var)


V.pred.sim     <- rbind(varIPDE,varIPD,varIP,varI)
V.pred.sim.rel <- apply(V.pred.sim,2,function(x) {x/max(x)})



####
####  PLOT THE FORECASTING UNCERTAINTY PARTITION -------------------------------
####
var_rel_preds <- as.data.frame(t(V.pred.sim.rel*100))
var_rel_preds$x <- 1:nrow(var_rel_preds)
my_cols <- c("#0A4D5B", "#139AB8", "#39B181","grey")
my_cols <- c("black", "grey55", "grey70","grey90")
variance_plot <- ggplot(data=var_rel_preds, aes(x=x))+
  geom_ribbon(aes(ymin=0, ymax=varIPDE), fill=my_cols[4])+
  geom_ribbon(aes(ymin=0, ymax=varIPD), fill=my_cols[3])+
  geom_ribbon(aes(ymin=0, ymax=varIP), fill=my_cols[2])+
  geom_ribbon(aes(ymin=0, ymax=varI), fill=my_cols[1])+
  ylab("Percent of uncertainty")+
  xlab("Forecast steps")+
  scale_x_continuous(breaks=seq(1,forecast_steps,by=1), 
                     labels=paste(seq(1,forecast_steps,by=1), "yrs"))+
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))+
  theme_few()



####
####  COMBINE PLOTS AND SAVE ---------------------------------------------------
####
plot_grid(calibration_plot, variance_plot, nrow = 2, labels = "AUTO")
ggsave(filename = "../figures/bison_combined.png",
       width = 4,
       height = 6,
       units = "in",
       dpi =200)