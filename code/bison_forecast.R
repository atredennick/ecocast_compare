################################################################################
##  bison_forecast.R: R script to fit a population growth model for YNP Bison,
##  forecast 10 new years, and partition the forecast variance.
##
##  Based on Dietze 2017, Ecological Applications
##  http://onlinelibrary.wiley.com/doi/10.1002/eap.1589/full
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: October 19, 2016
################################################################################

##  Clear everything...
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location...
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # only for RStudio



####
####  LOAD LIBRARIES ----
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(ggthemes)  # Pleasing themes for ggplot2
library(cowplot)   # Combining ggplots
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze)



####
####  SET MY PLOTTING THEME ----------------------------------------------------
####
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="white"),
        panel.background   = element_rect(fill = "#EFEFEF"),
        axis.text          = element_text(size=10, color="grey35", family = "Arial Narrow"),
        axis.title         = element_text(size=12, family = "Arial Narrow", face = "bold"),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color="black"),
        axis.line.y        = element_line(color="black"),
        strip.background   = element_blank(),
        strip.text         = element_text(size=10, color="grey15", family = "Arial Narrow"))



####
####  LOAD DATA ----------------------------------------------------------------
####
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1) 
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")
bison_dat <- bison_raw %>% 
  dplyr::select(-source) %>%     # drop the source column
  filter(year < 2011) %>%        # only use up to 2010 for now
  left_join(snow_ynp, by="year") # merge in SNOTEL data



####
####  PLOT BISON AND SNOW DATA -------------------------------------------------
####
plot_data <- bison_dat %>%
  dplyr::select(year, count.mean, count.sd, mean_snow_water_equiv_mm) %>%
  dplyr::rename(avg_swe = mean_snow_water_equiv_mm)

bison_growth_data <- bison_dat %>%
  dplyr::select(year, count.mean) %>%
  mutate(id = 1) %>% # constant id to work with ave()
  mutate(growth_rate = ave(count.mean, id, FUN=function(x) c(0, diff(log(x)))))

docolor <- "#278DAF"

bison_plot <- ggplot(plot_data, aes(x = year, y = count.mean))+
  geom_line(color = docolor, alpha = 0.6)+
  geom_point(size=1.5, color = docolor)+
  geom_errorbar(aes(ymin = count.mean-count.sd, ymax = count.mean+count.sd), width=0.5, size=0.5, color = docolor)+
  ylab("Number of bison")+
  xlab("Year")+
  my_theme

bison_growth <- ggplot(bison_growth_data, aes(x = year, y = growth_rate))+
  geom_line(color = docolor, alpha = 0.6)+
  geom_point(size=1.5, color = docolor)+
  ylab("Population growth rate (r)")+
  xlab("Year")+
  my_theme

snow_plot <- ggplot(plot_data, aes(x = year, y = avg_swe))+
  geom_line(color = docolor, alpha = 0.6)+
  geom_point(size=1.5, color = docolor)+
  ylab("Mean SWE (mm)")+
  xlab("Year")+
  my_theme

the_plots <- list(bison_plot, bison_growth, snow_plot)
suppressWarnings( # Ignore warning about 6 NA rows for errorbars where sd not reported
  plot_grid(plotlist = the_plots, labels = "AUTO", ncol = 3)
)
ggsave(filename = "../figures/bison_data_plots.png", height = 3, width = 10, units = "in", dpi = 120)


####
####  JAGS State-Space Model ---------------------------------------------------
####
r_mu_prior <- log(1.11) # lambda = 1.11 in Hobbs et al. 2015
r_sd_prior <- sd(log(rnorm(100000,1.11,0.024))) # sd_lambda = 0.024 in Hobbs et al. 2015
my_model <- "  
  model{

    #### Variance Priors
    sigma_proc ~ dgamma(0.01,0.01)
    tau_proc <- 1/sigma_proc^2
    
    #### Fixed Effects Priors
    r  ~ dnorm(0.1, 1/0.02^2) # intrinsic growth rate, informed prior
    b  ~ dnorm(0,0.0001)      # strength of density dependence (r/K)
    b1 ~ dnorm(0,0.0001)      # effect of snow
    
    #### Initial Conditions
    z[1] ~ dnorm(Nobs[1], tau_obs[1]) # varies around observed abundance at t = 1
    
    #### Process Model
    for(t in 2:npreds){
      mu[t] <- max( 1, log( z[t-1]*exp(r + b*z[t-1] + b1*x[t]) ) )
      z[t] ~ dlnorm(mu[t], tau_proc)
    }
    
    #### Data Model
    for(j in 2:n){
      Nobs[j] ~ dnorm(z[j], tau_obs[j])
    }
  
  }"



####
####  Fit Bison Forecasting Model ----------------------------------------------
####
##  For years without observation error, set to max observed standard deviation
##  TODO: Impute in the model?
na_sds                       <- which(is.na(bison_dat$count.sd)==T)
bison_dat[na_sds,"count.sd"] <- max(bison_dat$count.sd, na.rm=T)

##  Set up SWE forecasts, relative to scaling of observations
swe_mean     <- mean(bison_dat$mean_snow_water_equiv_mm)
swe_sd       <- sd(bison_dat$mean_snow_water_equiv_mm)
forecast_swe <- snow_ynp %>%
  filter(year > max(bison_dat$year)) %>%
  pull(mean_snow_water_equiv_mm)
scl_fut_swe  <- (forecast_swe - swe_mean) / swe_sd
scl_fut_swe  <- c(scl_fut_swe, rep(0,(10 - length(scl_fut_swe)))) # tack on 0s, avg swe, until reach 10 years

##  Prepare data list
mydat <- list(Nobs    = bison_dat$count.mean, # mean counts
              n       = nrow(bison_dat), # number of observations
              tau_obs = 1/bison_dat$count.sd^2, # transform s.d. to precision
              x       = c(as.numeric(scale(bison_dat$mean_snow_water_equiv_mm)),scl_fut_swe), # snow depth, plus forecast years
              npreds  = nrow(bison_dat)+10) # number of total predictions (obs + forecast)

##  Random variables to collect
out_variables <- c("r", "b", "b1", "sigma_proc", "z")

##  Send to JAGS
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=3) 
           update(mc3, n.iter = 10000) 
mc3.out <- coda.samples(model=mc3, 
                        variable.names=out_variables, 
                        n.iter=10000) 
# summary(mc3.out)$stat
# summary(mc3.out)$quantile

## Split output
out          <- list(params=NULL, predict=NULL, model=my_model,data=mydat)
mfit         <- as.matrix(mc3.out,chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
bison_dat[na_sds,"count.sd"] <- NA # set the unobserved std. devs. back to NA for plotting
predictions        <- rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = c(bison_dat$year, (max(bison_dat$year)+1):(max(bison_dat$year)+10)),
                                 observation = c(bison_dat$count.mean,rep(NA,10)),
                                 upper_observation = c(bison_dat$count.mean+bison_dat$count.sd,rep(NA,10)),
                                 lower_observation = c(bison_dat$count.mean-bison_dat$count.sd,rep(NA,10)),
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)

##  Check parameter chains for convergence and mixing
# plot(fitted_model$params)
# gelman.diag(fitted_model$params)
# heidel.diag(fitted_model$params)



####
####  PLOT POSTERIOR DISTRIBUTIONS OF PARAMETERS -------------------------------
####
post_params <- as.data.frame(as.matrix(fitted_model$params))
max_iters   <- nrow(post_params)
post_params <- post_params %>%
  mutate(iteration = 1:max_iters) %>%
  gather(key = parameter, value = estimate, -iteration) %>%
  mutate(prior = c(rnorm(max_iters,0,1000), # b prior
                   rnorm(max_iters,0,1000), # b1 prior
                   rnorm(max_iters,r_mu_prior,r_sd_prior), # r prior
                   runif(max_iters,0,10))) # sd prior

prior_col <- "#CF4C26"
ggplot(post_params, aes(x = estimate, y = ..density..))+
  geom_histogram(fill = docolor, color = "white", bins = 20)+
  geom_line(data = filter(post_params, parameter == "r"),
            aes(x = prior), 
            stat = "density", 
            color = prior_col)+
  facet_wrap(~parameter, scales = "free", ncol = 4)+
  ylab("Posterior density")+
  xlab("Parameter estimate")+
  my_theme
ggsave(filename = "../figures/bison_post_params.png", 
       height = 3, 
       width = 10, 
       units = "in", 
       dpi = 120)



####
####  PLOT DATA AND POSTERIOR PREDICTIONS --------------------------------------
####
pred_color <- "#CF4C26"
obs_color  <- "#278DAF"
calibration_plot <- ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction),
              fill=pred_color, 
              color=NA, 
              alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), 
                width=0.5, 
                color=obs_color, 
                size=0.2)+
  geom_point(aes(y=observation), color=obs_color, size=0.5)+
  geom_vline(aes(xintercept=max(bison_dat$year)), linetype=2,color="grey55")+
  ylab("Number of bison")+
  xlab("Year")+
  my_theme



####
####  Partition Forecast Uncertainty -------------------------------------------
####
##  Function for the ecological process (population growth)
iterate_process <- function(Nnow, xnow, r, b, b1, sd_proc) { 
  Ntmp <- Nnow*exp(r + b*Nnow + b1*xnow)
  Ntmp[which(Ntmp < 1)] <- 1 # catch "bad" values
  N    <- rlnorm(length(Nnow), log(Ntmp), sd_proc)
}


##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values and no process error.
forecast_steps <- 10
num_iters      <- 1000
x              <- scl_fut_swe
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
param_summary  <- summary(fitted_model$params)$quantile
b              <- param_summary[1,3]
b1             <- param_summary[2,3]
r              <- param_summary[3,3]
sd_proc        <- param_summary[4,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r = r, b = b, b1 = b1, sd_proc = 0)
  forecasts[,t] <- z
}
varI <- apply(forecasts,2,var)


##  Initial conditions and parameter uncertainty
x              <- scl_fut_swe
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
b              <- params[sample_params,1]
b1             <- params[sample_params,2]
r              <- params[sample_params,3]
sd_proc        <- param_summary[3,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r = r, b = b, b1 = b1, sd_proc = 0)
  forecasts[,t] <- z
}
varIP <- apply(forecasts,2,var)


##  Initial conditions, parameter, and process uncertainty
x              <- scl_fut_swe
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
b              <- params[sample_params,1]
b1             <- params[sample_params,2]
r              <- params[sample_params,3]
sd_proc        <- params[sample_params,4]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r = r, b = b, b1 = b1, sd_proc = sd_proc)
  forecasts[,t] <- z
}
varIPE <- apply(forecasts,2,var)


V.pred.sim.rel <- apply(rbind(varIPE,varIP,varI),2,function(x) {x/max(x)})



####
####  Plot the proportion of uncertainty by partition
####
var_rel_preds <- as.data.frame(t(V.pred.sim.rel*100))
var_rel_preds$x <- 1:nrow(var_rel_preds)
my_cols <- c("#0A4D5B", "#139AB8", "#39B181")
variance_plot <- ggplot(data=var_rel_preds, aes(x=x))+
  geom_ribbon(aes(ymin=0, ymax=varI), fill=my_cols[1])+
  geom_ribbon(aes(ymin=varI, ymax=varIP), fill=my_cols[2])+
  geom_ribbon(aes(ymin=varIP, ymax=varIPE), fill=my_cols[3])+
  ylab("Percent of uncertainty")+
  xlab("Forecast steps")+
  scale_x_continuous(breaks=seq(2,forecast_steps,by=2), 
                     labels=paste(seq(2,forecast_steps,by=2), "yrs"))+
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))+
  my_theme



####
####  COMBINE PLOTS AND SAVE ----
####
suppressWarnings( # ignore warning about missing values, we know they are missing
  plot_grid(calibration_plot, variance_plot, nrow = 2, labels = "AUTO")
)
ggsave(filename = "../figures/bison_combined.png", 
       width = 4, 
       height = 6, 
       units = "in", 
       dpi =120)

