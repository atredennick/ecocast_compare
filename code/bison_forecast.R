##  R script to fit a population growth model for YNP bison,
##  forecast 10 new years, and partition the forecast variance.
##
##  Based on "Ecological Forecasting" by M. Dietze (2017)
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: October 19, 2016
##


rm(list=ls(all.names = TRUE))


####
####  Load libraries
####
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(reshape2)
library(plyr)
library(rjags)
library(coda)
# library(devtools)
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR)



####
####  Set Up My Plotting Theme -------------------------------------------------
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
        axis.line.y        = element_line(color="black"))



####
####  Load Data ----------------------------------------------------------------
####
bison_raw <- read.csv("../data/YNP_bison_population_size.csv", row.names = 1)
bison_dat <- bison_raw[which(!is.na(bison_raw$count.sd)),2:ncol(bison_raw)]



####
####  JAGS State-Space Model ---------------------------------------------------
####
my_model <- "  
  model{

    #### Variance Priors
    tau_proc ~ dgamma(0.0001, 0.0001)
    sigma_proc <- 1/sqrt(tau_proc)
    
    #### Fixed Effects Priors
    r ~ dunif(0,5)
    K ~ dunif(1,10000)
    
    #### Initial Conditions
    N0 ~ dunif(1,1000)
    Nmed[1] <- log(max(1, N0 + r * N0 * (1 - N0 / K)))
    N[1] ~ dlnorm(Nmed[1], tau_proc)
    
    #### Process Model
    for(t in 2:npreds){
      Nmed[t] <- log(max(1, N[t-1] + r * N[t-1] * (1 - N[t-1] / K)))
      N[t] ~ dlnorm(Nmed[t], tau_proc) 
    }
    
    #### Data Model
    ##  SD observations
    for(t in 1:n){
      var_obs[t] <- sd_obs[t]*sd_obs[t]
      shape[t] <- N[t]*N[t]/var_obs[t]
      rate[t] <- N[t]/var_obs[t]
      lambda[t] ~ dgamma(shape[t], rate[t])
      Nobs[t] ~ dpois(lambda[t])
    }

  }"



####
####  Fit Bison Forecasting Model ----------------------------------------------
####

##  Prepare data list
mydat         <- list(Nobs = round(bison_dat$count.mean), 
                      n = nrow(bison_dat),
                      sd_obs = bison_dat$count.sd,
                      npreds = nrow(bison_dat)+10)
out_variables <- c("r","K","sigma_proc","N")

##  Send to JAGS
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=3)
           update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=50000)

## Split output
out          <- list(params=NULL, predict=NULL, model=my_model,data=mydat)
mfit         <- as.matrix(mc3.out,chains=TRUE)
pred.cols    <- union(grep("N[",colnames(mfit),fixed=TRUE),grep("Nmed[",colnames(mfit),fixed=TRUE))
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
####  Plot the calibration data and predictions
####
pred_color <- "#CF4C26"
obs_color  <- "#0A9AB8"
calibration_plot <- ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), width=0.5, color=obs_color, size=0.2)+
  geom_point(aes(y=observation), color=obs_color, size=0.5)+
  geom_vline(aes(xintercept=max(bison_dat$year)), linetype=2,color="grey55")+
  ylab("Number of bison")+
  xlab("Year")+
  my_theme
# ggsave(filename = "../figures/bison_calibration.png", width = 4, height = 3, units = "in", dpi=120)



####
####  Partition Forecast Uncertainty -------------------------------------------
####
##  Function for the ecological process (population growth)
iterate_process <- function(Nnow, r, K, sd_proc) { 
  Ntmp <- Nnow + r * Nnow * (1 - Nnow / K)
  Ntmp[which(Ntmp < 1)] <- 1
  N    <- rlnorm(length(Nnow), log(Ntmp), sd_proc)
}


##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values and no process error.
forecast_steps <- 10
num_iters      <- 1000
x              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
param_summary  <- summary(fitted_model$params)$quantile
K              <- param_summary[1,3]
r              <- param_summary[2,3]
sd_proc        <- param_summary[3,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, r = r, K = K, sd_proc = 0)
  forecasts[,t] <- x
}
varI <- apply(forecasts,2,var)


##  Initial conditions and parameter uncertainty
x              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
K              <- params[sample_params,1]
r              <- params[sample_params,2]
sd_proc        <- param_summary[3,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, r = r, K = K, sd_proc = 0)
  forecasts[,t] <- x
}
varIP <- apply(forecasts,2,var)


##  Initial conditions, parameter, and process uncertainty
x              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
K              <- params[sample_params,1]
r              <- params[sample_params,2]
sd_proc       <- params[sample_params,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, r = r, K = K, sd_proc = sd_proc)
  forecasts[,t] <- x
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
  scale_x_continuous(breaks=seq(2,forecast_steps,by=2), labels=paste(seq(2,forecast_steps,by=2), "yrs"))+
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))+
  my_theme




png("../figures/bison_combined.png", width = 4, height = 6, units = "in", res = 120)
out_plot <- grid.arrange(calibration_plot, variance_plot, ncol=1)
dev.off()
