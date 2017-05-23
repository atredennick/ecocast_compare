##  R script to simulate population counts through time for populations
##  with different growth rates. Then I fit a state-space model to the data,
##  forecast 10 new years, and partition the forecast variance. The idea is
##  to generate predictions for testable hypotheses linking life history to
##  forecast uncertainty.
##
##  Partition approach based on "Ecological Forecasting" by M. Dietze (2017)
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: May 23, 2017
##

##  Clear the workspace
rm(list=ls(all.names = TRUE))

## Set working dir to source file location (only works in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 


####
####  Load libraries
####
# library(devtools)
# install_github("atredennick/ecoforecastR") # get MY latest version
library(tidyverse)
library(dplyr)
library(ggthemes)
library(gridExtra)
library(rjags)
library(coda)
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
####  LOGISTIC POPULATION GROWTH FUNCTION (DISCRETE TIME) ----
####
sim_logistic <- function(timesteps=40,num_observations=4,N0=500,r=0.18,
                         K=2000,sd_proc_log=0.2,sd_obs=0.2){
  Nobs <- matrix(data=NA,nrow = timesteps,ncol = num_observations)
  N <- numeric(timesteps)
  N[1] <- Nobs[1,] <- rnorm(num_observations,N0,N0*sd_obs)
  for(t in 2:timesteps){
    N[t] <- N[t-1] + r*N[t-1] * (1 - (N[t-1]/K))
    if(N[t] < 1) { N[t] <- 1 }
    N[t] <- rlnorm(1,log(N[t]), sd_proc_log)
    sdnow <- sd_obs*N[t]
    Nobs[t,] <- rnorm(num_observations,N[t],sdnow)
  }
  Nobs_out <- as.data.frame(Nobs) %>%
    mutate(year = c(1:timesteps)) %>%
    gather(observer,abundance,-year) %>%
    mutate(true_abundance = rep(N,num_observations))
  return(Nobs_out)
}



####
####  SIMULATE LOGISTIC MODEL AND SAVE OUTPUT AS "DATA" ----
####
nsim <- sim_logistic(r=0.18)

abund_data <- nsim %>%
  group_by(year) %>%
  summarise(mean_abundance = mean(abundance),
            sdev_abundance = sd(abundance)) %>%
  ungroup()

ggplot(abund_data, aes(x=year, y=mean_abundance))+
  geom_line()+
  geom_errorbar(aes(ymin=mean_abundance-sdev_abundance, ymax=mean_abundance+sdev_abundance))+
  geom_point()+
  my_theme



####
####  JAGS State-Space Model ----
####
my_model <- "  
  model{

    #### Variance Priors
    tau_proc ~ dgamma(0.0001, 0.0001)
    sigma_proc <- 1/sqrt(tau_proc)
    
    #### Fixed Effects Priors
    r ~ dunif(0,2)
    K ~ dunif(100,10000)
    
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
####  FIT FORECASTING MODEL ----
####
##  Prepare data list
mydat         <- list(Nobs = round(abund_data$mean_abundance), 
                      n = nrow(abund_data),
                      sd_obs = abund_data$sdev_abundance,
                      npreds = nrow(abund_data)+10)
out_variables <- c("r","K","sigma_proc","N")

##  Send to JAGS
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=3)
           update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=10000)

## Split output
out          <- list(params=NULL, predict=NULL, model=my_model,data=mydat)
mfit         <- as.matrix(mc3.out,chains=TRUE)
pred.cols    <- union(grep("N[",colnames(mfit),fixed=TRUE),grep("Nmed[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

# par(mfrow=c(1,3))
# plot(density(mfit[,"r"]),xlab=expression(italic("r")), main="")
# lines(density(runif(nrow(mfit),0,2), adjust=2),lty=2)
# plot(density(mfit[,"K"]),xlab=expression(italic("K")), main="")
# lines(density(rnorm(nrow(mfit),2000,1/sqrt(0.00001)), adjust=2),lty=2)
# plot(density(mfit[,"sigma_proc"]),xlab=expression(sigma[p]), main="")
# lines(density(1/sqrt(rgamma(nrow(mfit),0.0001,0.0001)), adjust=2),lty=2)

## Collate predictions
predictions        <- rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = c(abund_data$year, (max(abund_data$year)+1):(max(abund_data$year)+10)),
                                 observation = c(abund_data$mean_abundance,rep(NA,10)),
                                 upper_observation = c(abund_data$mean_abundance+abund_data$sdev_abundance,rep(NA,10)),
                                 lower_observation = c(abund_data$mean_abundance-abund_data$sdev_abundance,rep(NA,10)),
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)



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
x              <- sample(predictions[,nrow(abund_data)], num_iters, replace = TRUE)
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
x              <- sample(predictions[,nrow(abund_data)], num_iters, replace = TRUE)
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
x              <- sample(predictions[,nrow(abund_data)], num_iters, replace = TRUE)
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
####  Plot the calibration data and predictions
####
pred_color <- "#CF4C26"
obs_color  <- "#0A9AB8"
calibration_plot <- ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), width=0.5, color=obs_color, size=0.2)+
  geom_point(aes(y=observation), color=obs_color, size=0.5)+
  geom_vline(aes(xintercept=max(abund_data$year)), linetype=2,color="grey55")+
  # annotate(geom="text", x=1988, y=6100,
  #          label=paste0("population growth rate = ",round(exp(param_summary[2,3]),2)),
  #          size=4, family = "Arial Narrow")+
  ylab("Number of bison")+
  xlab("Year")+
  my_theme



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




# png("../figures/bison_combined.png", width = 4, height = 6, units = "in", res = 120)
out_plot <- grid.arrange(calibration_plot, variance_plot, ncol=1)
# dev.off()
