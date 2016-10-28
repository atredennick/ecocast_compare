##  R script to fit a population growth model for YNP,
##  forecast 10 new years, and partition the forecast variance.
##
##  Based on Dietze et al. (forthcoming)
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
library(reshape2)
library(plyr)
library(rjags)
library(coda)
# library(devtools)
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR)



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
    for(t in 2:n){
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
####  Fit Bison Forecasting Model
####

##  Prepare data list
mydat         <- list(Nobs = round(bison_dat$count.mean), 
                      n = nrow(bison_dat),
                      sd_obs = bison_dat$count.sd)
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
prediction_df      <- data.frame(year = bison_dat$year,
                                 observation = bison_dat$count.mean,
                                 upper_observation = bison_dat$count.mean+bison_dat$count.sd,
                                 lower_observation = bison_dat$count.mean-bison_dat$count.sd,
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)

##  Check parameter chains for convergence and mixing
plot(fitted_model$params)
gelman.diag(fitted_model$params)
heidel.diag(fitted_model$params)



####
####  Plot the calibration data and predictions
####
pred_color <- "dodgerblue"
ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), width=0.5, color="grey45")+
  geom_point(aes(y=observation), size=0.5)+
  ylab("Number of bison")+
  xlab("Year")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dotted", color="grey65"))
ggsave(filename = "../figures/bison_calibration.png", width = 4, height = 3, units = "in", dpi=120)


# 
# ####
# ####  Partition Forecast Uncertainty
# ####
# nens <- 1
# nsteps <- 30
# meta.ds <- list()
# meta.ds[[1]] = 1:nens
# meta.ds[[2]] = length(bison_climate_dat$year):(length(bison_climate_dat$year)+nsteps-1)
# meta.ds[[3]] = c("Intercept")
# newdata = array(1,dim = c(nens,nsteps,1), dimnames = "Intercept")
# 
# ## just initial condition uncertainty
# FE_pred.I <- predict_dlm_lnorm(fit=fitted_model, newdata = newdata, n.iter=500,include="I", steps=nsteps, start.time = NULL)
# 
# ## initial conditions + parameters
# FE_pred.IP <- predict_dlm_lnorm(fitted_model,newdata = newdata, n.iter=500,include=c("I","P"), steps=nsteps, start.time = NULL)
# 
# ## full uncertainty
# FE_pred.IPE <- predict_dlm_lnorm(fitted_model,newdata = newdata, n.iter=500,include=c("I","P","E"), steps=nsteps, start.time = NULL)
# 
# 
# 
# ## FULL
# plot_ss(meta.ds[[2]],FE_pred.IPE,ylab="NEE",xlab="Day of Year")
# varIPE <- apply(as.matrix(FE_pred.IPE$predict),2,var)
# 
# ## IP
# ciIP <- apply(as.matrix(FE_pred.IP$predict),2,quantile,c(0.025,0.5,0.975))
# ciEnvelope(meta.ds[[2]],ciIP[1,],ciIP[3,],col="lightGreen")
# varIP <- apply(as.matrix(FE_pred.IP$predict),2,var)
# 
# ## I
# ciI <- apply(as.matrix(FE_pred.I$predict),2,quantile,c(0.025,0.5,0.975))
# ciEnvelope(meta.ds[[2]],ciI[1,],ciI[3,],col="violet")
# lines(meta.ds[[2]],ciI[2,],col="darkGreen",lwd=2)
# varI <- apply(as.matrix(FE_pred.I$predict),2,var)
# 
# V.pred.sim.rel <- apply(rbind(varIPE,varIP,varI),2,function(x) {x/max(x)})
# 
# 
# 
# ####
# ####  Plot the proportion of uncertainty by partition
# ####
# var_rel_preds <- as.data.frame(t(V.pred.sim.rel))
# var_rel_preds$x <- 1:nrow(var_rel_preds)
# ggplot(data=var_rel_preds, aes(x=x))+
#   geom_ribbon(aes(ymin=0, ymax=varI), fill="grey35")+
#   geom_ribbon(aes(ymin=varI, ymax=varIP), fill="coral")+
#   geom_ribbon(aes(ymin=varIP, ymax=varIPE), fill="skyblue")+
#   ylab("Proportion of uncertainty")+
#   xlab("Forecast steps")+
#   # scale_x_continuous(breaks=seq(1,nsteps,by=1))+
#   theme_few()
# ggsave(filename = "../figures/bison_forecast_uncertainty.png", width = 4, height = 3, units = "in", dpi=120)
