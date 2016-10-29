##  R script to fit a simple GLM for sagebrush percent cover,
##  forecast 10 new years, and partition the forecast variance.
##
##  Based on Dietze et al. (forthcoming)
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: October 19, 2016
##


rm(list=ls(all.names = TRUE))
VERBOSE <- FALSE


####
####  Load libraries
####
library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(rjags)
library(coda)
library(devtools)
# install_github("atredennick/ecoforecastR", force=TRUE) # get latest version
library(ecoforecastR)



####
####  Load Data, Aggregate to Yearly Values ------------------------------------
####
sage_raw      <- read.csv("../data/ARTR_quadratCover.csv")
sage_raw$year <- sage_raw$year+1900 # makes a calendar year

##  Average by pasture, then mean and sd over pastures
sage_pasture <

# Division by 100 converts to percent cover
sage_dat <- ddply(sage_raw, .(year), summarise,
                  tot_cover = round(mean(totCover)/100),
                  sd_cover = sd(totCover)/100) 

clim_dat <- read.csv("../data/idaho_climate.csv")

## Merge observation and climate data
sage_climate_dat <- merge(sage_dat, clim_dat)


####
####  JAGS State-Space Model ---------------------------------------------------
####
my_model <- "  
  model{
  
    #### Variance Priors
    #tau_proc ~ dgamma(0.0001, 0.0001)
    #sigma_proc <- 1/sqrt(tau_proc)
    shape_p ~ dunif(0, 100)
    
    #### Fixed Effects Priors
    b0 ~ dnorm(0, 0.001)
    b1 ~ dnorm(0, 0.001)
    b2 ~ dnorm(0, 0.001)
    
    #### Initial Conditions
    N0      ~ dunif(0,10)
    Nmed[1] <- b0 + b1*N0 + b2*x[1]
    N[1]    ~ dgamma(shape_p, shape_p / exp(Nmed[1])) 
    
    #### Process Model
    for(t in 2:n){
      Nmed[t] <- b0 + b1*N[t-1] + b2*x[t]
      N[t]    ~ dgamma(shape_p, shape_p / exp(Nmed[t])) 
    }
    
    #### Data Model
    for(t in 1:n){
      var_obs[t] <- sd_obs[t]*sd_obs[t]
      shape[t]   <- N[t]*N[t]/var_obs[t]
      rate[t]    <- N[t]/var_obs[t]
      lambda[t]  ~ dgamma(shape[t], rate[t])
      Nobs[t]    ~ dpois(lambda[t])
    }
  
  }"



####
####  Fit Sagebrush Forecasting Model ------------------------------------------
####
##  Prepare data list
mydat         <- list(Nobs = sage_climate_dat$tot_cover, 
                      n = nrow(sage_climate_dat),
                      sd_obs = sage_climate_dat$sd_cover,
                      x = sage_climate_dat$ppt1)
out_variables <- c("b0","b1","b2","N")

##  Send to JAGS
mc3     <- jags.model(file=textConnection(my_model), data=mydat, n.chains=3)
           update(mc3, n.iter = 10000)
mc3.out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=20000)

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
prediction_df      <- data.frame(year = sage_climate_dat$year,
                                 observation = sage_climate_dat$tot_cover,
                                 upper_observation = sage_climate_dat$tot_cover+sage_climate_dat$sd_cover,
                                 lower_observation = sage_climate_dat$tot_cover-sage_climate_dat$sd_cover,
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)

##  Check parameter chains for convergence and mixing
if(VERBOSE==TRUE){
  plot(fitted_model$params)
  gelman.diag(fitted_model$params)
  heidel.diag(fitted_model$params)
}




####
####  Plot the Calibration Data and Predictions --------------------------------
####
pred_color <- "dodgerblue"
ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), width=0.5, color="grey45")+
  geom_point(aes(y=observation), size=0.5)+
  scale_y_continuous(breaks = seq(0,25,5))+
  ylab("Sagebrush cover (%)")+
  xlab("Year")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dotted", color="grey65"))
ggsave(filename = "../figures/sagebrush_calibration.png", width = 4, height = 3, units = "in", dpi=120)

