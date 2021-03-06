##  R script to fit a simple GLM for sagebrush percent cover,
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
library(devtools)
# install_github("atredennick/ecoforecastR", force=TRUE) # get latest version
library(ecoforecastR)



####
####  Load data, aggregate to yearly values
####
sage_raw <- read.csv("../data/ARTR_quadratCover.csv")
sage_raw$year <- sage_raw$year+1900 # makes a calendar year
sage_dat <- ddply(sage_raw, .(year), summarise,
                  tot_cover = round(mean(totCover)/100)) # converts to percent cover, integer

clim_dat <- read.csv("../data/idaho_climate.csv")

## Merge observation and climate data
sage_climate_dat <- merge(sage_dat, clim_dat)


####
####  Fit forecasting GLM
####
my_model <- list(obs="tot_cover", fixed="~ppt1+TmeanSpr1", random=NULL, n.iter=5000)
fitted_model <- fit_dlm_pois(model=my_model, data=sage_climate_dat)
predictions <- rbind(fitted_model$predict[[1]],
                     fitted_model$predict[[2]],
                     fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})

prediction_df <- data.frame(year = sage_climate_dat$year,
                            observation = sage_climate_dat$tot_cover,
                            median_prediction = median_predictions,
                            upper_prediction = upper_predictions,
                            lower_prediction = lower_predictions)



####
####  Plot the calibration data and predictions
####
pred_color <- "dodgerblue"
ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_point(aes(y=observation))+
  ylab("Cover of sagebrush (%)")+
  xlab("Year")+
  theme_few()
ggsave(filename = "../figures/sage_calibration.png", width = 4, height = 3, units = "in", dpi=120)



####
####  Partition Forecast Uncertainty
####
nens <- 1
nsteps <- 10
meta.ds <- list()
meta.ds[[1]] = 1:nens
meta.ds[[2]] = length(sage_climate_dat$year):(length(sage_climate_dat$year)+nsteps-1)
meta.ds[[3]] = c("Intercept")
newdata = array(1,dim = c(nens,nsteps,1), dimnames = "Intercept")

## just initial condition uncertainty
FE_pred.I <- predict_dlm_pois(fit=fitted_model, newdata = newdata, n.iter=500,include="I", steps=nsteps, start.time = NULL)





