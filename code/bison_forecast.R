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
library(rnoaa)
library(ncdf4)
# install_github("atredennick/ecoforecastR", force=TRUE) # get latest version
library(ecoforecastR)



####
####  Load data, aggregate to yearly values
####
bison_raw <- read.csv("../data/YNP_bison_population_size.csv", row.names = 1)
bison_dat <- bison_raw[,2:ncol(bison_raw)]


##  Retrieve climate data from NOAA
# lat <- 44.4280
# lon <- 110.5885
# out <- ncdc(datasetid='GHCND', stationid='GHCND:USW00014895', datatypeid='PRCP', startdate = '2010-05-01', enddate = '2010-10-31', limit=500)
# 
# ##  Get forecasted snow
# gefs_variables()
# snow_depth = gefs("Snow_depth_surface_ens",lat,lon,raw=TRUE)

# clim_dat <- read.csv("../data/idaho_climate.csv")

## Merge observation and climate data
# bison_climate_dat <- merge(bison_dat, clim_dat)


####
####  Fit forecasting GLM
####
my_model <- list(obs="count.mean", fixed=NULL, random=NULL, n.iter=5000)
fitted_model <- fit_dlm_lnorm(model=my_model, data=bison_dat)
predictions <- rbind(fitted_model$predict[[1]],
                     fitted_model$predict[[2]],
                     fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})

prediction_df <- data.frame(year = bison_dat$year,
                            observation = bison_dat$count.mean,
                            median_prediction = exp(median_predictions),
                            upper_prediction = exp(upper_predictions),
                            lower_prediction = exp(lower_predictions))



####
####  Plot the calibration data and predictions
####
pred_color <- "dodgerblue"
ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_point(aes(y=observation))+
  ylab("Number of bison")+
  xlab("Year")+
  theme_few()
ggsave(filename = "../figures/bison_calibration.png", width = 4, height = 3, units = "in", dpi=120)



