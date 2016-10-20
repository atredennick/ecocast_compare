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


##  Load snotel data, aggregate to get average snow depth by year
snotel_raw <- read.csv("../data/west_yellowstone_snotel.csv", na.strings = "", stringsAsFactors = FALSE)
snotel_raw <- snotel_raw[2:nrow(snotel_raw), ] # removes id row
snotel_dat <- melt(snotel_raw, id.vars = c("Water.Year","Day"))
snotel_dat$value <- as.numeric(snotel_dat$value)
snotel_means <- ddply(snotel_dat, .(Water.Year), summarise,
                      avg_snow = mean(value))

## Merge observation and climate data
bison_climate_dat <- merge(bison_dat, snotel_means, by.x = "year", by.y="Water.Year")

##  Retrieve climate data from NOAA for forecasts




####
####  Fit forecasting GLM
####
my_model <- list(obs="count.mean", fixed=NULL, random=NULL, n.iter=10000)
fitted_model <- fit_dlm_lnorm(model=my_model, data=bison_climate_dat)
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

# params <- as.data.frame(rbind(fitted_model$params[[1]],
#                               fitted_model$params[[2]],
#                               fitted_model$params[[3]]))

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



####
####  Partition Forecast Uncertainty
####
nens <- 1
nsteps <- 10
meta.ds <- list()
meta.ds[[1]] = 1:nens
meta.ds[[2]] = max(bison_climate_dat$year):(max(bison_climate_dat$year)+nsteps)
meta.ds[[3]] = c("Intercept")
newdata = array(1,dim = c(nens,nsteps,1), dimnames = "Intercept")

## just initial condition uncertainty
FE_pred.I <- predict_dlm_lnorm(fit=fitted_model, newdata = newdata, n.iter=500,include="I", steps=nsteps, start.time = max(bison_climate_dat$year))

## initial conditions + parameters
FE_pred.IP <- predict_dlm_lnorm(fitted_model,n.iter=n.iter.fx,include=c("I","P"))

plot_ss(meta.ds[[2]],FE_pred.I,ylab="NEE",xlab="time")
plot_ss(meta.ds[[2]],FE_pred.IP,ylab="NEE",xlab="time")


