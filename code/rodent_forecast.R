##  R script to fit a population growth model for portal rodent(s),
##  forecast 10 new years, and partition the forecast variance.
##
##  Based on "Ecological Forecasting" by M. Dietze (2017)
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: May 19, 2017
##

##  Clear the workspace
rm(list=ls(all.names = TRUE))


####
####  LOAD LIBRARIES ----
####
library(RCurl)
library(tidyverse)
library(dplyr)
library(ggthemes)
library(rjags)
library(coda)
# library(devtools)
# install_github("atredennick/ecoforecastR") # get MY latest version
library(ecoforecastR)



####
####  READ IN AND FORMAT DATA ----
####
github_url   <- "https://raw.githubusercontent.com/weecology/PortalData/master/"
rodent_url   <- paste0(github_url,"Rodents/Portal_rodent.csv")
plots_url    <- paste0(github_url,"SiteandMethods/Portal_plots.csv")
rodents      <- read.csv(text=getURL(rodent_url))
portal_plots <- read.csv(text=getURL(plots_url))


rodents_agg <- rodents %>%
  mutate(num_inds = 1) %>%
  group_by(yr,plot,species) %>%
  summarise(total_inds = sum(num_inds)) %>%
  ungroup()
