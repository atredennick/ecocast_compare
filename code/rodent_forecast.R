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
# library(devtools)
# install_github("atredennick/ecoforecastR") # get MY latest version
# install_github("weecology/PortalDataSummaries") # get Weecology summarizing package
library(tidyverse)
library(dplyr)
library(ggthemes)
library(rjags)
library(coda)
library(ecoforecastR)
library(PortalDataSummaries)



####
####  READ IN AND FORMAT DATA ----
####
rodents <- abundance(path = "../data", 
                     incomplete = FALSE,
                     length = "longterm", 
                     shape="flat",
                     level = "Plot")

moons <- read.csv("../data/PortalData/Rodents/moon_dates.csv") %>%
  select(Period, CensusDate)

plot_info <- read.csv("../data/PortalData/SiteandMethods/Portal_plot_treatments.csv")

rodents <- rodents %>%
  left_join(moons, by = c("period" = "Period")) %>%
  rename(census_date = CensusDate) %>%
  separate(census_date, c("year","month","day"))

pp_data <- rodents %>%
  filter(species=="PP")

pp_agg <- pp_data %>%
  group_by(year,plot) %>%
  summarise(avg_abundance = mean(abundance),
            sdv_abundance = sd(abundance)) %>%
  ungroup()

ggplot(pp_agg, aes(x=as.numeric(year), y=avg_abundance, color=as.factor(plot))) +
  geom_line()+
  geom_point()


