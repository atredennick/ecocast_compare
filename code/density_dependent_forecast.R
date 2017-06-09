##  density_dependent_forecast.R: script to simulate time series of population
##  abundances across a range of density dependence. The time series are then
##  used to fit state-space forecasting models. Forecasts of population
##  abundance are made for ten years, and the forecast uncertainty is
##  partitioned using the methods from Dietze (2017).
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: June 6, 2017


##  Clear the workspace
rm(list=ls(all.names = TRUE))

##  Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # only for RStudio



####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(dplyr)
library(ecoforecastR)



####
####  FUNCTION FOR GOMPERTZ POPULATION GROWTH ----
####
# Function to simulate population time series using Gompertz population
# growth equation in log space. See Knape and de Valpine (2012, Ecol. Letts.)
# for details.
# In a nutshell, model is:
#       x(t+1) = a + c*x(t) + error(t)
#       y(t) = x(t) + eta(t)
# where x is log abundance, a is an intercept, c is the strength of density
# dependence, and error(t) is distributed as N(0,tau).

sim_gompertz <- function(generations=100, a=0, c=1, tau=sqrt(0.5), sigma=sqrt(0.5), Nstart=1){
  logN <- numeric(generations)
  logN[1] <- log(Nstart)
  for(t in 2:generations){
    logN[t] <- a + c*logN[t-1] + rnorm(1,0,tau)
  }
  logY <- logN + rnorm(generations,0,sigma)
  return(logY[21:100])
}



####
####  SIMULATE A BUNCH OF TIME SERIES ACROSS GRADIENT OF DD (c) ----
####
c_vector <- seq(-0.99,0.99,0.1)
pop_list <- list()
for(do_c in c_vector){
  pop_list[[as.character(do_c)]] <- sim_gompertz(c=do_c)
}



####
####  FUNCTION TO FIT STATE-SPACE GOMPERTZ MODEL ----
####
fit_ss_gompertz <- function(iters=5000, chains=2, data_list){
  gompertz_model <- "
  model{
    ### Process model
    z[1] ~ dunif(0,10)
    for(t in 2:pred_times){
      mu[t] = a + c*z[t-1]
      z[t] ~ dnorm(mu[t], tau)
    }
    
    ### Observation model
    for(i in 1:nobs){
      y[i] ~ dnorm(z[i], sigma)
    }
    
    ### Parameter model
    a ~ dnorm(0,0.001)
    c ~ dnorm(0,0.001) T(-0.99,0.99)
    tau ~ dgamma(0.001,0.001)
    sigma ~ dgamma(0.001,0.001)
    
    ### Generated quantity
    c_strength = 1-c
  }
  "
  
  model <- jags.model(textConnection(gompertz_model), data=data_list, n.adapt = 1000)
  update(model, n.iter=iters)
  output <- coda.samples(model=model,variable.names=c("c"), n.iter=iters)
  return(output)
}



####
####  FIT THE MODELS AND SAVE ----
####
pdf("test.pdf")
par(mfrow=c(4,5))
for(i in 1:20){
  dopop <- i
  tmp_data <- pop_list[[dopop]]
  data_list <- list(y = tmp_data,
                    nobs = length(tmp_data),
                    pred_times = length(tmp_data)+10)
  mcmc <- fit_ss_gompertz(data_list = data_list)
  # plot(tmp_data,type="l")
  plot(density(mcmc[[1]]), xlim=c(-1,1),main="",xlab=expression(hat(c)))
  abline(v = c_vector[dopop],col="red")
}
dev.off()




