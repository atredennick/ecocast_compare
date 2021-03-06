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
####  Load Data, Aggregate to Yearly Values ------------------------------------
####
counts          <- read.csv("../data/idaho_annuals_counts_v3.csv")
bromus_counts   <- subset(counts, species=="Bromus tectorum")
bromus_dat      <- ddply(bromus_counts, .(year), summarise,
                         avg_count = mean(count)) 
bromus_dat$year <- bromus_dat$year+1900
clim_dat        <- read.csv("../data/idaho_climate.csv")

## Merge observation and climate data
bromus_climate_dat <- merge(bromus_dat, clim_dat)



####
####  JAGS State-Space Model ---------------------------------------------------
####
my_model <- "  
  model{

    #### Variance Priors
    tau_proc ~ dgamma(0.0001, 0.0001)
    sigma_proc <- 1/sqrt(tau_proc)
    tau_obs ~ dgamma(0.0001, 0.0001)
    sigma_obs <- 1/sqrt(tau_obs)

    #### Fixed Effects Priors
    alpha ~ dnorm(0,0.001)
    beta ~ dnorm(0,0.001)
    betax ~ dnorm(0,0.001)
    
    #### Initial Conditions
    N0 ~ dunif(1,100)
    mu[1] <- alpha + beta*N0 #+ betax*x[1]
    N[1] ~ dnorm(mu[1], tau_proc)

    #### Process Model
    for(t in 2:npreds){
      mu[t] <- alpha + beta*N[t-1] #+ betax*x[t]
      N[t] ~ dnorm(mu[t], tau_proc)  
    }
    
    #### Data Model
    for(t in 1:n){
      Nobs[t] ~ dlnorm(N[t], tau_obs)
    }

  }"



####
####  Fit Sagebrush Forecasting Model ------------------------------------------
####
##  Prepare data list
mydat         <- list(Nobs = bromus_climate_dat$avg_count, 
                      n = nrow(bromus_climate_dat),
                      x = c(bromus_climate_dat$TmeanSpr1, rep(mean(bromus_climate_dat$TmeanSpr1),10)),
                      npreds = nrow(bromus_climate_dat)+10)
out_variables <- c("alpha","beta","sigma_proc","N")

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
predictions        <- exp(rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]]))
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = c(bromus_climate_dat$year, (max(bromus_climate_dat$year)+1):(max(bromus_climate_dat$year)+10)),
                                 observation = c(bromus_climate_dat$avg_count,rep(NA,10)),
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
pred_color <- "#CF4C26"
obs_color  <- "#0A9AB8"
ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_point(aes(y=observation), size=0.5, color=obs_color)+
  ylab("Density of cheatgrass (N / sq. meter)")+
  xlab("Year")+
  my_theme



####
####  Partition Forecast Uncertainty -------------------------------------------
####
##  Function for the ecological process (population growth)
iterate_process <- function(Nnow, alpha, beta, sd_proc) { 
  Ntmp <- alpha+beta*Nnow
  N    <- rnorm(length(Nnow), Ntmp, sd_proc)
}


##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values and no process error.
forecast_steps <- 10
num_iters      <- 1000
x              <- sample(predictions[,nrow(bromus_dat)], num_iters, replace = TRUE)
param_summary  <- summary(fitted_model$params)$quantile
alpha              <- param_summary[1,3]
beta              <- param_summary[2,3]
sd_proc        <- param_summary[3,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, alpha=alpha, beta=beta, sd_proc = 0)
  forecasts[,t] <- exp(x)
}
varI <- apply(forecasts,2,var)


##  Initial conditions and parameter uncertainty
x              <- sample(predictions[,nrow(bromus_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
alpha             <- params[sample_params,1]
beta              <- params[sample_params,2]
sd_proc       <- 0

forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, alpha=alpha, beta=beta, sd_proc = 0)
  forecasts[,t] <- exp(x)
}
varIP <- apply(forecasts,2,var)


##  Initial conditions, parameter, and process uncertainty
x              <- sample(predictions[,nrow(bromus_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
alpha              <- params[sample_params,1]
beta              <- params[sample_params,2]
sd_proc        <- params[sample_params,3]
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  x <- iterate_process(Nnow = x, alpha=alpha, beta=beta, sd_proc = sd_proc)
  forecasts[,t] <- exp(x)
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


