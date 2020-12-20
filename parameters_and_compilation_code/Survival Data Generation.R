library(coxed)
library(survival)
library(flexsurv)

### It is possible to specify your own hazard function! 

### You WILL have right-censored Data (meaning true control unit)

### Key assumption: Assume right-censored observations are "never treated" 

### So they can remain control units for the period of interest

## Here, durations = adoption time




hazardToPropScore <- function(hVec){
  
  ft <- rep(NA, length(hVec))

  ft[1] <- hVec[1]
  
  otherVec <- 1-hVec
  
  ft[2: length(ft)] <- hVec[2: length(ft)]*cumprod(otherVec)[1:(length(ft)-1)]
  
  return(ft)
  
}




my.hazard <- function(t){ 
  dnorm((log(t) - log(50))/log(10)) /
    (log(10)*t*(1 - pnorm((log(t) - log(50))/log(10))))
}


hazard <- function(t){
  
  
  .95*(t > 10)
  
}


# type="tvc" gives time-varying covariates


### Vary the sd argument to make the Xs more significant to the model
N <- 10000
Time <- 20

simdata <- sim.survdata(N=N, T=Time, num.data.frames=1, censor=0.8, xvars=2, sd=3,
                        censor.cond=TRUE, hazard.fun = hazard)

personPeriodFormatter <- function(simdata){

  unitIds <- seq(1, N, 1)
  
  actualData = simdata$data
  
  treatedUnit <- actualData$failed
  
  sampledTimes <- simdata$data$y
  
  sampledTimes[!treatedUnit] <- max(sampledTimes)+1
  
  
  Xs <- simdata$xdata
  
  unitTimeCov <- replicate(Time, Xs, simplify=FALSE)
  
  ### NOTE: Time is the third dimension here, so multiplication should be done with MARGIN=3
  covArray <- array(NA, dim=c(N, dim(Xs)[2], Time))
  
  for (i in 1:Time){ 
    
    covArray[,,i] <- as.matrix(unitTimeCov[[i]])
    
  }
  
  expandedUnitIds <- rep(unitIds, each=Time)
  
  expandedSampledTimes = rep(seq(1, Time, 1), length.out=dim(Xs)[1]*Time)
  
  expandedXs <- apply(Xs, MARGIN=2, FUN=rep, each=length(seq(1, Time, 1)))
  
  eventAtAndAfterIndicator <- W_maker(N=N, Time=Time, 
                                      ones_per_row = Time+1-sampledTimes)
  
  eventBeforeandAtIndicator <- t(sapply(pmin(Time, sampledTimes), FUN=function(x) 
    c(rep(1, x), rep(0, Time-x))))
  
  # * eventBeforeandAtIndicator
  
  eventAt <- eventAtAndAfterIndicator* eventBeforeandAtIndicator
  
  eventIndicator <- as.numeric(t(eventAt))
  
  
  survivalDataPersonPeriod <- data.frame(cbind(expandedUnitIds, expandedXs, expandedSampledTimes, eventIndicator))
  
  names(survivalDataPersonPeriod) <- c("id", paste( "X", seq(1, dim(Xs)[2], 1), sep=''), 'time', 'event')
  
  survivalDataPersonPeriod$time <- as.factor(survivalDataPersonPeriod$time)
  
  dataAndW <- list(survivalDataPersonPeriod, eventAtAndAfterIndicator, covArray)
  
  names(dataAndW) <- c('data', 'W', 'covArray')
  
  return(dataAndW)

}


#### Developing an appropriate model is next 

mod <- glm(event ~ .-1, data=survivalDataPersonPeriod[,-1])

plot(mod$fitted.values[1:20])

