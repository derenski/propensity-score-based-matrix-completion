

# For this analysis, time is in N (natural numbers)

failTimeSpikeGeo <- function(SpikeTime=10, p=.9, C=.2){
  
  IntegrationConst <- 1/(1-p)
  
  densFun <- function(t){p^(t-SpikeTime)*(t >= SpikeTime)/IntegrationConst}

  force(densFun)
  
  function(t){densFun(t)}
  
}

densityFormatter <- function(ft, C, Time){
  
  firstHalf <- ft(seq(1,Time,1))
  
  FirstHalfScaling <- C/sum(firstHalf)
  
  SecondHalfScaling <- (1-C)/(1-sum(firstHalf))
  
  forcedf <- force(ft)

  function(t){(FirstHalfScaling*(t <= Time)+ SecondHalfScaling*(t > Time))*ft(t)}

}




FtGenerator <- function(ft){
  
  cdf <- function(t){
    
    sapply(t, FUN=function(x) sum(ft(seq(1, x, 1))))
    
  }
  
  force(cdf)
  
  function(t){cdf(t)}
  
  
}


StGenerator <- function(Ft){
  
  StFun <- function(t){1-Ft(t)}
  
  force(StFun)
  
  function(t) StFun(t)
  
  
}


htGenerator <- function(ft){
  
  Ft <- FtGenerator(ft)
  
  St <- StGenerator(Ft)
  
  hazardRatio <- function(t){
    
    allVals <- ft(t)/St(t)
    
    return(allVals)
    
    }
  
  return(hazardRatio)
  
}

htTransformer <- function(ht, r){
  
  function(t){ht(t)*exp(r)/(1-(1-exp(r))*ht(t))}
  
  
}

betaGenerator <- function(K){
  
  unNormedVector <- runif(K, min=-1, max=1)
  
  return(array(unNormedVector/sum(unNormedVector), dim=c(K, 1)))
  
}


HtGenerator <- function(ht){
  
  Hfun <- function(t){
    
    sapply(t, FUN=function(x) sum(ht(seq(1, x, 1))))
    
  }
  
  force(Hfun)
  
  function(t){Hfun(t)}
  
}

baseHazardScaler <- function(r, Xs, ht, beta, ft, Time, C, St){
  
  relevantDensityVals <- ft(seq(1, Time-1, 1))
  
  relevantSurvivalVals <- St(seq(1, Time, 1))

  relevanthtVals <- ht(seq(1, Time-1, 1))
  
  productTerms <- t(apply(Xs, MARGIN=1, FUN=function(X) (1/(1+as.numeric(exp(array(X, dim=c(1, length(X))) 
                 %*% beta+r)-1)*relevanthtVals))))
  
  value <- mean(apply(productTerms, MARGIN=1, FUN=function(x) prod(x)))-((1-C)/relevantSurvivalVals[Time])
  
  return(value)
  
}

inverse = function (f, lower = 1, upper = 100) {
  
  func <- function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]$root
  
  actualEval <- function(ys){sapply(ys, FUN = func)}
  
  force(actualEval)
  
  function(ys){actualEval(ys)}  
  
}


FtForTimeGeneration <- function(X, beta, ht){
  
  FtAlt <- function(t) {
    
    if (t==1){
      
      return(0)
      
    }else{
    
    return(1-prod((1-ht(seq(1, t-1, 1)))/(1+(exp(as.numeric(X %*% beta))-1)*ht(seq(1, t-1, 1)))))
      
    }
      
    
  }
  
  force(FtAlt)
  
  function(ts) sapply(ts, FUN = FtAlt)
  
}

quantileFromCDF <- function(Ft){
  
  areWeDone <- FALSE
  
  finder <- function(p){
    
    areWeDone <- FALSE  
    
    i <- 1
    
    while (!areWeDone){
      
      foundIt <- p <= Ft(i)
      
      if ((!foundIt) &(i < 100)){
        
        i <- i+1
        
      }else{
        
        return(i)
        
      }
      
      
    }
    
  }
  
  return(finder)
  
}





## Next steps are to simulate data

coxSampler <- function(N, ft, C, Time, K){
  
  formatft <- densityFormatter(ft=ft, C=C, Time=Time)
  
  ourFt <- FtGenerator(formatft)
  
  ourSt <- StGenerator(ourFt)
  
  ourht <- htGenerator(formatft)
  
  H <- HtGenerator(ourht)
  
  Hinv <- inverse(H, lower=2, upper=500)
  
  ourBeta <- 3*betaGenerator(K)
  
  Xs <- abs(LaplacesDemon::rmvn(N, mu=rep(0, K)))
  
  rScaler <- uniroot(baseHazardScaler,Xs=Xs, ht=ourht, beta=ourBeta, 
                     ft=formatft, Time=Time, C=C, St=ourSt, interval=c(-100,100))[1]$root
  
  transformedht <- htTransformer(ht=ourht, r=rScaler)
  
  transformedHt <- HtGenerator(transformedht)
  
  Unifs <- runif(N)
  
  sampledTimes <- unlist(sapply(seq(1, N, 1), FUN=function(x)
    quantileFromCDF(FtForTimeGeneration(Xs[x, ], ourBeta, transformedht))(Unifs[x])))

  
  
  expandedSampledTimes = rep(seq(1, max(sampledTimes), 1), length.out=dim(Xs)[1]*max(sampledTimes))
  
  expandedXs <- apply(Xs, MARGIN=2, FUN=rep, each=length(seq(1, max(sampledTimes), 1)))
  
  eventAtAndAfterIndicator <- W_maker(N=N, Time=max(sampledTimes), 
                     ones_per_row = max(sampledTimes)+1-sampledTimes)
  
  eventBeforeandAtIndicator <- t(sapply(sampledTimes, FUN=function(x) 
    c(rep(1, x), rep(0, max(sampledTimes)-x))))
  
  eventAt <- eventAtAndAfterIndicator * eventBeforeandAtIndicator
  
  eventIndicator <- as.numeric(t(eventAt))
  
  
  survivalDataPersonPeriod <- data.frame(cbind(expandedXs, expandedSampledTimes, eventIndicator))
  
  names(survivalDataPersonPeriod) <- c(paste("X", seq(1, K, 1), sep=''), 'time', 'event')
  
  survivalDataPersonPeriod$time <- as.factor(survivalDataPersonPeriod$time)
  
  survivalDataCaseDuration <- data.frame(cbind(Xs, sampledTimes))
  
  
  names(survivalDataCaseDuration) <- c(paste("X", seq(1, K, 1), sep=''), 'adoption_time')
  
  
  finalOutput <- list(transformedht, survivalDataPersonPeriod, survivalDataCaseDuration, ourBeta)
  
  names(finalOutput) <- c('ht', 'survivalDataPersonPeriod', "survivalDataCaseDuration", 'beta')
  
  return(finalOutput)
  
}

firstft <- failTimeSpikeGeo(p=.7, C=.3)
modelParams <- coxSampler(N=500, ft=firstft, C=.2, Time=20, K=5)

(sum(modelParams$survivalDataCaseDuration$adoption_time <= 20)/
  length(modelParams$survivalDataCaseDuration$adoption_time))



survDat <- modelParams$survivalDataPersonPeriod


fitThing <- glm(event~ X1+time-1, data=survDat, family=binomial(link='logit'))



haz<-1/(1+exp(-coef(fitThing)))








