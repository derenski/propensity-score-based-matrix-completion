library(glmnet)
library(foreach)
library(doParallel)
library(dplyr)
library(reticulate)
library(jocre)

use_python("/usr/bin/python3", required = T)

server_name <- Sys.info()['nodename']

if (server_name == "Compy386"){

  source_python("/home/josh/GitHub/mnar_mc/mc_algorithms.py")

}else{

  source_python("/home/derenski/Github/mnar_mc/mc_algorithms.py")
  
}

treat.estimator <- function(Y, L.hat, W){ ## Efficient calculator of treatment estimator
  
  ## Y: Data matrix, L.hat: Estimated untreated matrix, W: Treatment indicator matrix
  
 upper.ones <- upper.tri(matrix(0, nrow=dim(Y)[2], ncol=dim(Y)[2]), diag = T)+0
 
 T.mat <- W %*% upper.ones
 
 treatment.estimator <- rep(NA, max(T.mat))
 
 for (i in 1:length(treatment.estimator)){
   
   treatment.estimator[i] <- mean(Y[T.mat==i]-L.hat[T.mat==i])
   
 } 
 
 return(treatment.estimator)
 
}


make_rho_mat <- function(rho,p){
  
  the_vec <- matrix(NA, nrow = p, ncol = p)
  
  for (i in 1:(p)){
    
    for (j in 1:p){
      
      the_vec[i,j] <- rho^abs(i-j)
      
    }
  }
  
  return(the_vec)
  
}

###### Helpers for parameter generation, and propensity score <=> weight conversion

propensity.score.to.weight <- function(prop.score.matrix, W){
  
  treated.score.weights <- 1/P_omega(prop.score.matrix, W=1-W)
  
  treated.score.weights[is.infinite(treated.score.weights)] <- 0
  
  untreated.score.weights <- 1/P_omega(1-prop.score.matrix, W=W)
  
  untreated.score.weights[is.infinite(untreated.score.weights)] <- 0
  
  final.weight.matrix <- treated.score.weights+untreated.score.weights
  
  return(final.weight.matrix)
  
}

propensityScoreToWeightConditional <- function(prop.score.matrix, W){
  
  treated.score.weights <- P_omega(W, W=1-W)
  
  untreated.score.weights <- P_omega(prop.score.matrix, W=W)/P_omega(1-prop.score.matrix, W=W)
  
  untreated.score.weights[is.infinite(untreated.score.weights) | is.nan(untreated.score.weights)] <- 0
  
  final.weight.matrix <- treated.score.weights+untreated.score.weights
  
  return(final.weight.matrix)
  
}

weight.to.propensity.score <- function(weight.matrix, W){
  
  prop.matrix <- matrix(NA, nrow=dim(weight.matrix)[1], ncol=dim(weight.matrix)[2])
  
  prop.matrix[W==1] <- 1/weight.matrix[W==1]
  
  prop.matrix[W==0] <- 1-(1/weight.matrix[W==0])
  
  return(prop.matrix)
  
}

###### Confidence Interval Methods


bs_percentile_method <- function(x, confidence_level){
  
    return(quantile(x, c((1-confidence_level)/2 , 
                         1-(1-confidence_level)/2), na.rm=T))
  
}


bs_z_method <- function(x, confidence_level){
  
  scaling_score <- qnorm(1-(1-confidence_level)/2)
  
  final_vec <- c(mean(x, na.rm=T)-scaling_score*sd(x, na.rm=T), 
                 mean(x, na.rm=T)+scaling_score*sd(x, na.rm=T))
  
  names(final_vec) <- c(paste(100*(1-confidence_level)/2, '%', sep=''), 
                        paste(100*(1-(1-confidence_level)/2), '%', sep=''))
  
  return(final_vec)
  
}

bs_cset_method <- function(x, confidence_level){
  
  thing <- cset(x, method="tost", alpha=1-confidence_level)
  
  conf_int <- t(thing$ci)
  
  rownames(conf_int) <- c(paste(100*(1-confidence_level)/2, '%', sep=''), 
                       paste(100*(1-(1-confidence_level)/2), '%', sep=''))
  
  return(conf_int)
  
}

bs_bc_method <- function(x, confidence_level){
  
  ## Produces 1-2alpha cis
  
  alpha <- (1-confidence_level)/2
  
  ecdf_func <- ecdf(x)
  
  param <- mean(x, na.rm=T)
  
  z0 <- qnorm(ecdf_func(param))
  
  z_low <- 2*z0+qnorm(alpha)
  
  z_high <- 2*z0+qnorm(1-alpha)
  
  conf_int <- c(quantile(x, pnorm(z_low)), quantile(x, pnorm(z_high)))
  
  names(conf_int) <- c(paste(100*(1-confidence_level)/2, '%', sep=''), 
                       paste(100*(1-(1-confidence_level)/2), '%', sep=''))
  
  return(conf_int)
  
}

bs_bca_method <- function(x, confidence_level){
  
  ## Produces 1-2alpha cis
  
  alpha <- (1-confidence_level)/2
  
  ecdf_func <- ecdf(x)
  
  param <- mean(x)
  
  z0 <- qnorm(ecdf_func(param))
  
  ## Acceleration constant just for expected value
  
  mu3 <- sum((x-param)^3)
  
  mu2 <- sum((x-param)^2)
  
  gamma <- mu3/(mu2^1.5)
  
  a <- (gamma/6)*sqrt(length(x))
  
  z_low <- z0 + (z0+qnorm(alpha))/(1-a*(z0+qnorm(alpha)))
  
  z_high <- z0 + (z0+qnorm(1-alpha))/(1-a*(z0+qnorm(1-alpha)))
  
  conf_int <- c(quantile(x, pnorm(z_low), na.rm=T), quantile(x, pnorm(z_high), na.rm=T))
  
  names(conf_int) <- c(paste(100*(1-confidence_level)/2, '%', sep=''), 
                        paste(100*(1-(1-confidence_level)/2), '%', sep=''))
  
  return(conf_int)
  
}




################## Bootstrap bias and calculation utilities

untreatedDataBootstrapSample <- function(Y, W, weight_matrix){

    controlUnits <- which(rowSums(W)==0)
    pretreatmentTimes <- which(colSums(W)==0)
    Y_untreated <- Y[controlUnits,]
    probsUntreated <- weight.to.propensity.score(weight_matrix, W=W)[controlUnits,]

    treatLengthEmpiricalDist <- table(rowSums(W))/dim(W)[1]
    
    WTraining <- array(0, dim=dim(Y_untreated))

    adoptionTimedist <- rev(treatLengthEmpiricalDist)
    
    
    completeAdoptionTimeDist <- rep(0, dim(Y_untreated)[2]+1)
    
    names(completeAdoptionTimeDist) <- 1:(dim(Y_untreated)[2]+1)

    names(adoptionTimedist) <- dim(Y_untreated)[2]-as.numeric(names(adoptionTimedist))+1

    completeAdoptionTimeDist[names(adoptionTimedist)] <- adoptionTimedist[names(adoptionTimedist)]

    numberToSamp <- ceiling(adoptionTimedist*dim(Y_untreated)[1])
    
    initialSamp <- rep(names(adoptionTimedist)[order(numberToSamp, decreasing = T)], 
                       times=numberToSamp[order(numberToSamp, decreasing = T)])
    
    lengthDifference <- length(initialSamp)-dim(WTraining)[1]
    
    if (lengthDifference != 0){
        
        if (lengthDifference > 0){
        
            initialSamp <- initialSamp[(lengthDifference+1):length(initialSamp)]
                
        }else{

            newSamp <- sample(initialSamp, size=lengthDifference, replace=TRUE)

            initialSamp <- c(initialSamp, newSamp)
        }
        
    }

    newTreatmentDesign <- as.numeric(sample(initialSamp, size=length(initialSamp), replace=FALSE))

    for (i in 1:dim(WTraining)[1]){

        WTraining[i, 1:dim(WTraining)[2] >= newTreatmentDesign[i]] <- 1

    }

    probAdoptTimeMatrix <- sapply(completeAdoptionTimeDist, FUN=rep, times=dim(WTraining)[1])[, 1:dim(WTraining)[2]]

    pijs <- t(apply(probAdoptTimeMatrix, MARGIN=1, FUN=cumsum))

    alteredProbAdopt <- (probsUntreated+
    pijs*(1-probsUntreated))

    alteredWeightsAdopt <- propensity.score.to.weight(alteredProbAdopt, W=WTraining)
    
    return(list(Y_untreated=Y_untreated, W_bootstrapped=WTraining,
               weights_bootstrapped=alteredWeightsAdopt))
    
    }

biasSampler <- function(Y, W, weight_matrix, numberOfBiasSamples=100){
    
    estimBiasValues <- array(NA, dim=c(numberOfBiasSamples, max(rowSums(W))))

    for (biasIteration in 1:numberOfBiasSamples){

        untreatedBootstrappedPieces = untreatedDataBootstrapSample(Y, W, weight_matrix)

        info <- completion_with_rank_estimation_validate_mu(Y=untreatedBootstrappedPieces$Y_untreated, 
                                                            W=untreatedBootstrappedPieces$W_bootstrapped,  
                                                            weight_matrix = untreatedBootstrappedPieces$weights_bootstrapped,
                                                                       initial_rank=40,
                                                                       tolerance=1e-04, 
                                                                       validation_max_iter=500,
                                                                       min_iter=10,
                                                                       max_iter=1000,
                                                                       mu_grid=0,
                                                                       K=5)

        estim <- treat.estimator(untreatedBootstrappedPieces$Y_untreated, info$L_hat, 
                                 untreatedBootstrappedPieces$W_bootstrapped)

        estimBiasValues[biasIteration, ] <- estim-0

        }
    
    return(estimBiasValues)
    
    }
bootstrapCI <- function(Y, W, weightMatrix, bootstrap_samps, method, ...){ ## Generate Boostrap estimates of treatment effect
  
  boot_ests <- array(NA, dim=c(bootstrap_samps, max(rowSums(W))))
  
  treatLengthDist <- rowSums(W)
  
  for (b in 1:bootstrap_samps){
    
    bSampIndicies <- c()
    
    for (treatLength in unique(treatLengthDist)){
      
      whatToSamp <- which(treatLengthDist==treatLength)
      
      if (length(whatToSamp)==1){
        
        bSampIndicies <- c(bSampIndicies, whatToSamp)
        
      }else{
        
        bSampIndicies <- c(bSampIndicies, sample(whatToSamp, size=length(whatToSamp), replace=T))
        
      }
      
    }
    
    Y_b_samp <- Y[bSampIndicies,]
    
    W_b_samp <- W[bSampIndicies,]
    
    weightMatrix_bsamp <- weightMatrix[bSampIndicies, ]
    
    info <- method(Y=Y_b_samp, W=W_b_samp, weightMatrix_bsamp, ...)
    
    estim <- treat.estimator(Y_b_samp, info$L_hat, W_b_samp)
    
    boot_ests[b, ] <- estim
    
  }
  
  return(boot_ests)
  
}


### Quick verification of bootstrap methods


we_testing <- FALSE

if (we_testing){
  
  num_boots_iters <- 500
  
  B <- 100
  
  good_coverage <- rep(NA, num_boots_iters)
  
  cis <- matrix(NA, nrow=2, ncol=num_boots_iters)
  
  for (n in 1:num_boots_iters){
    
    emp_data <- rep(NA, B)
    
    # the_data <- rnorm(1000, 3, 1)
    
    the_data <- rpois(1000, 3)
    
    for (b in 1:B){
      
      the.sample <- sample(the_data, length(the_data), replace=T)
      
      emp_data[b] <- mean(the.sample)
      
    }
    
    cis[,n] <- bc_method(emp_data, confidence_level=.95)
    
    good_coverage[n] <- (cis[1,n] <= 3 & 3 <= cis[2,n])
    
    print(n)
    
  }
  
  mean(good_coverage)
  
}

propensity_score_perturbation <- function(P, V){ ## Function for adding noise to propensity scores. 
  
  NMatrix <- ceiling(P*(1-P)/V)
  
  perturbedP <- rbinom(n=prod(dim(P)), size=as.numeric(NMatrix), prob=as.numeric(P))/NMatrix
  
  perturbedP <- matrix(perturbedP, nrow=dim(P)[1], ncol=dim(P)[2])
  
  perturbedP[V==0] <- P
  
  return(perturbedP)
  
}


#weight.to.propensity.score(final.weight.matrix, W=W)


prob_cloglog <- function(r, x){
  
  1-exp(-exp(r*x))
  
}

prob_logit <- function(r, x){
  
  1/(1+exp(x+r))
  
}
  
g_r <- function(func, r, x, C){ ## Controls proportion of TREATED units
  
  return(mean(  func(r,x)  -C)) ## Assumes positive covariates (which is very realistic in application)
  
}

parameterMakerBlock <- function(dist, N=100,Time=20, Time0=14, propTreat=.2, K=5,
                                  glmBeta=NULL, qFunc=qnorm, qFuncParms=list(0,1), 
                                  ...){ ## A helper function for generating parameters
  
  probFun <- function(r, xBeta){
    
    exp(xBeta+r)/(1+(exp(xBeta+r)))
    
    
  }
  
  rFun <- function(r, xBeta, C, lengthCutoff){
    
    vals <- probFun(r=r, xBeta=xBeta)
    
    return(mean(vals, na.rm=T)-C)
    
  }
  
  if (is.null(glmBeta)){
    
    glmBeta <- rep(1, K)
    
  }
  
  timeMarker <- seq(1, Time, 1)
  
  unitMarker <- seq(1, N, 1)
  
  unitTimeBase <- data.frame(cbind(rep(unitMarker, each=Time), rep(timeMarker, length.out=N*Time)),
                             stringsAsFactors = F)
  
  names(unitTimeBase) <- c('unit', 'time')
  
  Xs <- norta(n=N, corr_mat=diag(rep(1), K), qFunc=qFunc)
  
  XsforrEstimation <- norta(n=10000, corr_mat=diag(rep(1), K), qFunc=qFunc)
  
  r <- uniroot(rFun, xBeta= XsforrEstimation %*% glmBeta, C=propTreat,
               lower=-50, upper=50)$root
  
  fixedCovData <- data.frame(cbind(unitMarker, Xs),
                             stringsAsFactors = F)
  
  names(fixedCovData) <- c('unit', paste("X",seq(1, dim(Xs)[2], 1), sep=''))
  
  finalCovData <- merge(unitTimeBase, fixedCovData)
  
  unitTimeCov <- replicate(Time, Xs, simplify = F)
  
  ### NOTE: Time is the third dimension here, so multiplication should be done with MARGIN=3
  covArray <- array(NA, dim=c(N, dim(Xs)[2], Time))
  
  for (i in 1:Time){ 
    
    covArray[,,i] <- as.matrix(unitTimeCov[[i]])
    
  }
  
  propScore <- probFun(r=r, xBeta=Xs %*% glmBeta)
  
  doWeTreat <- rbern(N, p=propScore)
  
  propScoreMatrix <- t(sapply(propScore, FUN=function(x) rep(x, Time)))
  
  propScoreMatrix[,1:(Time0)] <- 0
  
  startTimes <- c(Time+1, Time0+1)[doWeTreat+1]
  
  W <- W_maker(N, Time, ones_per_row=Time+1-startTimes)
  
  trueWeightMatrix <- propensity.score.to.weight(propScoreMatrix, W)
  
  allUnitsAndTimes <- expand.grid(timeMarker, unitMarker)
  
  responseUnmatrixed <- apply(allUnitsAndTimes, MARGIN=1, FUN=function(x) W[x[2], x[1]])
  
  finalCovData$treated <- responseUnmatrixed
  
  dataForModeling <- data.frame(unitMarker, Xs, doWeTreat)
  
  names(dataForModeling)[c(1, dim(dataForModeling)[2])] <- c('unit', 'treated')
  
  mod <- glm(treated ~ ., data=dataForModeling[,-1], family=binomial)
  
  propScoreEst <- t(sapply(mod$fitted.values, FUN=rep, each=Time))
  
  propScoreEst[,1:(Time0)] <- 0
    
  weightMatrixEst <- propensity.score.to.weight(propScoreEst, W)
                              
  weightMatrixConditionalEst <- propensityScoreToWeightConditional(propScoreMatrix, W)
  
  finalStuff <- list(trueWeightMatrix, propScoreEst, weightMatrixEst, weightMatrixConditionalEst, W, covArray)
  
  names(finalStuff) <- c("trueWeightMatrix", "propScoreEst", "weightMatrixEst", "weightMatrixConditionalEst", "W", "covArray")
  
  return(finalStuff)
  
}


### TO still control probability of treatment, make the average of all these numbers=propTreat!
# adoption_distribution_matrix <- t(sapply(rbeta(1000, shape1=1, shape2=3), FUN= function(x) 
#  adoption_distribution_geo_decay(seq(1, Time, 1), decay=x, T0=10, Time=20)))

#average_adoption_curve <- apply(adoption_distribution_matrix, MARGIN=2, FUN=mean)
                              

parameterMakerStaggered <- function(N, Time=20, propTreat=.2, generatingModel='Bin', K, 
                                      lengthCutoff=10, glmBeta=NULL,
                                      qFunc=qnorm, qFuncParms = list(1,1), ...){ 
    
    
    
  
  if (is.null(glmBeta)){
    
    glmBeta <- rep(1, K)
    
  }
  
  ### generatingModel can be "Bin" or "Poi"
  fullyParametricProbModel <- function(N=N, Time=Time, generatingModel='Poi', propTreat=propTreat, K, theBeta,
                                       lengthCutoff=10,
                                       qFunc=qnorm){
    
    probFun <- function(r, xBeta){
      
      exp(xBeta+r)/(1+(exp(xBeta+r)))
      
      
    }
      
    rFunBin <- function(r, xBeta, C, lengthCutoff){
      
      vals <- (1-probFun(r=r, xBeta=xBeta))^lengthCutoff
      
      return(mean(vals, na.rm=T)-(1-C))
      
    }
    
    rFunPois <- function(r, xBeta, C, lengthCutoff){
      
      vals <- exp(-1*exp(xBeta+r))
      
      return(mean(vals, na.rm=T)-(1-C))
      
    }
    
    all_function_names <- as.character(ls.str())
    
    rFunNames <- all_function_names[str_detect(all_function_names, pattern="rFun")]
    
    chosenrFun <- get(rFunNames[str_detect(rFunNames,  pattern=generatingModel)])
    
    ### Truncated Poisson or binomial work here
    
    Xs <- norta(n=N, corr_mat=diag(rep(1), K), qFunc=qFunc)
    
    unitTimeCov <- replicate(Time, Xs, simplify = F)
    
    ### NOTE: Time is the third dimension here, so multiplication should be done with MARGIN=3
    covArray <- array(NA, dim=c(N, dim(Xs)[2], Time))
    
    for (i in 1:Time){ 
      
      covArray[,,i] <- as.matrix(unitTimeCov[[i]])
      
    }
    
    XsforEstimation <- norta(n=10000, corr_mat=diag(rep(1), K), qFunc=qFunc)
    
    r <- uniroot(chosenrFun, xBeta= XsforEstimation %*% theBeta, C=propTreat,
                 lengthCutoff=lengthCutoff,
                 lower=-500, upper=500)$root

    
    if (generatingModel=='Bin'){
    
    actualProbTreat <- probFun(r=r, xBeta=Xs %*% theBeta)
    
    treatmentLengths <- rbinom(n=N, size=lengthCutoff, prob=actualProbTreat)
    
    probAdopt <- t(sapply(actualProbTreat, FUN=function(x) rev(dbinom(seq(0, Time, 1), 
                                                                      size=lengthCutoff, prob=x))))
    
    }else{
      
      lambdaParams <- exp(Xs %*% theBeta+r)
      
      treatmentLengths <- pmin(rpois(n=N, lambda=lambdaParams), lengthCutoff)
      
      probAdopt <- t(sapply(lambdaParams, FUN=function(x) rev(dpois(seq(0, Time, 1), lambda=x))))
      
    }

    adoptionTime <- Time-treatmentLengths+1
    
    eventAtAndAfterIndicator <- W_maker(N=N, Time=Time, 
                                        ones_per_row = Time+1-adoptionTime)
    
    #eventBeforeandAtIndicator <- t(sapply(pmin(Time, sampledTimes), FUN=function(x) 
    #  c(rep(1, x), rep(0, Time-x))))
    
    successes <- treatmentLengths
    failures <- lengthCutoff-treatmentLengths
    
    didTheyAdopt <-treatmentLengths != 0
    
    caseDurationData <- data.frame(cbind(Xs, adoptionTime, didTheyAdopt, successes, failures))
    
    names(caseDurationData) <- c(paste('X', seq(1, K, 1), sep=''), 'y', 'failed', 'successes', 'failures')
    
    truePropScoreMat <- t(apply(probAdopt, MARGIN=1, FUN=cumsum))[,1:Time]
    
   # truePropScoreMat <- probAdopt[, 1:Time]
    
  #  truePropScoreMat <- sapply(1:N, FUN=function(x) c(rep(0, Time-lengthCutoff),
   #                                                rep(dbinom(treatmentLengths[x], size=lengthCutoff, prob=actualProbTreat[x]), lengthCutoff)))
    
    #truePropScoreMat <- matrix(truePropScoreMat, nrow=N, ncol=Time, byrow=T)
    
    finalOutput <- list(truePropScoreMat, caseDurationData, eventAtAndAfterIndicator, covArray)
    
    names(finalOutput) <- c("truePropScoreMat", "caseDurationData", "W", "covArray")
    
    return(finalOutput)
    
  }
  
  simulatedDataPack <- fullyParametricProbModel(N=N, Time=Time, propTreat=propTreat, K=K, theBeta=glmBeta, 
                                            generatingModel=generatingModel, lengthCutoff=lengthCutoff, 
                                            qFunc=qFunc)
  
  simulatedData <-  simulatedDataPack$caseDurationData
  
  truePropScores <-  simulatedDataPack$truePropScoreMat
  
  covArray <- simulatedDataPack$covArray
  
  if (generatingModel=='Bin'){
    
    NeededColumnPersonPeriodModeling <- !(names(simulatedData) %in% c('y', 'failed'))

    mod <- glm(cbind(successes, failures) ~ ., data=simulatedData[, NeededColumnPersonPeriodModeling], family=binomial)

    probAdoptEst <- t(sapply(mod$fitted.values, FUN=function(x) rev(dbinom(seq(0, Time, 1), 
                                                                         size=lengthCutoff, prob=x))))[,1:Time]
  
  }else{
    
    NeededColumnPersonPeriodModeling <- !(names(simulatedData) %in% c('y', 'failures', 'failed'))
    
    mod <- glm(successes ~ ., data=simulatedData[, NeededColumnPersonPeriodModeling], family=poisson)
    
    probAdoptEst <- t(sapply(mod$fitted.values, FUN=function(x) rev(dpois(seq(0, Time, 1), lambda=x))))[,1:Time]
    
  }
  
  propScoreEst <- t(apply(probAdoptEst, MARGIN=1, FUN=cumsum))
  
  W <- simulatedDataPack$W
  
  trueWeightMatrix <- propensity.score.to.weight(truePropScores, W)
  
  perturbedFinalWeightMatrix <- propensity.score.to.weight(propScoreEst, W)
                             
  weightMatrixConditionalEst <- propensityScoreToWeightConditional(truePropScores, W)
  
  
  finalStuff <- list(W, propScoreEst, trueWeightMatrix, perturbedFinalWeightMatrix, weightMatrixConditionalEst, covArray)
  
  names(finalStuff) <- c("W", "propScoreEst", "trueWeightMatrix", "weightMatrixEst", "weightMatrixConditionalEst", "covArray")

  return(finalStuff)
  
}

                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
parameterMakerStaggeredProvidedCovariates <- function(Time , modelCovariates, propTreat=.2, generatingModel='Bin',
                                      lengthCutoff=10, glmBeta=NULL, ...){   
    
  N <- dim(modelCovariates)[1] 
  
  K <- dim(modelCovariates)[2]  
  
  if (is.null(glmBeta)){
    
    glmBeta <- rep(1, K)
    
  }
  
  ### generatingModel can be "Bin" or "Poi"
  fullyParametricProbModel <- function(N=N, Time=Time, generatingModel='Poi', propTreat=propTreat, K=k, theBeta,
                                       lengthCutoff=10, modelCovariates){
    
    probFun <- function(r, xBeta){
      
      exp(xBeta+r)/(1+(exp(xBeta+r)))
      
      
    }
      
    rFunBin <- function(r, xBeta, C, lengthCutoff){
      
      vals <- (1-probFun(r=r, xBeta=xBeta))^lengthCutoff
      
      return(mean(vals, na.rm=T)-(1-C))
      
    }
    
    rFunPois <- function(r, xBeta, C, lengthCutoff){
      
      vals <- exp(-1*exp(xBeta+r))
      
      return(mean(vals, na.rm=T)-(1-C))
      
    }
    
    all_function_names <- as.character(ls.str())
    
    rFunNames <- all_function_names[str_detect(all_function_names, pattern="rFun")]
    
    chosenrFun <- get(rFunNames[str_detect(rFunNames,  pattern=generatingModel)])
    
    ### Truncated Poisson or binomial work here
    
    Xs <- modelCovariates
    
    unitTimeCov <- replicate(Time, Xs, simplify = F)
    
    ### NOTE: Time is the third dimension here, so multiplication should be done with MARGIN=3
    covArray <- array(NA, dim=c(N, dim(Xs)[2], Time))
    
    for (i in 1:Time){ 
      
      covArray[,,i] <- as.matrix(unitTimeCov[[i]])
      
    }
    
    XsforEstimation <- modelCovariates
    
    r <- uniroot(chosenrFun, xBeta= XsforEstimation %*% theBeta, C=propTreat,
                 lengthCutoff=lengthCutoff,
                 lower=-500, upper=500)$root

    
    if (generatingModel=='Bin'){
    
    actualProbTreat <- probFun(r=r, xBeta=Xs %*% theBeta)
    
    treatmentLengths <- rbinom(n=N, size=lengthCutoff, prob=actualProbTreat)
    
    probAdopt <- t(sapply(actualProbTreat, FUN=function(x) rev(dbinom(seq(0, Time, 1), 
                                                                      size=lengthCutoff, prob=x))))
    
    }else{
      
      lambdaParams <- exp(Xs %*% theBeta+r)
      
      treatmentLengths <- pmin(rpois(n=N, lambda=lambdaParams), lengthCutoff)
      
      probAdopt <- t(sapply(lambdaParams, FUN=function(x) rev(dpois(seq(0, Time, 1), lambda=x))))
      
    }

    adoptionTime <- Time-treatmentLengths+1
    
    eventAtAndAfterIndicator <- W_maker(N=N, Time=Time, 
                                        ones_per_row = Time+1-adoptionTime)
    
    #eventBeforeandAtIndicator <- t(sapply(pmin(Time, sampledTimes), FUN=function(x) 
    #  c(rep(1, x), rep(0, Time-x))))
    
    successes <- treatmentLengths
    failures <- lengthCutoff-treatmentLengths
    
    didTheyAdopt <-treatmentLengths != 0
    
    caseDurationData <- data.frame(cbind(Xs, adoptionTime, didTheyAdopt, successes, failures))
    
    names(caseDurationData) <- c(paste('X', seq(1, K, 1), sep=''), 'y', 'failed', 'successes', 'failures')
    
    truePropScoreMat <- t(apply(probAdopt, MARGIN=1, FUN=cumsum))[,1:Time]
    
   # truePropScoreMat <- probAdopt[, 1:Time]
    
  #  truePropScoreMat <- sapply(1:N, FUN=function(x) c(rep(0, Time-lengthCutoff),
   #                                                rep(dbinom(treatmentLengths[x], size=lengthCutoff, prob=actualProbTreat[x]), lengthCutoff)))
    
    #truePropScoreMat <- matrix(truePropScoreMat, nrow=N, ncol=Time, byrow=T)
    
    finalOutput <- list(truePropScoreMat, caseDurationData, eventAtAndAfterIndicator, covArray)
    
    names(finalOutput) <- c("truePropScoreMat", "caseDurationData", "W", "covArray")
    
    return(finalOutput)
    
  }
  
  simulatedDataPack <- fullyParametricProbModel(N=N, Time=Time, propTreat=propTreat, K=K, theBeta=glmBeta, 
                                            generatingModel=generatingModel, lengthCutoff=lengthCutoff, 
                                            modelCovariates=modelCovariates)
  
  simulatedData <-  simulatedDataPack$caseDurationData
  
  truePropScores <-  simulatedDataPack$truePropScoreMat
  
  covArray <- modelCovariates
  
  if (generatingModel=='Bin'){
    
    NeededColumnPersonPeriodModeling <- !(names(simulatedData) %in% c('y', 'failed'))

    mod <- glm(cbind(successes, failures) ~ ., data=simulatedData[, NeededColumnPersonPeriodModeling], family=binomial)

    probAdoptEst <- t(sapply(mod$fitted.values, FUN=function(x) rev(dbinom(seq(0, Time, 1), 
                                                                         size=lengthCutoff, prob=x))))[,1:Time]
  
  }else{
    
    NeededColumnPersonPeriodModeling <- !(names(simulatedData) %in% c('y', 'failures', 'failed'))
    
    mod <- glm(successes ~ ., data=simulatedData[, NeededColumnPersonPeriodModeling], family=poisson)
    
    probAdoptEst <- t(sapply(mod$fitted.values, FUN=function(x) rev(dpois(seq(0, Time, 1), lambda=x))))[,1:Time]
    
  }
  
  propScoreEst <- t(apply(probAdoptEst, MARGIN=1, FUN=cumsum))
  
  W <- simulatedDataPack$W
  
  trueWeightMatrix <- propensity.score.to.weight(truePropScores, W)
  
  perturbedFinalWeightMatrix <- propensity.score.to.weight(propScoreEst, W)
                             
  weightMatrixConditionalEst <- propensityScoreToWeightConditional(truePropScores, W)
  
  
  finalStuff <- list(W, propScoreEst, trueWeightMatrix, perturbedFinalWeightMatrix, weightMatrixConditionalEst, covArray)
  
  names(finalStuff) <- c("W", "propScoreEst", "trueWeightMatrix", "weightMatrixEst", "weightMatrixConditionalEst", "covArray")

  return(finalStuff)
  
}
                       
                             
                             
                             
                             
                             
                             
parameterMakerBlockProvidedCovariates <- function(Time , modelCovariates, propTreat=.2, Time0=90, glmBeta=NULL, ...){   
    
  N <- dim(modelCovariates)[1] 
  
  K <- dim(modelCovariates)[2]  
  
  if (is.null(glmBeta)){
    
    glmBeta <- rep(1, K)
    
  }
  
      probFun <- function(r, xBeta){
    
    exp(xBeta+r)/(1+(exp(xBeta+r)))
    
    
  }
  
  rFun <- function(r, xBeta, C, lengthCutoff){
    
    vals <- probFun(r=r, xBeta=xBeta)
    
    return(mean(vals, na.rm=T)-C)
    
  }
  
  timeMarker <- seq(1, Time, 1)
  
  unitMarker <- seq(1, N, 1)
  
  unitTimeBase <- data.frame(cbind(rep(unitMarker, each=Time), rep(timeMarker, length.out=N*Time)),
                             stringsAsFactors = F)
  
  names(unitTimeBase) <- c('unit', 'time')
  
  Xs <- modelCovariates
  
  r <- uniroot(rFun, xBeta= Xs %*% glmBeta, C=propTreat,
               lower=-500, upper=500)$root
  
  fixedCovData <- data.frame(cbind(unitMarker, Xs),
                             stringsAsFactors = F)
  
  names(fixedCovData) <- c('unit', paste("X",seq(1, dim(Xs)[2], 1), sep=''))
  
  finalCovData <- merge(unitTimeBase, fixedCovData)
  
  unitTimeCov <- replicate(Time, Xs, simplify = F)
  
  ### NOTE: Time is the third dimension here, so multiplication should be done with MARGIN=3
  covArray <- array(NA, dim=c(N, dim(Xs)[2], Time))
  
  for (i in 1:Time){ 
    
    covArray[,,i] <- as.matrix(unitTimeCov[[i]])
    
  }
  
  propScore <- probFun(r=r, xBeta=Xs %*% glmBeta)
  
  doWeTreat <- rbern(N, p=propScore)
  
  propScoreMatrix <- t(sapply(propScore, FUN=function(x) rep(x, Time)))
  
  propScoreMatrix[,1:(Time0)] <- 0
  
  startTimes <- c(Time+1, Time0+1)[doWeTreat+1]
  
  W <- W_maker(N, Time, ones_per_row=Time+1-startTimes)
  
  allUnitsAndTimes <- expand.grid(timeMarker, unitMarker)
  
  responseUnmatrixed <- apply(allUnitsAndTimes, MARGIN=1, FUN=function(x) W[x[2], x[1]])
  
  finalCovData$treated <- responseUnmatrixed
  
  dataForModeling <- data.frame(unitMarker, Xs, doWeTreat)
  
  names(dataForModeling)[c(1, dim(dataForModeling)[2])] <- c('unit', 'treated')
  
  mod <- glm(treated ~ ., data=dataForModeling[,-1], family=binomial)
  
  propScoreEst <- t(sapply(mod$fitted.values, FUN=rep, each=Time))
  
  propScoreEst[,1:(Time0)] <- 0
                              
  trueWeightMatrix <- propensity.score.to.weight(propScoreMatrix, W)
    
  weightMatrixEst <- propensity.score.to.weight(propScoreEst, W)
                              
  weightMatrixConditionalEst <- propensityScoreToWeightConditional(propScoreMatrix, W)
  
  finalStuff <- list(W,  propScoreEst, trueWeightMatrix, weightMatrixEst, weightMatrixConditionalEst, covArray)
                     
  names(finalStuff) <- c("W", "propScoreEst", "trueWeightMatrix", "weightMatrixEst", "weightMatrixConditionalEst", "covArray")

  return(finalStuff)
    
  }
                            
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
                             
### Calculates Xbeta for the outcome model
covariateEffectComponent <- function(covArray, betaEffect){
  
  if (dim(covArray)[2] != length(betaEffect)){
    
    stop('Incorrect dimension for betaEffect')
    
  }
  
  value <- apply(covArray, MARGIN=3, FUN=function(x) x %*% betaEffect)
  
  return(value)
  
}



#### The Effect Size Functions

delta_t_updown <- function(t, arg_max=3, y_max=7, cutoff=3, ...){ ## Max is at exp(mu-1)
  
  implied_mu <- log(arg_max)+1
  
  scaling_factor <- y_max*exp(.5+log(arg_max))
  
  f_t <- scaling_factor*(1/t)*exp(-.5*(log(t)-implied_mu)^2)
  
  f_t[cutoff < cutoff & t > arg_max] <- cutoff
  
  return(f_t)
  
}


delta_t_decay <- function(t, y_max=6, halfway_time=4, cutoff=3, ...){ 
  
  implied_lambda <- -log(2)/(halfway_time-1)
  
  f_t <- y_max*exp(implied_lambda*(t-1))
  
  f_t[f_t < cutoff] <- cutoff
  
  return(f_t)
  
}


delta_t_plateau <- function(t, y_max=5, halfway_time=3, ...){ 
  
  implied_lambda <- -log(2)/halfway_time
  
  f_t <-  y_max*(1-exp(implied_lambda*(t)))
  
  return(f_t)
  
}

delta_t_constant <- function(t, value=10, ...){ 
  
  return(rep(value, length(t)))
  
}

treated_matrix_creator <- function(x, f_of_t=delta_t_constant,...){
  
  treatment_row <- rep(0, length(x))
  
  treatment_times <- which(x==1)
  
  if (length(treatment_times)!=0){
    
    treatment_row[treatment_times] <- f_of_t(treatment_times-(min(treatment_times)-1),...)
    
  }
  
  return(treatment_row)
  
}

W_maker <- function(N, Time, ones_per_row){
  
  W <- t(sapply(ones_per_row, FUN=function(x) c(rep(0, Time-x),
                                                rep(1,  x)) ))
  
  return(W)
  
}

### A Copula Method for Generating 
### Multidimensional covariates

norta <- function(n, corr_mat, qFunc=qnorm){
  
  Zs <- rmvn(n=n, mu=rep(0, dim(corr_mat)[1]), Sigma=corr_mat)
  
  phi_Zs <- pnorm(Zs)
  
  desired_variables <- matrix(NA, nrow=n, ncol=dim(corr_mat)[1])
    
  for (row in 1:dim(desired_variables)[1]){
      
    desired_variables[row, ] <- qFunc(phi_Zs[row, ])
    
    }
      
    return(desired_variables)
    
}



### A Copula Method for Generating 
### Serial correlation in non-gaussian settings

nortaNoise <- function(number, corr_mat, desired_mean_matrix, 
                  distribution='gaussian', scalar_sigma=1,
                  df=3){
  
  Zs <- rmvn(n=number, mu=rep(0, dim(corr_mat)[1]), Sigma=corr_mat)
  
  phi_Zs <- pnorm(Zs)
  
  desired_variables <- matrix(NA, nrow=number, ncol=dim(corr_mat)[1])
  
  if (distribution == 'poisson'){
    
    for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- qpois(phi_Zs[row, ], 
                                        lambda=desired_mean_matrix[row, ])
      
    }
    
  }else if(distribution == 'scaled_gamma'){
    
    for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- qgamma(phi_Zs[row, ], 
                                         shape=desired_mean_matrix[row, ], rate=1)
      
    }
    
  }else if (distribution == 'exponential'){
    
    for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- qexp(phi_Zs[row, ], 
                                       rate=1/desired_mean_matrix[row, ])
      
    }
    
  } else if (distribution=='t'){
    
    if (df < 3){
      
      print("Warning: May Not Have Desired Covariance Structure")
      
    }
    
    for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- (qt(phi_Zs[row, ], df=df)
                                   +desired_mean_matrix[row, ])
      
    }
    
  }else if (distribution == 'gaussian'){
    
    for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- qnorm(phi_Zs[row, ], 
                                        mean=desired_mean_matrix[row, ], sd=scalar_sigma)
      
    }
    
  } 
  
  
  return(desired_variables)
  
}



S_tau <- function(x, tau=0){
  
  return(sign(x)*pmax(abs(x)-tau, 0))
  
}

P_omega <- function(A, W){ ## Helper function for projection
  ## Projects onto {(i,j) : W_{(i,j)} == 0}
  
  A[W!=0] <- 0
  
  return(A)
  
}


D_tau <- function(x, tau=0){
  
  svd_x <- svd(x)
  
  return(svd_x$u %*% S_tau(diag(svd_x$d), tau=tau) %*% t(svd_x$v))
  
}































##############################################################################################################
##############################################################################################################


##### L completion with automatic rank estimation (unknown rank)

#### NOTE: This function is never used in this code, as it is a special case of the weighted algorithm

completion_with_rank_estimation <- function(Y, W, r_init=40, 
                                            tolerance=1e-04, 
                                            min_iter=10,
                                            max_iter=1000,
                                            mu=1e-10){
  
  shrink_operator <- function(x, mu){ ## A helper function, the shrinkage operator
    
    if (x > 1*mu){
      
      return(x-mu)
      
    }else if (abs(x)< mu){
      
      return(0)
      
    }else{
      
      return(x+mu)
      
    }
    
  }
  
  N <- dim(Y)[1]
  
  Time <- dim(Y)[2]
  
  Us <- rmvn(n=N, mu=rep(0, r_init))
  
  Us <- apply(Us, MARGIN=2, FUN = function(x) x/norm(x, '2'))
  
  W_vec <- as.numeric(rmvn(n=1, mu=rep(0,r_init)))
  
  Vs <- rmvn(n=Time, mu=rep(0,r_init))
  
  Vs<- apply(Vs, MARGIN=2, FUN = function(x) x/norm(x, '2'))
  
  r <- r_init
  
  L_k <- P_omega(A=Y, W=W)
  
  iterating <- TRUE
  
  iter_number <- 1
  
  while (iterating){
    
    L_k_r <- L_k
    
    if(iter_number==1){
      
      parts_to_update <- 1:length(W_vec)
      
    }else{
      
      parts_to_update <- which(W_vec > 0)
      
    }
    
    for (r_number in parts_to_update){ ## Note: weights don't influence update of u,v (because we normalize)
      
      Us[, r_number] <- (L_k_r %*% Vs[, r_number])
      
      Us[, r_number] <- Us[, r_number]/norm(Us[, r_number], '2')
      
      Vs[, r_number] <- (t(L_k_r) %*% Us[, r_number])
      
      Vs[, r_number] <- Vs[, r_number]/norm(Vs[, r_number], '2')
      
      W_vec[r_number] <- max(0, shrink_operator(x=sum((as.matrix(Us[, r_number]) 
                                                       %*%  t(Vs[, r_number])) * L_k_r), mu=mu))
      
      L_k_r <- L_k_r - W_vec[r_number] * (Us[, r_number] %*% t(Vs[, r_number]))
      
    }
    
    L_k_plus_1 <- L_k
    
    Z <- L_k-L_k_r
    
    L_k_plus_1[W!=0] <- Z[W!=0]
    
    threshold_number <- (10^-3)*(sum(W==0)/prod(dim(W)))*sum(W_vec)
    
    W_vec[W_vec < threshold_number] <- 0
    
    first_num <- (norm(P_omega(L_k_plus_1
                               -Z, W=W), 'F')
                  /norm(P_omega(L_k_plus_1, W=W)))
    ## Is rank r approx close to original matrix?
    
    
    second_num <- (norm(L_k_plus_1-L_k, 'F')
                   /norm(L_k_plus_1))
    
    ## Is updating changing much? 
    
    if (((first_num < tolerance) | (second_num  < tolerance) |
         (iter_number >= max_iter)) & !(iter_number < min_iter)){
      
      L_k_final <- Z
      
      iterating <- FALSE
      
      break
      
    }else{
      
      L_k <- L_k_plus_1
      
      iter_number <- iter_number+1
      
      next
      
    }
    
    
  }
  
  threshold_number <- (10^-3)*(sum(W==0)/prod(dim(W)))*sum(W_vec)
  
  R_star <- sum(W_vec > threshold_number)
  
  final_output <- list(Z, R_star)
  
  names(final_output) <- c("L_hat", 'rank_estimate')
  
  return(final_output)
  
}


#################

## A weighted version of the above algorithm

completion_with_rank_estimation_weighted <- function(Y, W, weight_matrix=NULL, r_init=40, 
                                                     tolerance=1e-04, 
                                                     min_iter=10,
                                                     max_iter=1000,
                                                     mu=1e-10){
  ### Do NOT root the weight matrix!! 
  
  if(is.null(weight_matrix)){
    
    weight_matrix <- matrix(1, nrow=dim(Y)[1], ncol=dim(Y)[2])
    
  }
    
  shrink_operator <- function(x, mu){ ## A helper function, the shrinkage operator
    
    if (x > 1*mu){
      
      return(x-mu)
      
    }else if (abs(x)< mu){
      
      return(0)
      
    }else{
      
      return(x+mu)
      
    }
    
  }
  
  N <- dim(Y)[1]
  
  Time <- dim(Y)[2]
  
  Us <- rmvn(n=N, mu=rep(0, r_init))
  
  Us <- apply(Us, MARGIN=2, FUN = function(x) x/norm(x, '2'))
  
  W_vec <- as.numeric(rmvn(n=1, mu=rep(0,r_init)))
  
  Vs <- rmvn(n=Time, mu=rep(0,r_init))
  
  Vs<- apply(Vs, MARGIN=2, FUN = function(x) x/norm(x, '2'))
  
  r <- r_init
  
  L_k <- P_omega(A=Y, W=W)
  
  iterating <- TRUE
  
  iter_number <- 1
  
  while (iterating){
    
    L_k_r <- L_k
    
    if(iter_number==1){
      
      parts_to_update <- 1:length(W_vec)
      
    }else{
      
      parts_to_update <- which(W_vec > 0)
      
    }
    
    for (r_number in parts_to_update){ ## Note: weights don't influence update of u,v (because we normalize)
      
      ## Two algorithms should be identical when matrix is just 1's (unless initialized differently)
      
      Us[, r_number] <- (1/rowSums(weight_matrix))*((weight_matrix*L_k_r) %*% Vs[, r_number])
      
      Us[, r_number] <- Us[, r_number]/norm(Us[, r_number], '2')
      
      Vs[, r_number] <-(1/colSums(weight_matrix))*(t(weight_matrix*L_k_r) %*% Us[, r_number])
      
      Vs[, r_number] <- Vs[, r_number]/norm(Vs[, r_number], '2')
      
      mult_thing <- (as.matrix(Us[, r_number]) 
                     %*%  t(Vs[, r_number]))
      
      special_number <- sum(weight_matrix*(mult_thing^2))
      
      ### Weights are getting too big
      
      # mu/special_number
      
      W_vec[r_number] <- max(0, shrink_operator(x=sum(mult_thing * L_k_r*weight_matrix)/special_number, 
                         mu=mu))
      
    #  print(mean(W_vec))
      
      L_k_r <- L_k_r - W_vec[r_number] * (Us[, r_number] %*% t(Vs[, r_number]))
      
    }
    
    L_k_plus_1 <- L_k
    
    Z <- L_k-L_k_r
    
    L_k_plus_1[W!=0] <- Z[W!=0]
    
  #  threshold_number <- (10^-3)*(sum(W==0)/prod(dim(W)))*sum(W_vec)
    
  #  W_vec[W_vec < threshold_number] <- 0
    
    
    ### Weighting 
    
   # first_num <- (norm(P_omega(weight_matrix*(L_k_plus_1
    #                           -Z), W=W), 'F')
  #                /norm(P_omega(L_k_plus_1, W=W)))
    ## Is rank r approx close to original matrix?
    
    
   # second_num <- (norm(weight_matrix*(L_k_plus_1-L_k), 'F')
   #                /norm(L_k_plus_1))

      first_num <- (norm(P_omega(weight_matrix*(L_k_plus_1-Z), W=W), 'F')
                /norm(P_omega(weight_matrix, W=W), 'F'))
    
    
     second_num <- (norm(weight_matrix*(L_k_plus_1-L_k), 'F')/norm(weight_matrix, 'F'))
                    
    
    ## Is updating changing much? 
    
    if (((first_num < tolerance) | (second_num  < tolerance) |
         (iter_number >= max_iter)) & !(iter_number < min_iter)){
      
      L_k_final <- Z
      
      iterating <- FALSE
      
      break
      
    }else{
      
      L_k <- L_k_plus_1
      
      iter_number <- iter_number+1
      
      next
      
    }
    
    
  }
  
  # (10^-3)
  
  threshold_number <- (10^-3)*(sum(W==0)/prod(dim(W)))*sum(W_vec)
  
  R_star <- sum(W_vec > threshold_number)
  
  final_output <- list(Z, R_star)
  
  names(final_output) <- c("L_hat", 'rank_estimate')
  
  return(final_output)
  
}


####################################################





#### Completion with rank estimation, validating mu

completion_with_rank_estimation_validate_mu <- function(Y, W, weight_matrix=NULL,
                                                        initial_rank=20,
                                                        tolerance=1e-03, 
                                                        validation_max_iter=300,
                                                        min_iter=10,
                                                        max_iter=1000,
                                                        mu_grid=10^seq(-2, 4,1),
                                                        K=4){
  if(is.null(weight_matrix)){ ## If weight matrix not specified, defaults to unweighted
    
    weight_matrix <- matrix(1, nrow=dim(Y)[1], ncol=dim(Y)[2])
    
  }
  
  if (sum(is.na(Y)) > 0){
    
    Y[is.na(Y)] <- 0
    
  }
  
  doing_validation <- length(mu_grid) > 1
  
  if (doing_validation){

    model_errors <- foreach(potential_mu=mu_grid,
                                       .packages=c('LaplacesDemon'),
                                       .combine = 'c',
                  .export = c('completion_with_rank_estimation_weighted',
                              'P_omega', 'propensity.score.to.weight', 
                              'weight.to.propensity.score')) %dopar% {
    
      mean_squared_errors_this_k <- rep(NA, K)
    
      for (fold in 1:K){ ## Are they doing cv correctly?
      
        card_O <- sum(W==0)
        
        prob_chosen <- 1-(card_O/prod(dim(Y)))
        
       # prob_chosen <- weight.to.propensity.score(weight_matrix, W=W)
        
        choose_samp <- matrix(rbern(prod(dim(W)), prob_chosen), nrow=dim(W)[1], ncol=dim(W)[2])
        
        W_train <- W+choose_samp*(1-W)
        
        cv_prop_matrix <- weight.to.propensity.score(weight_matrix, W=W)+prob_chosen*(1-weight.to.propensity.score(weight_matrix, W=W))
        
        cv_weight_matrix <- propensity.score.to.weight(cv_prop_matrix, 
                                                       W=W_train)*(1-W_train)
        
        completion_this_mu <- completion_with_rank_estimation_weighted(Y, W_train, 
                                                    weight_matrix=cv_weight_matrix,
                                                    r_init=initial_rank, 
                                                    tolerance=tolerance, 
                                                    min_iter=min_iter,
                                                    max_iter=validation_max_iter,
                                                    mu=potential_mu)
      
        L_k <- completion_this_mu$L_hat
        
        ## *cv_weight_matrix

        mean_squared_errors_this_k[fold] <- (norm((1-W)*sqrt(cv_weight_matrix)*(W_train)*(Y-L_k), 'F')
                                             /sqrt(sum((1-W)*(W_train))))
      
      }
    
      mean(mean_squared_errors_this_k)
      
    
    }
    
    chosen_mu <- mu_grid[which.min(model_errors)]
    
  }else{
    
    chosen_mu <- mu_grid
    
  }
  
  the_final_weight_matrix <- weight_matrix*(1-W)
  
  ### YOU CHANGE THE WEIGHTS HERE 
  #the_final_weight_matrix <- propensity.score.to.weight(the_final_weight_matrix, 
  #                           W=matrix(0, nrow=dim(W)[1], ncol=dim(W)[2]))*(1-W)
  
  #the_final_weight_matrix <- weight_matrix
  
  completion_chosen_mu <- completion_with_rank_estimation_weighted(Y, W, 
                                                 weight_matrix=the_final_weight_matrix,
                                                                   r_init=initial_rank, 
                                                                   tolerance=tolerance, 
                                                                   min_iter=min_iter,
                                                                   max_iter=max_iter,
                                                                   mu=chosen_mu)
  
  final_output <- list(completion_chosen_mu$L_hat, 
                       completion_chosen_mu$rank_estimate,
                       chosen_mu)
  
  names(final_output) <- c("L_hat", 'rank_estimate', 'chosen_mu')
  
  return(final_output)
  
}


##### Athey-based method
### Implementation of Athey et al's matrix recovery method. 

matrix_completion_causal <- function(Y, W, num_iter=1000, K=5, 
                            lambda_grid=c(0, 10^seq(-20, 0, 1)),
                            tolerance=1e-03){
  
  ### Cross Validation Section 
  
  card_O <- sum(W==0)
  
  fold_size <- floor((card_O^2)/prod(dim(Y)))
  
  mean_squared_errors_lambdas <- foreach(potential_lambda=lambda_grid,
          .packages=c('LaplacesDemon'),
          .combine = 'c',
          .export = c('P_omega', 'D_tau', 'S_tau')) %dopar%{
    
    mean_squared_errors_this_k <- rep(NA, K)
    
    for (fold in 1:K){ ## Are they doing cv correctly?
      
      untreated_cells <- which(W==0)
      
      train_fold <- sample(untreated_cells, size=fold_size, replace=F)
      
      card_O_k <- length(train_fold)
      
      val_fold <- untreated_cells[!(untreated_cells %in% train_fold)]
      
      P_O_Y <- P_omega(Y, W=W)
      
    #  P_O_Y[(W==1)] <- 0
      
      P_O_Y[(untreated_cells %in% val_fold)] <- 0
      
      L_k <- P_O_Y
      
      first_L_k <- L_k
      
      for (iter_num in 1:num_iter){
        
        P_O_perp_L_k <- P_omega(L_k, W=1-W)
        
     #   P_O_perp_L_k[W==0] <- 0
        
        P_O_perp_L_k[(untreated_cells %in% train_fold)] <- 0
        
        L_k_plus_1 <- D_tau(x=P_O_Y+P_O_perp_L_k, tau=.5*potential_lambda*fold_size)
        
        L_difference <- mean((L_k_plus_1[W==1]-L_k[W==1])^2)
        
        if (L_difference > tolerance){
          
          L_k <- L_k_plus_1
          
          next
          
        }else{
          
          L_k <- L_k_plus_1
          
          break
          
        }
        
        
      }
      
      mean_squared_errors_this_k[fold] <- mean((L_k[val_fold]-Y[val_fold])^2, 
                                               na.rm=T)
      
    }
    
    mean(mean_squared_errors_this_k, na.rm=T)
    
  }
  
  cv_lambda <- lambda_grid[which.min(mean_squared_errors_lambdas)]
  
  card_O <- sum(W==0)
  
  P_O_Y <- P_omega(Y, W=W)
  
 # P_O_Y[W==1] <- 0
  
  L_k <- P_O_Y
  
  for (iter_num in 1:num_iter){
    
    P_O_perp_L_k <- P_omega(L_k, W=1-W)
    
    # P_O_perp_L_k[W==0] <- 0
    
    L_k_plus_1 <- D_tau(x=P_O_Y+P_O_perp_L_k, tau=.5*cv_lambda *card_O)
    
    L_difference <- mean((L_k_plus_1[W==1]-L_k[W==1])^2)
    
    if (L_difference > tolerance){
      
      L_k <- L_k_plus_1
      
      next
      
    }else{
      
      L_k <- L_k_plus_1
      
      break
      
    }
    
  }
  
  output <- list(L_k, cv_lambda)
  
  names(output) <- c("L_hat", "cv_lambda")
  
  return(output)
  
}






weightedSoftImpute_validate_lambda <- function(Y, W, weight_matrix, num_iter=1000, K=5, 
                                               lambda_grid=c(0, 100),
                                               tolerance=1e-03){
  
  ### Cross Validation Section 
  
  card_O <- sum(W==0)
  
  fold_size <- floor((card_O^2)/prod(dim(Y)))
  
  mean_squared_errors_lambdas <- c()
  
  untreated_units <- which(apply(W, MARGIN=1, FUN=function(x) all(x==0)))
  
  for (potential_lambda in lambda_grid){
                                           
    mean_squared_errors_this_k <- rep(NA, K)
                                           
    for (fold in 1:K){ ## Are they doing cv correctly?
      
      prob_chosen <- 1-(card_O/prod(dim(Y)))
      
      #prob_chosen <- weight.to.propensity.score(weightMat, W=W)
      
      choose_samp <- matrix(rbern(prod(dim(W)), prob_chosen), nrow=dim(W)[1], ncol=dim(W)[2])
      
      W_train <- W+choose_samp*(1-W)
      
      cv_prop_matrix <- weight.to.propensity.score(weight_matrix, W=W)+prob_chosen*(1-weight.to.propensity.score(weight_matrix, W=W))
      
      cv_weight_matrix <- propensity.score.to.weight(cv_prop_matrix, 
                                                     W=W_train)
      
      L_k <- weighted_softimpute(X=np_array(Y), M=np_array(1-W_train), W=np_array(cv_weight_matrix), lmbda=potential_lambda,
                                 apg_max_iter=as.integer(num_iter), apg_eps=tolerance,
                                 apg_use_restart=TRUE)

      mean_squared_errors_this_k[fold] <- (norm((1-W)*cv_weight_matrix*(W_train)*(Y-L_k), 'F')
                                           /sqrt(sum((1-W)*(W_train))))
      
    } 
    
    mean_squared_errors_lambdas <- c(mean_squared_errors_lambdas,
                                     mean(mean_squared_errors_this_k, na.rm=T))
    
    }
  
  cv_lambda <- lambda_grid[which.min(mean_squared_errors_lambdas)]
  
  L_k <- weighted_softimpute(X=np_array(Y), M=np_array(1-W), W=np_array(weight_matrix), lmbda=cv_lambda,
                             apg_max_iter=as.integer(num_iter), apg_eps=tolerance,
                             apg_use_restart=TRUE)
  
  dimnames(L_k) <- dimnames(Y)                               
                                 
  output <- list(L_k, cv_lambda)
  
  names(output) <- c("L_hat", "cv_lambda")
  
  return(output)
  
}

















### Competitors for our method: 





### Completion with factor model

completion_factor_model <- function(Y,W, propScoreMat, numFactors=3){
  
  Y_proj <- P_omega(Y, W)
  
  theMeans <- matrix(rowSums(Y_proj)/(dim(W)[2]-rowSums(W)), nrow=dim(Y_proj)[1], ncol=1)
  
  unmeand_cov <- (Y_proj %*% t(Y_proj))-(theMeans%*%t(theMeans))
  
  times_where_both_observed <- (1-W) %*% t(1-W) 
  
  altered_cov_est <- unmeand_cov/times_where_both_observed
  
  cov_est_for_PCA <- (1/dim(W)[1])*altered_cov_est
  
  eigen_decomp <- prcomp(cov_est_for_PCA, center=F, tol=1e-20)
  
  lambda_i <- (unname(eigen_decomp$rotation[])*sqrt(dim(Y)[1]))[, 1:numFactors]

  firstPart <- t((1-W)/(1-propScoreMat)*Y)

  firstPart[is.na(firstPart)] <- 0
  
  F_t <- firstPart %*% lambda_i/dim(Y)[1]
  
  L_hat <- lambda_i %*% t(F_t)
  
  return(L_hat)

  }














