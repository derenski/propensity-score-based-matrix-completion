library(ggplot2)
library(LaplacesDemon)
library(glmnet)
library(foreach)
library(doParallel)
library(quadprog)
library(openxlsx)
library(stringr)
library(dplyr)
library(reticulate)
library(jocre)
library(knitr)
library(reshape)

use_python("/usr/bin/python3", required = T)

### For setting working directory
##setwd("../../../../Causal Inference with Propensity Scores/simulations_and_reports/parameters_and_compilation_code")


server_name <- Sys.info()['nodename']

params <- read.xlsx("parameters_and_descriptions.xlsx",
                            sheet="parameter_data", rowNames = T)

source('causal_inference_methods_code_corrupted.R')

mse_and_se_of_mse <- function(error_mat){
  
  squared_errors <- error_mat
  
  the_mse <- mean(squared_errors, na.rm=T)
  
  se_mse <- sqrt((sd(squared_errors, na.rm=T)^2)/prod(dim(squared_errors)))
  
  final_stuff <- c(the_mse, se_mse)
  
  names(final_stuff) <- c("mse", "se_mse")
  
  return(final_stuff)
  
}

### For setting working directory
##setwd("../../../../Causal Inference with Propensity Scores/simulations_and_reports/parameters_and_compilation_code")


twoEntriesToValParenth <- function(twoVals){ ### For formatting distance tables to LATEX
  
  paste(twoVals[1]," (",twoVals[2],")", sep="")
  
}

number_of_Es <- as.numeric(params["number_of_Es", ])

bootstrap_samps <- as.numeric(params["bootstrap_samps", ])

N <- as.numeric(params["N", ])

Time <- as.numeric(params["Time", ])

K <- as.numeric(params["K", ])

generatingModel <- params['generatingModel',]

lengthCutoff <- as.numeric(params['lengthCutoff',])

matrix_type <- c("tall", 'wide')[(Time > N)+1]

Time0 <- as.numeric(params["Time0", ])

rho_parameter <- as.numeric(params["rho_parameter", ])

betaNorm <- as.numeric(params["beta_norm", ])

tau <- as.numeric(params["tau", ])

propTreat <- as.numeric(params["prop_treat", ])

sigma_squared <- as.numeric(params["sigma_squared", ])

min_iter <- as.numeric(params["min_iter", ])

max_iter <- as.numeric(params["max_iter", ])

tolerance <- as.numeric(params["tolerance", ])

L_scaling <- as.numeric(params["L_scaling", ])

arg_max <- as.numeric(params["arg_max", ])

y_max <- as.numeric(params["y_max", ])

halfway_time <-as.numeric( params["halfway_time", ])

cutoff <- as.numeric(params["cutoff", ])

design <- params["design", ]

all_function_names <- as.character(ls.str())

treatment_function_names <- all_function_names[str_detect(all_function_names, pattern="delta_t_")]

bs_method_names <- all_function_names[str_detect(all_function_names, pattern="bs_")]

prob_func_names <- all_function_names[str_detect(all_function_names, pattern="prob_(?!func)")]

max_lag <- as.numeric(params["max_lag", ])

group_gap <- as.numeric(params["group_gap", ])

desired_coverage <- as.numeric(params["desired_coverage", ])

treatment_function <- get(treatment_function_names[str_detect(treatment_function_names, 
                                                              pattern=params["treatment_effect_function", ])])

bootstrap_method <- get(bs_method_names[str_detect(bs_method_names, 
                                                            pattern=params["bootstrap_method", ])])

correct_bias <- as.logical(params["correct_bias", ])

link_func <- params["link_func", ]

prob_func <- get(prob_func_names[str_detect(prob_func_names, 
                                 pattern=link_func)])
##################### Fixed Parameter Simulations
##################################################
###################################################
#################################################

mus.chosen <- rep(NA, number_of_Es)


#betaEffect <- rep(10, K)
#betaForGLM <- rep(3, K)

#covDist <- abs(rnorm(100000))^1.5

### Okay scenarios: 

#### betaEffect <- rep(10, K), betaForGLM <- rep(1, K)
## covDist <- abs(rnorm(100000))^1.5

# betaEffect <- rep(5, K)
#betaForGLM <- rep(2, K)
# covDist <- abs(rnorm(100000))^.5


#betaEffect <- rep(5, K) 
#betaForGLM <- rep(1, K)
#covDist <- abs(rnorm(100000))

#### When covariates have outsized effect on resulting time series, IPS, FACTOR are best
### In this case, treat/untreated covariate distributions overlap, 
#betaEffect <- rep(50, K) 
#betaForGLM <- rep(1, K)
#covDist <- abs(rnorm(100000))

#betaEffect <- rep(2, K) 
#betaForGLM <- rep(.05, K)
#covDist <- abs(rnorm(100000))

#betaEffect <- rep(5, K) 
#betaForGLM <- rep(50, K)
#covDist <- abs(rnorm(100000))
#### When treated/untreated groups produce very different time series, weightedSoftImpute and FACTOR are best


betaEffect <- (betaNorm/K)*rep(1, K) 
betaForGLM <- rep(.5, K)
covDist <- rexp(10000, rate=1)

#rnorm(10000, mean=c(-2,2), sd=.75)

### Key assumption: no hard separation between treated, untreated groups

qarb <- function(emp_dist){
  
  quantile_func <- function(p){
  
    initial_dist <- quantile(emp_dist, p, na.rm=TRUE)
    
    initial_dist[is.infinite(initial_dist)] <- 0
    
    return(initial_dist)
  
  }
  
  return(quantile_func)
  
}

#covDist <- rlaplacep(100000, mu=c(-1, 1), tau=1)

### Log standard-normal is a good distribution for covariates


theQuantileFunction <- qarb(covDist)
if (design=="staggered_adoption"){ 
  
  observed.params <- parameterMakerStaggered(N=N, Time=Time, propTreat=propTreat, 
                                               generatingModel = "Bin",
                                               K=K, lengthCutoff=lengthCutoff,
                                               glmBeta=betaForGLM,
                                               qFunc=theQuantileFunction)
  
  # norm(observed.params$trueWeightMatrix-observed.params$weightMatrixEst, 'F')/sqrt(prod(dim(observed.params$weightMatrixEst)))

  xBeta <- covariateEffectComponent(covArray=observed.params$covArray, betaEffect)
  
  W <- observed.params$W
  
}else if (design=="block_treatment"){
  
  # mean=c(-.5*group_gap, .5*group_gap)
  
  observed.params <- parameterMakerBlock(N=N, Time=Time, Time0=Time0, propTreat=propTreat, 
                                           k=K,
                                           link_func=link_func,
                                           glmBeta=betaForGLM,
                                           qFunc=theQuantileFunction, 
                                           prob_func = prob_func)
  
  xBeta <- covariateEffectComponent(covArray=observed.params$covArray, betaEffect)
  
  W <- observed.params$W
  
}

final.weight.matrix <- observed.params$weightMatrixEst

tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, 
                      f_of_t=treatment_function, arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, value=tau))

treatment_times <- sort(unique(dim(W)[2]-rowSums(W)+1))

treatment_times <- treatment_times[treatment_times <= dim(W)[2]]

delta_t <- treatment_function(1:(Time-(min(treatment_times)-1)),
                              arg_max=arg_max, 
                              y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, value=tau)

autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])


norm((final.weight.matrix-W)*(1-W), 'F')

sig_to_noise_ratios <- c()

desired_clusters <- 8

max_available_clusters <- detectCores()

cl <- makeCluster(min(c(max_available_clusters, desired_clusters))-3)

registerDoParallel(cl)

  # *rnorm(N*Time)
  
L <- matrix(rnorm(N*Time, mean=5, sd=3), nrow=N, ncol=Time)
  
# L <- matrix(rexp(N*Time, .1), nrow=N, ncol=Time)

#L <- matrix(5, nrow=N, ncol=Time)  

### Single run comparison

  # set.seed(3729)
  
  errors_this_L_mc_nnm <- rep(NA, number_of_Es)
  
  errors_this_L_weightedSoftImpute <- rep(NA, number_of_Es)
  
  errors_this_L_completionFactorModel <- rep(NA, number_of_Es)

  errors_this_L_r1Comp <- rep(NA, number_of_Es)
  
  errors_this_L_weightedR1Comp <- rep(NA, number_of_Es)
  
  errors_this_L_oracle <- rep(NA, number_of_Es)
  
  ### To get factor estimator to perform badly, increase proportion of treated units 

  for (j in 1:number_of_Es){
    
    print(j)
    
    Y <- nortaNoise(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*observed.params$W+xBeta, 
                 distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))

    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=4, 
                                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                                            tolerance=1e-04)
    
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    final.weight.matrix.normalized <- final.weight.matrix/sum(final.weight.matrix)
    
    
    ### Optional 1bitMC approach for propensity score estimation: 
    
    ## one_bit_MC_fully_observed(M=W, link=std_logistic_function, 
    ## link_gradient=grad_std_logistic_function, tau=1, gamma=1)
    
    #### max_rank=None, min_value=None, max_value=None
    ### Need to validate lambda, use your validation method to do so. 
    
    
  #  allWSIErrors <- c()
    
 #   for (lamb in seq(0, 2000, 100)){
    
 #   weightedSoftImputeInfo <- weighted_softimpute(X=np_array(Y), M=np_array(W), W=np_array(final.weight.matrix), lmbda=lamb,
  #                                              apg_max_iter=as.integer(1000), apg_eps=1e-6,
  #                                                apg_use_restart=TRUE)
    
  #  L_weightedSoftImpute <- weightedSoftImputeInfo
    
 #   tau_estimate_weightedSoftImpute <- treat.estimator(Y=Y, L.hat=L_weightedSoftImpute, W=W)
    
  #  allWSIErrors <- c(allWSIErrors, mean(abs(tau_estimate_weightedSoftImpute-delta_t)^2))
    
   # print(lamb)
    
 #   }
    
    weightedSoftImputeInfo <- weightedSoftImpute_validate_lambda(Y, W, weight_matrix=final.weight.matrix, num_iter=1000, K=5, 
                               lambda_grid=seq(0, 2000, 100), tolerance=1e-03)
    
    L_weightedSoftImpute <- weightedSoftImputeInfo$L_hat
    
    L_completionFactorModel <- completion_factor_model(Y=Y, W=W, propScoreMat = observed.params$propScoreEst,
                                                       numFactors=rankMatrix(mc_nnm_info$L_hat)[1])
    
    ### Keep in mind that the weight matrix you feed in is not the exact matrix used for weighting (See paper)
      
      
      
      
    r1Comp_info <- completion_with_rank_estimation_validate_mu(Y=Y, W=W,  
                                                               weight_matrix = 
                                                               array(1, dim=dim(Y)),
                                                               initial_rank=40,
                                                               tolerance=1e-04, 
                                                               validation_max_iter=500,
                                                               min_iter=10,
                                                               max_iter=1000,
                                                               mu_grid=c(0, 10^seq(-3,3,1)),
                                                               K=6)  
      
    L_r1Comp <- r1Comp_info$L_hat 
      
      
      
    weightedR1Comp_info <- completion_with_rank_estimation_validate_mu(Y=Y, W=W,  
                                                               weight_matrix = final.weight.matrix,
                                                               initial_rank=40,
                                                               tolerance=1e-04, 
                                                               validation_max_iter=500,
                                                               min_iter=10,
                                                               max_iter=1000,
                                                               mu_grid=c(0, 10^seq(-3,3,1)),
                                                               K=6)
    # c(0, 10^seq(-4, 3, 1))
    
    
    mus.chosen[j] <- weightedR1Comp_info$chosen_mu
    
    L_weightedR1Comp <- weightedR1Comp_info$L_hat 
    
    ## Oracle in the sense that we know the untreated counterfactual
    
    ## What are we estimating now?? It is no longer L, because that is just a baseline matrix
    
    ## We are actually performing completion on Y (to get untreated counterfactuals)
    
    ## Why does using low-rank approximations work?
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
    
    tau_estimate_weightedSoftImpute <- treat.estimator(Y=Y, L.hat=L_weightedSoftImpute, W=W)
    
    tau_estimate_completionFactorModel <- treat.estimator(Y=Y, L_completionFactorModel, W=W)
      
      
    tau_estimate_r1Comp <- treat.estimator(Y=Y, L.hat=L_r1Comp, W=W)
    
    tau_estimate_weightedR1Comp <- treat.estimator(Y=Y, L.hat=L_weightedR1Comp, W=W)
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L+xBeta, W=W)
    
    error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    #error_tau_weightedSoftImpute <- min(allWSIErrors)
    
    error_tau_weightedSoftImpute <- mean(abs(tau_estimate_weightedSoftImpute-delta_t)^2)
    
    error_tau_completionFactorModel <- mean(abs(tau_estimate_completionFactorModel-delta_t)^2)
      
    error_tau_r1Comp <- mean(abs(tau_estimate_r1Comp-delta_t)^2)
    
    error_tau_weightedR1Comp <- mean(abs(tau_estimate_weightedR1Comp-delta_t)^2)
    
    error_tau_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)
    
    errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
    
   # noiseChooser <- which.min(allWSIErrors)+floor(rnorm(1, 1))
    
  #  if (noiseChooser > length(allWSIErrors)){
      
  #    noiseChooser <- length(allWSIErrors)
      
  #  }else if (noiseChooser < 0){
      
  #  noiseChooser <- 1
      
  #  }
    
  # errors_this_L_weightedSoftImpute[j] <- allWSIErrors[noiseChooser]
    
    errors_this_L_weightedSoftImpute[j] <- error_tau_weightedSoftImpute
    
    errors_this_L_completionFactorModel[j] <- error_tau_completionFactorModel
      
    errors_this_L_r1Comp[j] <- error_tau_r1Comp

    errors_this_L_weightedR1Comp[j] <- error_tau_weightedR1Comp
    
    errors_this_L_oracle[j] <- error_tau_oracle
    
  }
  
  


error_data <- data.frame(rbind(round(mse_and_se_of_mse(errors_this_L_mc_nnm), 3), 
                               round(mse_and_se_of_mse(errors_this_L_weightedSoftImpute), 3),
                               round(mse_and_se_of_mse(errors_this_L_completionFactorModel), 3),
                               round(mse_and_se_of_mse(errors_this_L_r1Comp), 3),
                               round(mse_and_se_of_mse(errors_this_L_weightedR1Comp), 3),
      round(mse_and_se_of_mse(errors_this_L_oracle), 3)), stringsAsFactors = F)

error_data <- cbind(c("MC-NNM", "Weighted softImpute", "FACTOR", 'R1Comp',
                      "Weighted R1Comp", "Oracle"),error_data )

names(error_data) <- c("Method", "MSE", "SE")

effect_plot <- (ggplot(NULL, aes(x=1:length(delta_t), y=delta_t)) 
                + geom_line() + theme_bw() + xlab("Time") + ylab(expression(delta[t]))
                +ggtitle(""))


errorDisplayName <- twoEntriesToValParenth(names(error_data)[2:3])

error_data[, errorDisplayName] <- apply(error_data[,2:dim(error_data)[2]], MARGIN=1, 
                                                    FUN=function(x) twoEntriesToValParenth(x))

errorForDisplay <- error_data[, c(1,4)]



texErrorTable <- kable(errorForDisplay, "latex")





#### Leave none arguments alone unless you intend to specify them


############################################
############################################
###########################################
##### Bootstrapping

### Uses same parameters as in point estimate section                                        
                                        
                                        
boot_ests <- array(NA, dim=c(bootstrap_samps, max(rowSums(W)), number_of_Es))

delta_t_estimates <- array(NA, dim=c(number_of_Es, length(delta_t)))
  
for (e in 1:number_of_Es){
    
  Y <- nortaNoise(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*observed.params$W+xBeta, 
                 distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
    
 # boot_ests[, , e] <- bootstrappedEstimateGenerator(Y=Y, W=W, weightMatrix=final.weight.matrix, 
  #                                                  B1=B1, B2=B2)
    
    
    #boot_ests[, , e] <- bootstrapCI(Y=Y, W=W, weightMatrix=final.weight.matrix, 
    #            bootstrap_samps=bootstrap_samps, method=completion_with_rank_estimation_validate_mu,
    #           initial_rank=40, tolerance=1e-04, validation_max_iter=500,
    #            min_iter=10, max_iter=1000, mu_grid=0, K=5)
    
    
  boot_ests[, , e] <- bootstrapCI(Y=Y, W=W, weightMatrix=final.weight.matrix, 
                bootstrap_samps=bootstrap_samps, method=completion_with_rank_estimation_validate_mu,
               initial_rank=40, tolerance=1e-04, validation_max_iter=500,
                min_iter=10, max_iter=1000, mu_grid=0, K=5)
    
  if (correct_bias){
    
      theBias <- biasSampler(Y=Y, W=W, weight_matrix=final.weight.matrix, 
                             numberOfBiasSamples=500)


      boot_ests[, , e] <- boot_ests[, , e]-theBias[sample(1:dim(theBias)[1],
                                                          size=dim(boot_ests)[1], 
                                                          replace=T),]
      
      }
    
  fullEstimateinfo <- completion_with_rank_estimation_validate_mu(Y=Y, W=W,  
                                                               weight_matrix = final.weight.matrix,
                                                               initial_rank=40,
                                                               tolerance=1e-04, 
                                                               validation_max_iter=500,
                                                               min_iter=10,
                                                               max_iter=1000,
                                                               mu_grid=0,
                                                               K=5)
    
    
    
    fullEstimator <- treat.estimator(Y, fullEstimateinfo$L_hat, W)
    
    delta_t_estimates[e,] <-   fullEstimator
    
    print(e)
    
    }
#### Scenario setup


# bootstrap_method
#### BOOTSTRAP DIAGNOSTICS  
                    
### FROM BOOTSTRAP TO FULL DATA              
options(stringsAsFactors=F)                    
                    
                    
bootSimNumber <- 1


deltaTConsidered <- delta_t_estimates[bootSimNumber,]
                    
bootstrapMethodConsidered <-  boot_ests[,,bootSimNumber]                    
                    
iterationBootEstimates <- data.frame(melt(bootstrapMethodConsidered))   
                    
names(iterationBootEstimates) <- c("Iteration", "Time", "Value")                     
         
                    
                    
iterationBootEstimates <- iterationBootEstimates %>% mutate(Iteration = as.numeric(
    str_extract(as.character(Iteration), "[0-9]+$")),
    Time = as.numeric(str_extract(as.character(Time), "[0-9]+$")) )                   
                    
                    
delta_t_estimate_data <- cbind.data.frame(1:length(deltaTConsidered), deltaTConsidered)                    
                    
names(delta_t_estimate_data) <- c("Time", "Value")

ggplot(iterationBootEstimates, aes(x=Time, y=Value, color=factor(Iteration))) + geom_line()+
    geom_line(data=delta_t_estimate_data, aes(x=Time, y=Value, color='black'), lwd=2, color='black')+ theme_bw(base_size=20) + theme(legend.position = "none")                     
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
    
                    
                    
### FROM FULL DATA TO TRUTH                    
                    
#iterationBootAverages <- data.frame(melt(delta_t_estimates[,,methodConsidered])) 
iterationBootAverages <- data.frame(melt(delta_t_estimates[, ]))                     
names(iterationBootAverages) <- c("Iteration", "Time", "Value")   
                    
iterationBootAverages <- iterationBootAverages %>% mutate(Iteration = as.numeric(
    str_extract(as.character(Iteration), "[0-9]+$")),
    Time = as.numeric(str_extract(as.character(Time), "[0-9]+$")) )                    

delta_t_data <- cbind.data.frame(1:length(delta_t), delta_t)
                    
names(delta_t_data) <- c("Time", "Value")

ggplot(iterationBootAverages, aes(x=Time, y=Value, color=factor(Iteration))) + geom_line()+
    geom_line(data=delta_t_data, aes(x=Time, y=Value, color='black'), lwd=2, color='black')+ theme_bw(base_size=20) + theme(legend.position = "none") 

                
### FROM BOOTSTRAP TO TRUTH                    

ggplot(iterationBootEstimates, aes(x=Time, y=Value, color=factor(Iteration))) + geom_line()+
    geom_line(data=delta_t_data, aes(x=Time, y=Value, color='black'), lwd=2, color='black')+ theme_bw(base_size=20) + theme(legend.position = "none")                          
                    
                    
extra_mses <- iterationBootAverages %>% group_by(Iteration) %>% summarize(mean((Value-delta_t)^2))            
                    
                    
                    
                    
                    
                    
### z_method works well consistently, percentile, bc, and bca method works really well for large sample size

#estimatedCoverages <- array(NA, dim=c(length(dimnames(boot_ests)[[4]]), length(delta_t)),
#                           dimnames=list(dimnames(boot_ests)[[4]], 1:length(delta_t)))
               
#for (methodConsidered in dimnames(estimatedCoverages)[[1]]){                    
                    
    if (! str_detect(params["bootstrap_method", ], 'cset')){

     # ci_array <- apply(boot_ests[ , , , methodConsidered], MARGIN=c(2, 3), FUN=bootstrap_method ,desired_coverage)
        
     ci_array <- apply(boot_ests[, , ], MARGIN=c(2, 3), FUN=bootstrap_method ,
                       1-(1-desired_coverage)/length(delta_t))

    }else{

     # ci_array <- plyr::aaply(boot_ests[ , , , methodConsidered],  3,bs_cset_method ,desired_coverage, .drop=FALSE)
      ci_array <- plyr::aaply(boot_ests[, , ],  3,bs_cset_method ,
                              1-(1-desired_coverage)/length(delta_t), .drop=FALSE)
      ci_array <- aperm(ci_array, perm=c(2,3,1))

    }

    betweenTwoNumbers <- function(number, numberSandwich){

        return(numberSandwich[1] <= number & number <= numberSandwich[2])

    }

    we_do_good <- matrix(NA, nrow=dim(ci_array)[3], ncol=dim(ci_array)[2])

    for (iterationIndex in 1:dim(ci_array)[3]){

        whereWeHit <- c()

        for (colIndex in 1:dim(ci_array)[2]){

            theBounds <- ci_array[, colIndex, iterationIndex]

            whereWeHit <- as.vector(c(whereWeHit, betweenTwoNumbers(delta_t[colIndex], theBounds)))

        }

        we_do_good[iterationIndex,] <- whereWeHit

    }                    

    estimated_coverage <- colMeans(we_do_good)
    
 #   estimatedCoverages[methodConsidered, ] <- estimated_coverage
    
#    }

confidenceDataFrame <- data.frame(t(estimated_coverage))

rownames(confidenceDataFrame)[1] <- 'Estimated Coverage'

names(confidenceDataFrame) <- paste("Time", 1:dim(we_do_good)[2], sep=' ')            

coverageTable <- kable(confidenceDataFrame, "latex")

## Coverage should be about 1-alpha-O(1/sqrt(N))

## How many untreated units are there? How many treated units are there? 

## You must consider the bootstrap bias from each sample individually. 

timesGoodEachIteration <- rowSums(we_do_good)                    
 
chosenExample <- which(timesGoodEachIteration == max(timesGoodEachIteration))[1]                    
                    
ciExample <- cbind.data.frame(1:length(delta_t), t(ci_array [,,chosenExample]))
names(ciExample) <- c("Time", "Lower", "Upper")
                    
ciExamplePlot <- ggplot(ciExample, aes(x=Time, y = Upper)) +
  geom_line(aes(y = Lower), linetype='dashed') + 
  geom_line(aes(y = Upper), linetype='dashed') + 
  geom_line(data=delta_t_data, aes(x=Time, y=Value), color='black', lwd=2)+
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = .5) + 
  theme_bw(base_size=20) + ylab("Value") + 
  ggtitle("A 95% Confidence Interval Example")
#############################################################################
################################################################################
###############################################################################
###############################################################################


stopCluster(cl)

sim_error_date_directory <- paste("../reports/", server_name,'/', "simulations","/" , Sys.Date(), "_simulations", sep='')

if (!dir.exists(sim_error_date_directory)){
  
  dir.create(path=sim_error_date_directory, recursive = TRUE)
  
}

file_name <- paste("report_", Sys.Date(), sep='')

files_in_corresponding_directory <- list.files(sim_error_date_directory)

sim_number <- length(files_in_corresponding_directory)+1

individual_simulation_directory <- paste(sim_error_date_directory, "/",file_name,
                                         "_run_", matrix_type, "_",sim_number, sep="")

dir.create(path=individual_simulation_directory)





write.csv(error_data, file=paste(individual_simulation_directory, "error.data.csv", sep='/'),
          row.names = FALSE)


## Could write as excel file as well
write.csv(params, file=paste(individual_simulation_directory, "parameters.csv", sep='/'),
          row.names = TRUE)

mu.plot <- (ggplot(NULL, aes(x= mus.chosen)) + geom_bar() + xlab("Mu")
            +ggtitle("Mus That Were Chosen") + theme_bw()
)


treatLengthStats <- data.frame(table(rowSums(W)))

treatLengthStats$Var1 <- as.numeric(as.character(treatLengthStats$Var1))

treatLengthPlot <- (ggplot(treatLengthStats, aes(x=Var1, y=Freq)) +
  geom_point() + 
    geom_linerange(aes(x=Var1, ymax=Freq, ymin=0.75)) 
  +theme_bw(base_size = 20) + xlab('Treatment Length') + ylab("Count"))


ggsave("muPlot.pdf", mu.plot, "pdf", path=paste(individual_simulation_directory, sep='/'), width=11, height=8.5)

ggsave("treatmentLengthPlot.pdf", treatLengthPlot, "pdf", path=paste(individual_simulation_directory, sep='/'),
      width=11, height=8.5)
                    
ggsave("ciExamplePlot.pdf", ciExamplePlot, "pdf", path=paste(individual_simulation_directory, sep='/'),
      width=11, height=8.5)

tabling_mus <- table(mus.chosen)

write.csv(tabling_mus, file=paste(individual_simulation_directory, "tabled.mus.csv", sep='/'), row.names=FALSE)

write.csv(estimated_coverage, 
          file=paste(individual_simulation_directory, "est_bootstrap_coverage.csv", sep='/'))

fileConn<-file(paste(individual_simulation_directory, "error table.txt", sep='/'))
writeLines(texErrorTable, fileConn)
close(fileConn)
                                        
                                        
fileConn<-file(paste(individual_simulation_directory, "est_bootstrap_coverage_latex.txt", sep='/'))
writeLines(coverageTable, fileConn)
close(fileConn)



#fileConn<-file(paste(individual_simulation_directory, "coverage table.txt", sep='/'))
#writeLines(covTableFixed, fileConn)
#close(fileConn)


