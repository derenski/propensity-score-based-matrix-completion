library(plyr)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(glmnet)


variableExpander <- function(keysAndVariable, unitKey, timeKey){ ### Expands into unit/time format
  
  orderedData <- keysAndVariable %>% arrange(get(unitKey), get(timeKey))

  outcomeName <- names(keysAndVariable)[!(names(keysAndVariable) %in% c(unitKey, timeKey))]
  
  outcomeMatrixForm <- matrix(NA, nrow=length(unique(keysAndVariable[, unitKey])), 
                              ncol=length(unique(keysAndVariable[, timeKey])), byrow=T)
  
  rownames(outcomeMatrixForm) <- unique(keysAndVariable[, unitKey])
  
  colnames(outcomeMatrixForm) <- str_replace_all(unique(keysAndVariable[, timeKey])[order(unique(keysAndVariable[, timeKey]))], pattern="-", replace="")
  
  for (index in 1:length(keysAndVariable[, unitKey])){
    
    outcomeMatrixForm[keysAndVariable[, unitKey][index], str_replace_all(keysAndVariable[, timeKey][index], pattern="-", replace="")] <- keysAndVariable[, outcomeName][index]
    
  }
  
  return(outcomeMatrixForm)
  
}


source('../causal_inference_methods_code_corrupted.R')

options(stringsAsFactors = F)

covidData <- read.csv("/home/josh/GitHub/covid-19-data/us-states.csv")

covidData <- covidData %>% filter(state %in% state.name) %>% mutate(date=ymd(date)) %>% arrange(state, date)

everyonesStartDate <- covidData %>% dplyr::group_by(state) %>% dplyr::summarize(startDate=min(date))






partyData <- read.csv("../../../California Smoking Data/political_party_data.csv")

popData <- read.delim("../../../California Smoking Data/Population_Data_Annual.txt")

lockdownData <- read.csv("../../../California Smoking Data/Covid Response Data/lockdown data.csv")

lockdownData <- data.frame(apply(lockdownData, MARGIN=2, FUN=function(x) str_replace(x, pattern="\\(.*\\)", replace="")))

openUpData <- read.csv("../../../California Smoking Data/Covid Response Data/lockdown lifting data.csv")

lockdownData <- lockdownData %>% mutate(State.of.emergency.declared=mdy(paste(State.of.emergency.declared, ", 2020", sep="")),
                                        Stay.at.home.ordered=mdy(paste(Stay.at.home.ordered, ", 2020", sep="") ))

underlockVec <- c()

for (i in 1:dim(covidData)[1]){
  
  underlockVec <- c(underlockVec, covidData[i,'date'] >= 
                      (lockdownData %>% filter(State.territory == covidData[i,'state']))$Stay.at.home.ordered)
  
  underlockVec[is.na(underlockVec)] <- 0
  
}


covidData <- covidData %>% mutate(under_lockdown = underlockVec)

openUpData <- openUpData %>% mutate(Date.enacted = mdy(Date.enacted), Date.lifted=mdy(Date.lifted))

relevantCovidData <- covidData %>% filter(date < min(openUpData$Date.lifted))


relevantPartyData <- partyData %>% filter(year==2019)

relevantPartyData <- relevantPartyData[, -1]

nebraska2019Data <- list(2019, "Nebraska", 0, 1, 18, 49, 18+49, 0, 0, 0, 0, 3, 3)

hawaii2019Data <-  list(2019, "Hawaii", 1, 0, 46, 5, 46+5, 24, 1, 25, 2, 0, 2)

relevantPartyData <- rbind.data.frame(relevantPartyData, nebraska2019Data, hawaii2019Data)

meltedPopData <- reshape2::melt(popData, id.vars='DATE')

names(meltedPopData) <- c("year", "abb", "population")

meltedPopData <- meltedPopData %>% mutate(year = year(ymd(year)), abb = str_replace(abb, "POP$", ""),
                                          population=population*1000)

statesWeNeed = cbind(state.abb, state.area, 
                     as.character(state.region), as.character(state.division))

nameAreaTable <- data.frame(cbind(state.name, statesWeNeed))

names(nameAreaTable) <- c('state', 'abb', "area", 'region', 'division')

nameAreaTable$area <- as.numeric(nameAreaTable$area)

popAreaTable <- (meltedPopData %>% inner_join(nameAreaTable, by = "abb") %>% mutate(pop_density = population/area))

relevantPopAreaData <- popAreaTable %>% filter(year==2019)

relevantPopAreaData <- relevantPopAreaData[, !(names(relevantPopAreaData) %in% "year")]



covidPopPartyData <- (relevantCovidData %>% inner_join(relevantPopAreaData, by='state') %>% 
                        inner_join(relevantPartyData, by='state'))

W <- variableExpander(keysAndVariable=covidPopPartyData[,c('state', 'date', 'under_lockdown')], unitKey = 'state', timeKey = 'date')

W[is.na(W)] <- 0


Y <- variableExpander(keysAndVariable=covidPopPartyData[,c('state', 'date', 'cases')], unitKey = 'state', timeKey = 'date')

Y[is.na(Y)] <- 0



##### Modeling Adoption of stay-at-home order

stayAtHomeDuration <- rowSums(W)

stayHomeCovariates <- covidPopPartyData[!duplicated(covidPopPartyData$state),c('state', "pop_density", "stateHouseDemocrats", "region")]


stayHomeCovariates <- stayHomeCovariates %>% mutate(order_duration=stayAtHomeDuration )

mod <- glm(order_duration~., data=stayHomeCovariates[,!(names(stayHomeCovariates) %in% 'state') ], family=poisson())

probAdoptEst <- t(sapply(mod$fitted.values, FUN=function(x) rev(dpois(seq(0, dim(Y)[2], 1), lambda=x))))[,1:dim(Y)[2]]

propScoreEst <- t(apply(probAdoptEst, MARGIN=1, FUN=cumsum))

ourWeightMatrix <- propensity.score.to.weight(propScoreEst, W)

ourWeightMatrix[ourWeightMatrix > 5] <- 5


mc_nnm_info <- matrix_completion_causal(Y=Y, W=W , num_iter=1000, K=4, 
                                        lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                                        tolerance=1e-04)

L_mc_nnm <- mc_nnm_info$L_hat

weightedSoftImputeInfo <- weightedSoftImpute_validate_lambda(Y, W=W , weightMat=ourWeightMatrix/sum(ourWeightMatrix), num_iter=1000, K=5, 
                                                             lambda_grid=seq(0, 2000, 100), tolerance=1e-03)

L_weightedSoftImpute <- weightedSoftImputeInfo$L_hat

L_completionFactorModel <- completion_factor_model(Y=Y, W=W, propScoreMat = propScoreEst,
                                                   numFactors=rankMatrix(mc_nnm_info$L_hat)[1])

r1comp_info <- completion_with_rank_estimation_validate_mu(Y=Y, W=W,  
                                                           weight_matrix = ourWeightMatrix/sum(ourWeightMatrix),
                                                           initial_rank=24,
                                                           tolerance=1e-04, 
                                                           validation_max_iter=500,
                                                           min_iter=100,
                                                           max_iter=1000,
                                                           mu_grid=0,
                                                           K=5)

effectEstMCNNM <- treat.estimator(Y=Y, L.hat=L_mc_nnm , W=W)

effectEstWeightedSoftImpute <- treat.estimator(Y=Y, L.hat=L_weightedSoftImpute, W=W)

effectEstFactorModel <-  treat.estimator(Y=Y, L.hat=L_completionFactorModel, W=W)

effectEstR1Comp <-  treat.estimator(Y=Y, L.hat=r1comp_info$L_hat, W=W)

plot(effectEstFactorModel[1:10])

treatLengthDist <- rowSums(W)

numberHaveEachYear <- colSums(W)





