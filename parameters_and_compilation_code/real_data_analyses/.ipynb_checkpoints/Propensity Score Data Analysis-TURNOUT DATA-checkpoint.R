library(plyr)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(glmnet)
library(LaplacesDemon)
library(knitr)
library(usmap)
## setwd("./Real Data Analyses")

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

### "smoking_data.csv"

turnoutData <- read.csv("../../../California Smoking Data/turnout.csv")

adoptionYear_mailin <- turnoutData %>% group_by(abb) %>% dplyr::summarize(thing = year[min(which(policy_mail_in==1))])






### Smoking tax implemented in 1988

#turnoutData <- read.csv("../../California Smoking Data/smoking_data.csv")

partyData <- read.csv("../../../California Smoking Data/political_party_data.csv")

partyData <- partyData %>% filter(!(state %in% c("North Dakota", "Nebraska", "Hawaii", "Alaska", "Kentucky") ))


popData <- read.delim("../../../California Smoking Data/Population_Data_Annual.txt")

meltedPopData <- reshape2::melt(popData, id.vars='DATE')

names(meltedPopData) <- c("year", "abb", "population")

meltedPopData <- meltedPopData %>% mutate(year = year(ymd(year)), abb = str_replace(abb, "POP$", ""),
                                          population=population*1000)

badStates <- c("HI", "AK", "NE", "ND", "KY")

statesWeNeed = cbind(state.abb[!(state.abb %in% badStates)], state.area[!(state.abb %in% badStates)], 
                     as.character(state.region)[!(state.abb %in% badStates)], as.character(state.division)[!(state.abb %in% badStates)])

nameAreaTable <- data.frame(cbind(unique(partyData$state), statesWeNeed))

names(nameAreaTable) <- c('state', 'abb', "area", 'region', 'division')

nameAreaTable$area <- as.numeric(nameAreaTable$area)


popAreaTable <- (meltedPopData %>% inner_join(nameAreaTable, by = "abb") %>% mutate(pop_density = population/area))

data2015 <- popAreaTable %>% filter(year==2015)

### Filtered Party data


## TO identify missingsness, look for where total is zero. 



neededYears <- unique(turnoutData$year)

filteredPartyData <- unique(partyData %>% filter(year %in% neededYears))

keys <- c('year', 'state')

X <- as.data.table(filteredPartyData)

filteredPartyDataUniqued <- X[,lapply(.SD,mean),keys]

##### Connecting population density to political party data

partyPopData <- filteredPartyDataUniqued %>% inner_join(popAreaTable, by=c('year', 'state'))

byRegion <- partyPopData %>% group_by(region, division, year) %>% dplyr::summarise(propStateHouseDem = sum(stateHouseDemocrats)/sum(stateHouseTotal))


# "North Central" "Northeast"     "South"         "West" 

regOfInterest <- byRegion %>% filter(region=='West')

bigDemPlot <- (ggplot(byRegion, aes(x=year, y=propStateHouseDem, col=division)) + geom_line() + theme_bw(base_size=20)
+ facet_wrap(~ region, scales='free') + ggtitle("Democrat Representation in State Houses")
  +ylab("Proportion of Representatives who are Democrat") + xlab("Year"))
  
DivOfInterest <- 'Mountain'
west <- byRegion %>% filter(division == DivOfInterest)
aPlot <- (ggplot(west, aes(x=year, y=propStateHouseDem)) + geom_line() + theme_bw(base_size=20)
  + ggtitle(paste("Democrat Power in State Houses\n","(",DivOfInterest," Division)", sep='')) + ylab("Proportion of Reps who are Democrat") + xlab("Year"))

#ggsave(paste("Democrats in", DivOfInterest, "State Houses.pdf"), plot=aPlot, width=11, height=8.5)
  
(ggplot(partyPopData, aes(x=stateHouseDemocrats, y=stateHouseRepublicans)) + geom_point() + theme_bw(base_size=20)
  )





### Modeling the adoption of policies 

modelingDataWithYear <- ((partyPopData %>% inner_join(turnoutData) %>% mutate(stateHousePropDemocrat = stateHouseDemocrats/stateHouseTotal,
                                                                            stateSenatePropDemocrat = stateSenateDemocrats/stateHouseTotal,
                                                                            usHousePropDemocrat = usHouseDemocrats/usHouseTotal))           
)

#### Imputing Minnesota


##### State House Republican Data

modelingDataWithYear$stateHouseRepublicans[modelingDataWithYear$stateHouseTotal==0] <- NA

missingStateHouseRepublicanData <- variableExpander(modelingDataWithYear[,c("abb", "year", "stateHouseRepublicans")], unitKey="abb", 
                 timeKey="year")





WmissingStateHouseRepublicanData <- array(0, dim=dim(missingStateHouseRepublicanData))

WmissingStateHouseRepublicanData[is.na(missingStateHouseRepublicanData)] <- 1

r1comp_infoMissingStateHouseRepublicanData <- completion_with_rank_estimation_validate_mu(Y=missingStateHouseRepublicanData, 
                                                                                          W=WmissingStateHouseRepublicanData,
                                                           initial_rank=24,
                                                           tolerance=1e-04, 
                                                           validation_max_iter=500,
                                                           min_iter=100,
                                                           max_iter=1000,
                                                           mu_grid=0,
                                                           K=5)

missingStateHouseRepublicanData["MN",]

### data starts at 1952

imputedMinnesotaStateHouseRepublicans <- floor(r1comp_infoMissingStateHouseRepublicanData$L_hat["MN",as.character(
  seq(1920, 1952,4))])



modelingDataWithYear$stateHouseRepublicans[is.na(modelingDataWithYear$stateHouseRepublicans)] <- imputedMinnesotaStateHouseRepublicans




#### State House Democrat Data


modelingDataWithYear$stateHouseDemocrats[modelingDataWithYear$stateHouseTotal==0] <- NA

missingstateHouseDemocratData <- variableExpander(modelingDataWithYear[,c("abb", "year", "stateHouseDemocrats")], unitKey="abb", 
                                                    timeKey="year")

WmissingstateHouseDemocratData <- array(0, dim=dim(missingstateHouseDemocratData))

WmissingstateHouseDemocratData[is.na(missingstateHouseDemocratData)] <- 1

r1comp_infoMissingstateHouseDemocratData <- completion_with_rank_estimation_validate_mu(Y=missingstateHouseDemocratData, 
                                                                                        W=WmissingstateHouseDemocratData,
                                                                                        initial_rank=24,
                                                                                        tolerance=1e-04, 
                                                                                        validation_max_iter=500,
                                                                                        min_iter=100,
                                                                                        max_iter=1000,
                                                                                        mu_grid=0,
                                                                                        K=5)

### data starts at 1952

imputedMinnesotastateHouseDemocrats <- floor(r1comp_infoMissingstateHouseDemocratData$L_hat["MN",as.character(
  seq(1920, 1952,4))])



modelingDataWithYear$stateHouseDemocrats[is.na(modelingDataWithYear$stateHouseDemocrats)] <- imputedMinnesotastateHouseDemocrats




modelingDataWithYear$stateHouseTotal[(modelingDataWithYear$state=="Minnesota") & modelingDataWithYear$stateHouseTotal==0] <- (
  modelingDataWithYear$stateHouseRepublicans[(modelingDataWithYear$state=="Minnesota") & modelingDataWithYear$stateHouseTotal==0] + 
    modelingDataWithYear$stateHouseDemocrats[(modelingDataWithYear$state=="Minnesota") & modelingDataWithYear$stateHouseTotal==0])








##### State Senate Republican Data

modelingDataWithYear$stateSenateRepublicans[modelingDataWithYear$stateSenateTotal==0] <- NA

missingStateSenateRepublicanData <- variableExpander(modelingDataWithYear[,c("abb", "year", "stateSenateRepublicans")], unitKey="abb", 
                                                    timeKey="year")





WmissingStateSenateRepublicanData <- array(0, dim=dim(missingStateSenateRepublicanData))

WmissingStateSenateRepublicanData[is.na(missingStateSenateRepublicanData)] <- 1

r1comp_infoMissingStateSenateRepublicanData <- completion_with_rank_estimation_validate_mu(Y=missingStateSenateRepublicanData, 
                                                                                          W=WmissingStateSenateRepublicanData,
                                                                                          initial_rank=24,
                                                                                          tolerance=1e-04, 
                                                                                          validation_max_iter=500,
                                                                                          min_iter=100,
                                                                                          max_iter=1000,
                                                                                          mu_grid=0,
                                                                                          K=5)

### data starts at 1952

imputedMinnesotaStateSenateRepublicans <- floor(r1comp_infoMissingStateSenateRepublicanData$L_hat["MN",as.character(
  seq(1920, 1952,4))])



modelingDataWithYear$stateSenateRepublicans[is.na(modelingDataWithYear$stateSenateRepublicans)] <- imputedMinnesotaStateSenateRepublicans




#### State Senate Democrat Data


modelingDataWithYear$stateSenateDemocrats[modelingDataWithYear$stateSenateTotal==0] <- NA

missingstateSenateDemocratData <- variableExpander(modelingDataWithYear[,c("abb", "year", "stateSenateDemocrats")], unitKey="abb", 
                                                  timeKey="year")


WmissingstateSenateDemocratData <- array(0, dim=dim(missingstateSenateDemocratData))

WmissingstateSenateDemocratData[is.na(missingstateSenateDemocratData)] <- 1

r1comp_infoMissingstateSenateDemocratData <- completion_with_rank_estimation_validate_mu(Y=missingstateSenateDemocratData, 
                                                                                          W=WmissingstateSenateDemocratData,
                                                                                          initial_rank=24,
                                                                                          tolerance=1e-04, 
                                                                                          validation_max_iter=500,
                                                                                          min_iter=100,
                                                                                          max_iter=1000,
                                                                                          mu_grid=0,
                                                                                          K=5)

### data starts at 1952

imputedMinnesotastateSenateDemocrats <- floor(r1comp_infoMissingstateSenateDemocratData$L_hat["MN",as.character(
  seq(1920, 1952,4))])



modelingDataWithYear$stateSenateDemocrats[is.na(modelingDataWithYear$stateSenateDemocrats)] <- imputedMinnesotastateSenateDemocrats




stateLocationData <- unique(modelingDataWithYear[,c("abb", "region", "division")])

modelingDataWithYear$stateSenateTotal[(modelingDataWithYear$state=="Minnesota") & modelingDataWithYear$stateSenateTotal==0] <- (
      modelingDataWithYear$stateSenateRepublicans[(modelingDataWithYear$state=="Minnesota") & modelingDataWithYear$stateSenateTotal==0] + 
    modelingDataWithYear$stateSenateDemocrats[(modelingDataWithYear$state=="Minnesota") & modelingDataWithYear$stateSenateTotal==0])









modelingDataNumericOnly <- modelingDataWithYear[, !(names(modelingDataWithYear) %in% c('state', 'abb', 
        'division', 'region', 'population', 'turnout', 'stateHouseRepublicans', 'stateSenateRepublicans', 'usHouseRepublicans',
        'republicanGovernor'))]


## Make treatment matrices for each of the voting policies

treatmentDurations <-(modelingDataWithYear %>% group_by(abb) %>% dplyr::summarize(duration_edr=sum(policy_edr), duration_mail_in=sum(policy_mail_in),
                                                                          duration_motor=sum(policy_motor)))


nAndTData = table(modelingDataWithYear$year)


neededYears <- as.numeric(names(nAndTData))

Ws <- alply(as.matrix(treatmentDurations[,-1]), 2, W_maker, N=as.numeric(nAndTData[1]), Time=length(nAndTData))

names(Ws) = names(treatmentDurations[,-1])

rownames(Ws$duration_edr) <- treatmentDurations$abb

colnames(Ws$duration_edr) <- names(nAndTData)

rownames(Ws$duration_mail_in) <- treatmentDurations$abb

colnames(Ws$duration_mail_in) <- names(nAndTData)

rownames(Ws$duration_motor) <- treatmentDurations$abb

colnames(Ws$duration_motor) <- names(nAndTData)

### Model this better

propScores <- matrix(NA, nrow=as.numeric(nAndTData[1]), ncol=length(nAndTData), byrow=T)

timeIndex <- 1

W <- Ws$duration_mail_in

if (F){

  for (theYear in unique(modelingDataNumericOnly$year)){
    
    thisData <- modelingDataNumericOnly %>% filter(year==theYear)
    
    model <- glm(policy_mail_in~policy_edr+policy_motor + stateHousePropDemocrat+pop_density
                 , data=thisData[, -1] , family=binomial())
    
    propScores[,timeIndex] <- model$fitted.values
    
    timeIndex <- timeIndex+1
    
  }
  
}

#### Alternative Modeling Strategy

areaUnderTimeSeries <- function(series){ ### Assumes equally spaced points
  
  individualAreas <- c()
  
  for (val in 1:(length(series)-1)){
    
    individualAreas <- c(individualAreas, .5*(series[val]+series[val+1]))
    
  }
  
  return(sum(individualAreas))
  
}

Y <- matrix(modelingDataWithYear$turnout, nrow=as.numeric(nAndTData[1]), ncol=length(nAndTData), byrow=T)

dataConstantVariablesOnly <- modelingDataWithYear %>% group_by(state) %>% dplyr::summarize(state_house_demo_score = areaUnderTimeSeries(stateHouseDemocrats/stateHouseTotal),
                                                                                           state_senate_demo_score = areaUnderTimeSeries(stateSenateDemocrats/stateSenateTotal),
                                                                                           prop_dem_gov = sum(democratGovernor)/length(democratGovernor),
    pop_density_score = mean(pop_density))


dataConstantVariablesOnly <- cbind(dataConstantVariablesOnly, treatmentDurations)

W <- Ws$duration_mail_in

mod <- glm(duration_mail_in~., data=dataConstantVariablesOnly[,!(names(dataConstantVariablesOnly) %in% 
                                                                 c('abb', 'state') )], 
           family=poisson())

probAdoptEst <- t(sapply(mod$fitted.values, FUN=function(x) rev(dpois(seq(0, dim(Y)[2], 1), lambda=x))))[,1:dim(Y)[2]]

propScores <- t(apply(probAdoptEst, MARGIN=1, FUN=cumsum))

ourWeightMatrix <- propensity.score.to.weight(propScores, W)

ourWeightMatrix <- propensity.score.to.weight(propScores , Ws$duration_mail_in)

# ourWeightMatrix[ourWeightMatrix  > 30] <- 30

desired_clusters <- 30

max_available_clusters <- detectCores()

cl <- makeCluster(min(c(max_available_clusters, desired_clusters))-3)

registerDoParallel(cl)

set.seed(5894532)


mc_nnm_info <- matrix_completion_causal(Y=Y, W=W , num_iter=1000, K=4, 
                                        lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                                        tolerance=1e-04)

L_mc_nnm <- mc_nnm_info$L_hat

weightedSoftImputeInfo <- weightedSoftImpute_validate_lambda(Y, W=W , weightMat=ourWeightMatrix/sum(ourWeightMatrix), num_iter=1000, K=5, 
                                                             lambda_grid=seq(0, 2000, 100), tolerance=1e-03)

L_weightedSoftImpute <- weightedSoftImputeInfo$L_hat

L_completionFactorModel <- completion_factor_model(Y=Y, W=W, propScoreMat = propScores,
                                                   numFactors=rankMatrix(mc_nnm_info$L_hat)[1])

r1comp_info <- completion_with_rank_estimation_validate_mu(Y=Y, W=W,  
                                                           weight_matrix = ourWeightMatrix,
                                                           initial_rank=40,
                                                           tolerance=1e-04, 
                                                           validation_max_iter=5000,
                                                           min_iter=100,
                                                           max_iter=10000,
                                                           mu_grid=0,
                                                           K=5)

effectEstMCNNM <- treat.estimator(Y=Y, L.hat=L_mc_nnm , W=W)

effectEstWeightedSoftImpute <- treat.estimator(Y=Y, L.hat=L_weightedSoftImpute, W=W)

effectEstFactorModel <-  treat.estimator(Y=Y, L.hat=L_completionFactorModel, W=W)

effectEstR1Comp <-  treat.estimator(Y=Y, L.hat=r1comp_info$L_hat, W=W)

stateLevelMeanIncreaseData <- c()

sum(W)/prod(dim(W))

for (aState in rownames(W)){
  
  rowOfInterest <- which(rownames(W)==aState)
  
  newW <- W
  
  newW[-rowOfInterest,] <- 0
  
  theEffect <- treat.estimator(Y=Y, L.hat=r1comp_info$L_hat, W=newW)
  
  stateLevelMeanIncreaseData <- rbind.data.frame(stateLevelMeanIncreaseData, 
                               cbind.data.frame(aState, mean(theEffect), sd(theEffect)/length(theEffect) ))
  
}

names(stateLevelMeanIncreaseData) <- c("abbr", 'Increase', "se_increase")

stateLevelMeanIncreaseData$Increase[is.nan(stateLevelMeanIncreaseData$Increase)] <- NA


fipsAndAbb <- statepop[, c("fips", "abbr")]

geoCodedIncreaseData <- (fipsAndAbb %>% right_join(stateLevelMeanIncreaseData, by='abbr')
                         %>% inner_join(stateLocationData, by=c('abbr'='abb')))

geoCodedIncreaseData %>% group_by(division) %>% dplyr::summarize(`Mean Increase` = mean(Increase, na.rm=T), `SD Increase`=
                                                                   sd(Increase, na.rm=T)/sum(!is.na(Increase)))

#### Visualizations for increase/decrease by state

us_map <- plot_usmap(data = geoCodedIncreaseData, values = "Increase", color = "red") + 
  scale_fill_continuous(name = "Increase (%)", label = scales::comma) +
  theme(legend.position = "right") + ggtitle("Turnout Increase due to Mail In Voting Across the US")



allEffectsTogether <- cbind(effectEstMCNNM, effectEstWeightedSoftImpute, effectEstFactorModel, effectEstR1Comp)[1:10,]

meanIncreases <- array(apply(allEffectsTogether, MARGIN=2, FUN= mean), dim=c(1, 4))

meanIncreaseTable <- rbind(c("NC-NNM", "Weighted Soft Impute", "Factor Model Estimator", "Weighted Rank-One Completion"),
                           meanIncreases)

rownames(meanIncreaseTable) <- c("Method", "Mean Turnout Increase (%)")

turnoutIncreaseTable <- kable(meanIncreases , 'latex')

increasePlotR1Comp <- (ggplot(NULL, aes(x=seq(1, 10, 1), effectEstR1Comp [1:10])) + geom_line() 
  + theme_bw(base_size=20) + xlab("General Elections Since Implementation") 
  + ylab("Average Increase in Turnout (%)") + ggtitle("Turnout Since Vote-by-Mail Implementation") + 
    scale_x_continuous(breaks = seq(1,10,1)))


turnoutSpecificAreas <- modelingDataWithYear %>%group_by(year, region, division) %>% dplyr::summarise(mean_turnout = mean(turnout))

big_turnout_plot <- (ggplot(turnoutSpecificAreas, aes(x=year, y=mean_turnout, col=division)) 
  + geom_line() + facet_wrap(~ region, scales='free') + theme_bw(base_size=20))



estimatedEffect <- stateLevelMeanIncreaseData %>% filter( abbr %in% c("CA", "CO", "OR"))

estimatedEffect <- (estimatedEffect %>% 
                      mutate(`Estimated Mean Increase (SE)`=paste(round(Increase, 3), 
                                                                  "% (",round(se_increase, 3), "%)", sep="") ))

estimatedEffect$`Literature Increase (SE)` <- c("-2.7% (.37%)", "10.938% (11.35%)", "4.4% (7%)")

compTable <- kable(estimatedEffect %>% dplyr::select(abbr, `Estimated Mean Increase (SE)`, `Literature Increase (SE)`),
                   "latex")

names(compTable)[1] <- 'State'



##### Summaries of adoption structure

adoptionYears <- as.numeric(colnames(W)[apply(W, MARGIN = 1, FUN=function(x) min(which(x==1)))])

numberEachYear <- table(adoptionYears)

allYears <- colnames(W)

numberEachYear[allYears[!(allYears %in% names(numberEachYear))]] <- 0

numberEachYear <- numberEachYear[order(as.numeric(names(numberEachYear)))]

adoptionTs <- ts(as.numeric(numberEachYear), frequency=.25, start=1920)

VBMAdoptionPlot <- (ggplot(NULL, aes(x=time(adoptionTs), y=adoptionTs)) + geom_line() + theme_bw(base_size=20) +
xlab("Year") + ylab("Number of States") + ggtitle("Vote By Mail Adoption over Time"))



bootstrap_samps <- 300

boot_ests <- matrix(NA, nrow=bootstrap_samps, ncol=length(effectEstR1Comp))

set.seed(734)


args_r1Comp <- list(
  initial_rank=40,
  tolerance=1e-04, 
  validation_max_iter=500,
  min_iter=10,
  max_iter=1000,
  mu_grid=0,
  K=5)

boot_ests_r1Comp <- bootstrapCI(Y=Y, W=W, ourWeightMatrix=ourWeightMatrix, bootstrap_samps=300, 
                         method=completion_with_rank_estimation_validate_mu, method_args = args_r1Comp )


Cis <- apply(rowMeans(boot_ests_r1Comp[,1:10]), MARGIN=2, FUN=bs_percentile_method ,.95)


Cis

bs_percentile_method(rowMeans(boot_ests), .95)


truehist(rowMeans(boot_ests))
abline(v=mean(effectEstR1Comp))


outputDirectory <- "./plots_and_tables"

if (!dir.exists(outputDirectory)){
  
  dir.create(outputDirectory)
  
}

ggsave(plot=increasePlotR1Comp, filename=paste(outputDirectory, "turnout increase r1comp.pdf", sep='/'), 
       units='in', width=8, height=8)

ggsave(plot=us_map, filename=paste(outputDirectory, "Average Increase by State.pdf", sep='/'), 
       units='in', width=8, height=8)

ggsave(VBMAdoptionPlot, filename=paste(outputDirectory, "VBM Adoption.pdf", sep='/'), 
       units='in', width=8, height=8)

ggsave(paste("Democrat Representation in State Houses.pdf"), plot=bigDemPlot, width=11, height=8.5)

fileConn<-file(paste(outputDirectory, "increase table.txt", sep='/'))
writeLines(turnoutIncreaseTable, fileConn)
close(fileConn)

fileConn<-file(paste(outputDirectory, "State comparison table.txt", sep='/'))
writeLines(compTable, fileConn)
close(fileConn)




