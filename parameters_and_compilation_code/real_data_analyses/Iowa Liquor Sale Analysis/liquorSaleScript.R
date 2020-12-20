library(ggplot2)
library(MASS)
library(usmap)
library(dplyr)
library(counties)
library(stringr)
library(reshape)
library(LaplacesDemon)
library(glmnet)

source('../../causal_inference_methods_code_corrupted.R')


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

seriesWithDrunkDrivingDecay <- function(aSeries, finalProp=.96, decayRate=.1, treatmentIndexer){ ### Creates treatment effect for implementation of drunk driving law
    
    if (length(aSeries) != length(treatmentIndexer)){
        
        stop("Series and Treatment Indexer Must be the Same Length.")
    }
    
    if (all(treatmentIndexer==0)){
        
        return(aSeries)
        
    }else{
    
        dataWhereTransformNeeded <- aSeries[treatmentIndexer==1]
    
        transformedPart <- finalProp*dataWhereTransformNeeded+(1-finalProp)*dataWhereTransformNeeded*exp(-1*decayRate*(1:length(dataWhereTransformNeeded)))
        
        aSeries[treatmentIndexer==1] <- transformedPart
        
        return(aSeries)
        
    }
    
}



saleData=read.csv("IowaLiquorSaleData.csv", stringsAsFactors=FALSE)
saleData = saleData %>% mutate(revenue_per_order = order_revenue/number_of_orders,
                               liquor_per_order = order_amount/number_of_orders,
                               year_month=order_year+round((order_month-1)/12, 3))

saleData  <- saleData  %>% filter(order_year < 2020)


allYearMonths <- unique(saleData[,c('order_year', "order_month")])

saleData$fips <- fips('iowa', saleData$county)



desired_clusters <- 30

max_available_clusters <- detectCores()

cl <- makeCluster(min(c(max_available_clusters, desired_clusters))-3)

registerDoParallel(cl)




#### Impute missing data for Adams county
expandedAmountPerOrderData <- variableExpander(saleData %>% select(c("county", 'year_month',
                                                                     "liquor_per_order")), 
                                     unitKey = 'county',
                timeKey='year_month')


WMissingDataAmountPerOrder <- is.na(expandedAmountPerOrderData)

WMissingDataAmountPerOrder <- WMissingDataAmountPerOrder+0


imputedAmountPerOrderData <- completion_with_rank_estimation_validate_mu(
        Y=expandedAmountPerOrderData, W=WMissingDataAmountPerOrder,  
        weight_matrix =NULL,
                                                               initial_rank=100,
                                                               tolerance=1e-04, 
                                                               validation_max_iter=5000,
                                                               min_iter=100,
                                                               max_iter=10000,
                                                               mu_grid=10^seq(-3,3, 1),
                                                               K=5)

imputedAmountPerOrder <- imputedAmountPerOrderData$L_hat






expandedRevenuePerOrderData <- variableExpander(saleData %>% select(c("county", 'year_month', "revenue_per_order")), 
                                      unitKey = 'county', timeKey='year_month')

WMissingDataRevenuePerOrder <- is.na(expandedRevenuePerOrderData)

WMissingDataRevenuePerOrder <- WMissingDataRevenuePerOrder+0


imputedRevenuePerOrderData <- completion_with_rank_estimation_validate_mu(
        Y=expandedRevenuePerOrderData, W=WMissingDataRevenuePerOrder,  
        weight_matrix =NULL,
                                                               initial_rank=100,
                                                               tolerance=1e-04, 
                                                               validation_max_iter=5000,
                                                               min_iter=100,
                                                               max_iter=10000,
                                                               mu_grid=10^seq(-3,3, 1),
                                                               K=5)

imputedRevenuePerOrder <- imputedRevenuePerOrderData$L_hat
 
 
 
 
expandedOrderCountData <- variableExpander(saleData %>% select(c("county", 'year_month', "number_of_orders")), 
                                      unitKey = 'county', timeKey='year_month')

WMissingDataOrderCount <- is.na(expandedOrderCountData)

WMissingDataOrderCount <- WMissingDataOrderCount+0


imputedOrderCountData <- completion_with_rank_estimation_validate_mu(
    Y=expandedOrderCountData, W=WMissingDataOrderCount,  
    weight_matrix =NULL,
                                                           initial_rank=100,
                                                           tolerance=1e-04, 
                                                           validation_max_iter=5000,
                                                           min_iter=100,
                                                           max_iter=10000,
                                                           mu_grid=10^seq(-3,3, 1),
                                                           K=5)

imputedOrderCount <- round(abs(imputedOrderCountData$L_hat), 0) 



completedOrderAmountData <- melt(imputedAmountPerOrder*imputedOrderCount)

names(completedOrderAmountData) <- c('county', 'year_month', 'order_amount')


completedRevenueData <- melt(imputedRevenuePerOrder*imputedOrderCount)

names(completedRevenueData) <- c('county', 'year_month', 'order_revenue')


completedOrderCountData <- melt(imputedOrderCount)

names(completedOrderCountData) <- c('county', 'year_month', 'number_of_orders')


saleData <- unique(saleData %>% select(c("county", "fips"))) %>% 
inner_join(completedOrderAmountData, by=c("county")) %>% 
inner_join(completedRevenueData) %>% 
inner_join(completedOrderCountData) %>% 
mutate(revenue_per_order=order_revenue/number_of_orders,
      amount_per_order=order_amount/number_of_orders)


iowaCountyData <- read.csv("IowaCountyData.csv", stringsAsFactors = F)
names(iowaCountyData) <- tolower(names(iowaCountyData))


iowaCountyData$areakm2 <- as.numeric(str_trim(str_replace(str_extract(iowaCountyData$area,
                        pattern ="(?<=\\().*(?=km2)"), pattern="\\,", replace="")))

 
iowaCountyData <- iowaCountyData %>% mutate(pop_density = population/areakm2,
                                           county = str_trim(str_replace(tolower(county), pattern='county',
                                                               replace=""))) 
 
saleDataWithCountyDemos <- saleData %>% left_join(iowaCountyData)



otherCensusData <- life %>% inner_join(demographics %>% select(-c("population"))) %>% inner_join(income) %>% 
filter(year==2012) %>% select(-c('year'))



saleDataWithCountyDemos <- saleDataWithCountyDemos %>% inner_join(otherCensusData, 
                                                                  by=c("fips"='county_fips') )


dataFrom2012 <- saleDataWithCountyDemos %>% filter(floor(year_month)==2012)


meanByCounty = saleData %>% group_by(county) %>% summarize(mean_order_amount=
                                                          mean(order_amount),
                                                          mean_order_revenue=mean(order_revenue))

meanByCounty$fips <- fips('iowa', meanByCounty$county)
############## SET PARAMETERS HERE
scoreBeta <- c(0, 1)

rateOfDecay <- .05
finalProp <- .6
Time0 <- 84


covariates <- saleDataWithCountyDemos %>% select(c("county", "revenue_per_order", 'amount_per_order')) %>% 
group_by(county) %>% summarize(mean_rev_per_order=log(mean(revenue_per_order, na.rm=T)), 
                               mean_amount_per_order=mean(amount_per_order, na.rm=T)) %>% 
select(mean_rev_per_order, mean_amount_per_order)


fullOutcomeData <- variableExpander(saleDataWithCountyDemos[,c("amount_per_order", 'county',
                                                           'year_month')], unitKey='county',
                                timeKey='year_month')

simulatedTreatmentData <- parameterMakerBlockProvidedCovariates(Time=dim(fullOutcomeData)[2],
                                                                modelCovariates=as.matrix(covariates), Time0=Time0, propTreat=.3,
                                     glmBeta=scoreBeta)

observedData<- t(sapply(1:dim(fullOutcomeData)[1], FUN=function(x) seriesWithDrunkDrivingDecay(aSeries=fullOutcomeData[x,], 
                                                                                          finalProp=.6, decayRate=.05, 
                                    treatmentIndexer = simulatedTreatmentData$W[x,])))

rownames(observedData) = rownames(fullOutcomeData)

treatDesignMatrix <- simulatedTreatmentData$W

rownames(treatDesignMatrix) <- rownames(fullOutcomeData)


adoptData <- melt(treatDesignMatrix) %>% group_by(X1) %>% summarize(status = c('Did not Adopt', 'Did Adopt')[max(value)+1]) %>% 
dplyr::rename(county=X1) %>% inner_join(meanByCounty)

adoptionMap <- plot_usmap(include='IA', region='counties',
           data = adoptData, values = "status") +
  scale_fill_manual(values = c("Did not Adopt" = "white", "Did Adopt" = "blue"), name = "Safe and Sober Ammendment Adoption") + 
  theme(legend.position = "right") + 
ggtitle("Iowa\'s Adoption of the Safe and Sober Ammendment, by County.")


                        
treatProfile <- colMeans(observedData[which(rowSums(simulatedTreatmentData$W) > 0),])
untreatProfile <- colMeans(observedData[which(rowSums(simulatedTreatmentData$W) == 0),])
neededTimes <- as.numeric(names(untreatProfile))
                        
                        
averageUnitData <- cbind.data.frame(rep(neededTimes, 2), c(unname(treatProfile), unname(untreatProfile)), rep(c('Treated', 'Untreated'), each=length(neededTimes)))

names(averageUnitData) <- c('Time', 'Log Amount per Order', 'Unit Type')                        
                        
treatAndUntreatedTogether <- ggplot(averageUnitData %>% filter(Time >= 2017), aes(x=Time, y=`Log Amount per Order`, col=`Unit Type`)) + geom_line(lwd=2) + 
theme_bw(base_size=20) + geom_vline(xintercept = averageUnitData$Time[Time0], lwd=1.5) + ggtitle("Iowa Liquor Data (2017-2019)")                        
                        
                        


set.seed(529371)

numSims <- 300

for (simNumber in 1:numSims){

    
    
    simulatedTreatmentData <- parameterMakerBlockProvidedCovariates(Time=dim(fullOutcomeData)[2],
                                                                modelCovariates=as.matrix(covariates), Time0=Time0, propTreat=.3,
                                     glmBeta=scoreBeta)
            
    if (simNumber ==1){
        
        effects <- array(NA, dim=c(5, dim(fullOutcomeData)[2]-Time0, numSims), dimnames=list(c('MC-NNM', 
                                            "Weighted softImpute", 'FACTOR', 'Weighted R1Comp', 'Truth'),
                                                                           paste('', 1:(dim(fullOutcomeData)[2]-Time0), sep=''),
                                                                           paste('iteration', 1:numSims, sep='_')))
        
    }


    observedData <- log(t(sapply(1:dim(fullOutcomeData)[1], FUN=function(x) seriesWithDrunkDrivingDecay(
        aSeries=fullOutcomeData[x,], finalProp=finalProp, decayRate=rateOfDecay,
        treatmentIndexer = simulatedTreatmentData$W[x,]))))
                             
    rownames(observedData) = rownames(observedData)                         
                             
    mc_nnm_info <- matrix_completion_causal(Y=observedData, W=simulatedTreatmentData$W, num_iter=1000, K=4, 
                                            lambda_grid=c(0, 10^seq(-4,2,1), seq(2,5,1)),
                                            tolerance=1e-04)

    weightedSoftImputeInfo <- weightedSoftImpute_validate_lambda(Y=observedData, W=simulatedTreatmentData$W, 
                                        weight_matrix=simulatedTreatmentData$weightMatrixEst, num_iter=1000,
                                        K=5, lambda_grid=seq(0, 2000, 100), tolerance=1e-03)

    L_completionFactorModel <- completion_factor_model(Y=observedData, W=simulatedTreatmentData$W, 
                                                       propScoreMat = simulatedTreatmentData$propScoreEst,
                                                       numFactors=rankMatrix(mc_nnm_info$L_hat)[1])

    R1CompInfo <- completion_with_rank_estimation_validate_mu(
                                      Y=observedData, W=simulatedTreatmentData$W , 
        weight_matrix=simulatedTreatmentData$trueWeightMatrix,
                                      initial_rank=40,
                                      tolerance=1e-04, 
                                      validation_max_iter=5000,
                                      min_iter=100,
                                      max_iter=10000,
                                      mu_grid=0,
                                      K=5)
                         
    effects['Weighted R1Comp', , simNumber]  <- treat.estimator(observedData, R1CompInfo$L_hat, simulatedTreatmentData$W)
    
    effects['FACTOR', , simNumber] <- treat.estimator(observedData, L_completionFactorModel, simulatedTreatmentData$W)
    
    effects['Weighted softImpute', , simNumber] <-  treat.estimator(observedData, weightedSoftImputeInfo$L_hat, 
                                                                    simulatedTreatmentData$W)

    effects['MC-NNM', , simNumber] <- treat.estimator(observedData, mc_nnm_info$L_hat, simulatedTreatmentData$W) 
    
    effects['Truth', , simNumber] <- treat.estimator(observedData, log(fullOutcomeData), simulatedTreatmentData$W) 
    
    
    }


meltedEffectData <- melt(effects[,,])

names(meltedEffectData) <- c('Method', 'Time', 'Iteration', 'Effect')
meltedEffectData$Method <- factor(meltedEffectData$Method, 
                                  levels=c("MC-NNM", "Weighted softImpute", 
                                           "FACTOR", "Weighted R1Comp", 'Truth'))


effectDataPaired <- meltedEffectData %>% inner_join(meltedEffectData %>% filter(Method=='Truth') %>% select(c("Time", "Iteration", "Effect"))
                                , by=c('Time', 'Iteration') ) %>% filter(Method != 'Truth')

names(effectDataPaired)[4:5] <- c('Effect', 'TrueEffect')



errorTable <- (effectDataPaired %>% group_by(Method) %>% summarize(
    MSE = round(mean((Effect-TrueEffect)^2), 2), SE = round(sd((Effect-TrueEffect)^2)/sqrt(n()), 2))
 %>% arrange(Method))



mseTableForLatex <- errorTable %>% mutate(`MSE (SE)` = paste(MSE, ' (', SE, ')', sep='')) %>% select(Method, `MSE (SE)`)

mseTableForLatex <- knitr::kable(mseTableForLatex, format='latex')




set.seed(529371)

numSims <- 300

bootSims <- 300

chosenLength <- NULL

for (simNumber in 1:numSims){
    
    startTime <- Sys.time()
    
    simulatedTreatmentData <- parameterMakerBlockProvidedCovariates(Time=dim(fullOutcomeData)[2],
                                                                modelCovariates=as.matrix(covariates), Time0=Time0, propTreat=.3,
                                     glmBeta=scoreBeta)

    if (simNumber ==1){
        
        bootEstData <- array(NA, dim=c(bootSims, dim(fullOutcomeData)[2]-Time0, numSims), 
                             dimnames=list(paste("boot_iteration", 1:bootSims, sep="_"),
                                                                           paste('', 1:(dim(fullOutcomeData)[2]-Time0), sep=''),
                                                                           paste('iteration', 1:numSims, sep='_')))
        
        trueEffects <- array(NA, dim=c(numSims, dim(fullOutcomeData)[2]-Time0), dimnames=list(paste('iteration', 1:numSims, sep='_'), 
                                                                             paste('', 1:(dim(fullOutcomeData)[2]-Time0), sep='')))
        
    }


    observedData <- log(t(sapply(1:dim(fullOutcomeData)[1], FUN=function(x) seriesWithDrunkDrivingDecay(
        aSeries=fullOutcomeData[x,], finalProp=finalProp, decayRate=rateOfDecay,
        treatmentIndexer = simulatedTreatmentData$W[x,]))))
    
    
    boot_ests_R1Comp <- bootstrapCI(Y=observedData, W=simulatedTreatmentData$W, 
    weightMatrix=simulatedTreatmentData$weightMatrixEst, bootstrap_samps=bootSims, 
    method=completion_with_rank_estimation_validate_mu,
    initial_rank=40,
                                  tolerance=1e-04, 
                                  validation_max_iter=5000,
                                  min_iter=100,
                                  max_iter=1000,
                                  mu_grid=0,
                                  K=5)
    
   # estimBiases <- biasSampler(Y=observedData, W=simulatedTreatmentData$W, 
    #    weight_matrix=simulatedTreatmentData$weightMatrixEst, numberOfBiasSamples=500)
    
    
   # biasSample <- estimBiases[sample(1:dim(estimBiases)[1], size=dim(boot_ests_R1Comp)[1], replace=T),]

    
   # boot_ests_R1Comp <- boot_ests_R1Comp+biasSample
    
    trueEffects[simNumber,] <- treat.estimator(observedData, log(fullOutcomeData), simulatedTreatmentData$W) 
          
    
    bootEstData[,,simNumber] <- boot_ests_R1Comp
                                 
    endTime <- Sys.time()
                                 
    print(simNumber)                            
    
    print(difftime(endTime, startTime, units='mins'))
    
    }


ci_array <- apply(bootEstData[, , ], MARGIN=c(2, 3), FUN=bs_percentile_method,
                       1-(1-.95)/1)




betweenTwoNumbers <- function(number, numberSandwich){

        return(numberSandwich[1] <= number & number <= numberSandwich[2])

    }

weDoGood <- matrix(NA, nrow=dim(ci_array)[3], ncol=dim(ci_array)[2])

for (iteration in 1:dim(bootEstData)[3]){
    
    whereWeHit <- c()
    
    for (colIndex in 1:dim(ci_array)[2]){

            theBounds <- ci_array[, colIndex, iteration]

            whereWeHit <- as.vector(c(whereWeHit, betweenTwoNumbers(trueEffects[iteration, colIndex], theBounds)))

        }
    
    weDoGood[iteration,] <- whereWeHit    
    
}





confidenceDataFrame <- data.frame(t(colMeans(weDoGood)))

rownames(confidenceDataFrame)[1] <- 'Estimated Coverage'

names(confidenceDataFrame) <- paste("Time", 1:dim(weDoGood)[2], sep=' ')

confidenceDataFrameLatex <- knitr::kable(confidenceDataFrame, format='latex')

chosenIteration <- which.min(which(rowSums(weDoGood) == max(rowSums(weDoGood))))

treatTimes <- 1:dim(ci_array)[2]

chosenIterationData <- ci_array[,,chosenIteration ]

chosenIterationTruth <- trueEffects[chosenIteration, ]

chosenTruthData <- cbind.data.frame(treatTimes, chosenIterationTruth)

colnames(chosenTruthData) <- c("Time", "Value")

bootEstimateData <- cbind.data.frame(treatTimes, colMeans(chosenIterationData))

colnames(bootEstimateData) <- c("Time", "Value")

ChosenCi <- apply(chosenIterationData, MARGIN=2, FUN=bs_percentile_method ,1-(1-.95)/length(chosenIterationTruth))


names(ChosenCi) <- treatTimes


ciExample <- cbind.data.frame(treatTimes, t(ChosenCi))
names(ciExample) <- c("Time", "Lower", "Upper")
                    
ciExamplePlot <- ggplot(ciExample, aes(x=factor(Time), y = Upper, group=1)) +
  geom_line(aes(y = Lower), linetype='dashed') + 
  geom_line(aes(y = Upper), linetype='dashed') + 
  geom_line(data=chosenTruthData , aes(x=factor(Time), y=Value), color='black', lwd=2)+

  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = .5) + 
  theme_bw(base_size=20)+ xlab("Time") + ylab("Difference in log Average Amount Per order (Liters)") + 
  ggtitle("Estimated Effect of Safe and Sober Law,\nwith a 95% Bootstrapped Confidence Interval")




if (!exists("./plotsAndTables")){
    
    dir.create('./plotsAndTables')
    
}

outputDir <- './plotsAndTables'




makeFilename <- function(directory, fileName){
    
    paste(directory, fileName, sep='/')
}

ggsave(makeFilename(outputDir, "ci_example_alcohol.pdf"), ciExamplePlot, width=11, height=8.5)

ggsave(makeFilename(outputDir, "adoption_map.pdf"), adoptionMap, width=11, height=8.5)
                                 
ggsave(makeFilename(outputDir, "treated_and_untreated.pdf"), treatAndUntreatedTogether, width=11, height=8.5)                                
                                 


fileConn<-file(makeFilename(outputDir, "effect_estimate_evaluation_alcohol.txt"))
writeLines(mseTableForLatex, fileConn)
close(fileConn)


fileConn<-file(makeFilename(outputDir, "estimated_coverage_alcohol.txt"))
writeLines(confidenceDataFrameLatex , fileConn)
close(fileConn)




stopCluster(cl)


