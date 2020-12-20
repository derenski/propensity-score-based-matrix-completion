library(plyr)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(glmnet)
library(reticulate)

use_python("/usr/bin/python3", required = T)

server_name <- Sys.info()['nodename']

if (server_name == "Compy386"){

  source_python("/home/josh/GitHub/mnar_mc/mc_algorithms.py")

}else{

  source_python("/home/derenski/Github/mnar_mc/mc_algorithms.py")
  
}

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

### California Adopts plicy in 1988


smokingData <- read.csv("../../../California Smoking Data/smoking_data.csv",
                       row.names=1)


partyData <- read.csv("../../../California Smoking Data/political_party_data.csv")


combinedData <- smokingData %>% inner_join(partyData)



longSmokingData = melt(smokingData, id.vars=c('state', 'year'))

justSales <- longSmokingData %>% filter(variable=="cigsale")

justSales <- justSales[, -3]



smokingY <- variableExpander(justSales, unitKey="state", timeKey="year")

