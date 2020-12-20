library(plyr)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(glmnet)

source('../causal_inference_methods_code_corrupted.R')

options(stringsAsFactors = F)

### California Adopts plicy in 1988


smokingData <- read.csv("../../../California Smoking Data/smoking_data.csv")


partyData <- read.csv("../../../California Smoking Data/political_party_data.csv")


combinedData <- smokingData %>% inner_join(partyData)













