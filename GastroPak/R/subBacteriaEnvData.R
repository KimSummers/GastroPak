# subBacteriaEnvData
#
# Get data for bacteria type, including mean and sd values
#
# file subBacteriaData
#
# inputs
# 	countData     - Count data to get sub data from
#   bacteriaType  - Include data with this bacteria type
#   sampleTypes   - Type of sample to get data for

# Version    Author       Date      Affiliation
# 1.00       J K Summers  10/01/24  Wellington Lab - School of Life Sciences - University of Warwick

subBacteriaEnvData <- function(countData, bacteriaType, sampleTypes) {

  subsetData <- countData[countData$Bacteria == bacteriaType, ]
  groupLength <- 3 * length(sampleTypes)

  firstCols <- c(which(colnames(subsetData) == "SamplingCode"),
                 which(colnames(subsetData) == "SamplingSite"),
                 which(colnames(subsetData) == "Replicate"))
  secondCols <- c(which(colnames(subsetData) == "River"),
                  which(colnames(subsetData) == "Location"))

  sumData <- subData(subsetData, groupLength, sampleTypes, "Env", "SampleType",
                     firstCols, secondCols)

  sumData$SampleType <- factor(sumData$SampleType, levels = unique(sumData$SampleType))
  sumData <- sumData[order(sumData$SamplingSite), ]
  colnames(sumData)[(length(firstCols) + 2):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "River", "Location")

  return(sumData)
}
