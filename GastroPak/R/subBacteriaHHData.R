# subBacteriaHHData
#
# Get data for bacteria type, including mean and sd values for household data
#
# file subBacteriaHHData
#
# inputs
# 	countData     - Count data to get sub data from
#   bacteriaType  - Include data with this bacteria type
#   sampleTypes   - Type of sample to get data for

# Version    Author       Date      Affiliation
# 1.00       J K Summers  10/01/24  Wellington Lab - School of Life Sciences - University of Warwick

subBacteriaHHData <- function(countData, bacteriaType, sampleTypes) {

  sumData <- NULL

  for (iSampleType in 1:length(sampleTypes))
  {
    subsetData <- countData[countData$Bacteria == bacteriaType &
                              countData$SampleType == sampleTypes[iSampleType], ]
    groupLength <- 3

    firstCols <- c(which(colnames(subsetData) == "SampleID"),
                   which(colnames(subsetData) == "SamplingSite"))
    secondCols <- c(which(colnames(subsetData) == "Household"))

    sumData <- rbind(sumData, subData(subsetData = subsetData, groupLength = groupLength,
                                      groupTypes = sampleTypes[iSampleType], campaignType = "HH",
                                      dataType = "SampleType", firstCols = firstCols,
                                      secondCols = secondCols))
  }

  sumData$SampleType <- factor(sumData$SampleType, levels = unique(sumData$SampleType))
  colnames(sumData)[(length(firstCols) + 2):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "Household")

  return(sumData)
}
