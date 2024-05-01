# subSampleEnvData
#
# Get data for sample type, including mean and sd values
#
# file subSampleEnvData
#
# inputs
# 	countData     - Count data to get sub data from
#   sampleType    - Include data with this sample type
#   bacteriaType  - Type of bacteria to get data for

# Version    Author       Date      Affiliation
# 1.00       J K Summers  10/01/24  Wellington Lab - School of Life Sciences - University of Warwick

subSampleEnvData <- function(countData, sampleType, bacteriaTypes) {

  subsetData <- countData[countData$SampleType == sampleType, ]
  groupLength <- 3 * length(bacteriaTypes)

  firstCols <- c(which(colnames(subsetData) == "SamplingCode"),
                 which(colnames(subsetData) == "SamplingSite"),
                 which(colnames(subsetData) == "Replicate"))
  secondCols <- c(which(colnames(subsetData) == "River"),
                  which(colnames(subsetData) == "Location"))

  sumData <- subData(subsetData, groupLength, bacteriaTypes, "Env", "Bacteria",
                     firstCols, secondCols)

  sumData$Bacteria <- factor(sumData$Bacteria, levels = unique(sumData$Bacteria))
  sumData <- sumData[order(sumData$SamplingSite), ]
  colnames(sumData)[(length(firstCols) + 2):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "River", "Location")

  return(sumData)
}
