# subCompHHData
#
# Get data for sample type, including mean and sd values
#
# file subCompHHData
#
# inputs
# 	countData     - Count data to get sub data from
#   sampleType    - Include data with this sample type

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/02/24  Wellington Lab - School of Life Sciences - University of Warwick

subCompHHData <- function(countData, sampleType) {

  subsetData <- subPlateData[subPlateData$SampleType == sampleType, ]
  groupLength <- 3

  firstCols <- c(which(colnames(subsetData) == "SampleID"),
                 which(colnames(subsetData) == "SamplingSite"),
                 which(colnames(subsetData) == "Replicate"),
                 which(colnames(subsetData) == "Media"))
  secondCols <- c(which(colnames(subsetData) == "Household"))

  sumData <- subData(subsetData, groupLength, NULL, "HH", "Media", firstCols,
                     secondCols)

  sumData$Media <- factor(sumData$Media, levels = unique(sumData$Media))
  sumData <- sumData[order(sumData$SamplingSite), ]
  colnames(sumData)[(length(firstCols) + 1):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "Household")

  return(sumData)
}
