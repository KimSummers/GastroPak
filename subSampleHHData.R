# subSampleHHData
#
# Get data for sample type, including mean and sd values
#
# file subSampleHHData
#
# inputs
# 	countData     - Count data to get sub data from
#   sampleType    - Include data with this sample type
#   bacteriaType  - Type of bacteria to get data for

# Version    Author       Date      Affiliation
# 1.00       J K Summers  10/01/24  Wellington Lab - School of Life Sciences - University of Warwick

subSampleHHData <- function(countData, sampleType, bacteriaTypes) {

  subData <- subPlateData[subPlateData$SampleType == sampleType, ]
  groupLength <- 3 * length(bacteriaTypes)

  firstCols <- c(which(colnames(subData) == "SampleID"),
                 which(colnames(subData) == "SamplingSite"),
                 which(colnames(subData) == "Replicate"))
  secondCols <- c(which(colnames(subData) == "Household"))

  for (iRow in (1:(nrow(subData) / groupLength)))
  {

    for (iBacteria in 1:length(bacteriaTypes))
    {
      meanVal <- mean(subData$MeanCfu[(subData$SampleID ==
                                         subData$SampleID[iRow * groupLength]) &
                                        (subData$Bacteria == bacteriaTypes[iBacteria])],
                      na.rm = TRUE)
      sdVal <- sd(subData$MeanCfu[(subData$SampleID ==
                                     subData$SampleID[iRow * groupLength]) &
                                    (subData$Bacteria == bacteriaTypes[iBacteria])],
                  na.rm = TRUE)
      sumDataRow <- cbind(subData[iRow * groupLength, firstCols],
                          bacteriaTypes[iBacteria], meanVal, sdVal,
                          subData[iRow * groupLength, secondCols])

      colnames(sumDataRow)[length(firstCols) + 1] <- "Bacteria"

      if ((iRow == 1) & (iBacteria == 1))
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

    }

  }

  sumData$meanVal[sumData$meanVal == 0] <- NaN
  sumData$Bacteria <- factor(sumData$Bacteria, levels = unique(sumData$Bacteria))
  sumData <- sumData[order(sumData$SamplingSite), ]
  colnames(sumData)[(length(firstCols) + 2):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "Household")

  return(sumData)
}
