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

  subData <- subPlateData[subPlateData$Bacteria == bacteriaType, ]
  groupLength <- 3 * length(sampleTypes)

  firstCols <- c(which(colnames(countData) == "SampleID"),
                 which(colnames(countData) == "SamplingSite"))
  secondCols <- c(which(colnames(countData) == "Household"))

  for (iRow in (1:(nrow(subData) / groupLength)))
  {

    for (iSample in 1:length(sampleTypes))
    {
      meanVal <- mean(subData$MeanCfu[(subData$SampleID ==
                                        subData$SampleID[iRow * groupLength]) &
                                        (subData$SampleType == sampleTypes[iSample])],
                      na.rm = TRUE)
      sdVal <- sd(subData$MeanCfu[(subData$SampleID ==
                                     subData$SampleID[iRow * groupLength]) &
                                    (subData$SampleType == sampleTypes[iSample])],
                  na.rm = TRUE)
      sumDataRow <- cbind(subData[iRow * groupLength, firstCols],
                          sampleTypes[iSample], meanVal, sdVal,
                          subData[iRow * groupLength, secondCols])

      colnames(sumDataRow)[length(firstCols) + 1] <- "SampleType"

      if ((iRow == 1) & (iSample == 1))
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

    }

  }

  sumData$meanVal[sumData$meanVal == 0] <- NaN
  sumData$SampleType <- factor(sumData$SampleType, levels = unique(sumData$SampleType))
  colnames(sumData)[(length(firstCols) + 2):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "Household")

  return(sumData)
}
