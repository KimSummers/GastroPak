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

  subData <- subPlateData[subPlateData$Bacteria == bacteriaType, ]
  groupLength <- 3 * length(sampleTypes)

  firstCols <- c(which(colnames(subData) == "SamplingCode"),
                 which(colnames(subData) == "SamplingSite"),
                 which(colnames(subData) == "Replicate"))
  secondCols <- c(which(colnames(subData) == "River"),
                  which(colnames(subData) == "Location"))

  for (iRow in (1:(nrow(subData) / groupLength)))
  {

    for (iSample in 1:length(sampleTypes))
    {
      meanVal <- mean(subData$MeanCfu[(subData$SamplingSite ==
                                         subData$SamplingSite[iRow * groupLength]) &
                                        (subData$SampleType == sampleTypes[iSample])],
                      na.rm = TRUE)
      sdVal <- sd(subData$MeanCfu[(subData$SamplingSite ==
                                     subData$SamplingSite[iRow * groupLength]) &
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
  sumData <- sumData[order(sumData$SamplingSite), ]
  colnames(sumData)[(length(firstCols) + 2):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "River", "Location")

  return(sumData)
}
