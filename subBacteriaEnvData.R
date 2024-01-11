# subBacteriaData
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

subBacteriaData <- function(countData, bacteriaType, sampleTypes) {

  subData <- subPlateData[subPlateData$Bacteria == bacteriaType, ]
  groupLength <- 3 * length(bacteriaTypes)

  firstCols <- c(which(colnames(subData) == "SamplingCode"),
                 which(colnames(subData) == "SamplingSite"),
                 which(colnames(subData) == "Replicate"))
  secondCols <- c(which(colnames(subData) == "River"),
                  which(colnames(subData) == "Location"))

  for (iRow in (1:(nrow(subData) / groupLength)))
  {

    for (iBacteria in 1:length(bacteriaTypes))
    {
      meanVal <- mean(subData$MeanCfu[(subData$SamplingSite ==
                                         subData$SamplingSite[iRow * groupLength]) &
                                        (subData$Bacteria == bacteriaTypes[iBacteria])],
                      na.rm = TRUE)
      sdVal <- sd(subData$MeanCfu[(subData$SamplingSite ==
                                     subData$SamplingSite[iRow * groupLength]) &
                                    (subData$Bacteria == bacteriaTypes[iBacteria])],
                  na.rm = TRUE)
      sumDataRow <- cbind(subData[iRow * groupLength, firstCols],
                          bacteriaTypes[iBacteria], meanVal, sdVal,
                          subData[iRow * groupLength, secondCols])

      colnames(sumDataRow)[4] <- "Bacteria"

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
  colnames(sumData)[5:8] <- c("MeanCfu", "StdDev", "River", "Location")

  return(sumData)
}
