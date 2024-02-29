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

  subData <- subPlateData[subPlateData$SampleType == sampleType, ]
  groupLength <- 3

  firstCols <- c(which(colnames(subData) == "SampleID"),
                 which(colnames(subData) == "SamplingSite"),
                 which(colnames(subData) == "Replicate"),
                 which(colnames(subData) == "Media"))
  secondCols <- c(which(colnames(subData) == "Household"))

  for (iRow in (1:(nrow(subData) / groupLength)))
  {
    meanVal <- mean(subData$MeanCfu[(subData$SampleID ==
                                       subData$SampleID[iRow * groupLength]) &
                                      (subData$Media == subData$Media[iRow * groupLength])],
                    na.rm = TRUE)
    sdVal <- sd(subData$MeanCfu[(subData$SampleID ==
                                   subData$SampleID[iRow * groupLength]) &
                                  (subData$Media == subData$Media[iRow * groupLength])],
                na.rm = TRUE)
    sumDataRow <- cbind(subData[iRow * groupLength, firstCols], meanVal, sdVal,
                        subData[iRow * groupLength, secondCols])

    if (iRow == 1)
    {
      sumData <- sumDataRow
    }else
    {
      sumData <- rbind(sumData, sumDataRow)
    }

  }

  sumData$meanVal[sumData$meanVal == 0] <- NaN
  sumData$Media <- factor(sumData$Media, levels = unique(sumData$Media))
  sumData <- sumData[order(sumData$SamplingSite), ]
  colnames(sumData)[(length(firstCols) + 1):ncol(sumData)] <-
    c("MeanCfu", "StdDev", "Household")

  return(sumData)
}
