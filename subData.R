# subData
#
# Get data for bacteria type, including mean and sd values
#
# file subData
#
# inputs
# 	subsetData    - Data to get summary data from
#   groupLength   - How many data points in each group
#   groupTypes    - What groups of data are required
#   campaignType  - Household or environmental
#   dataType      - What parameter data is grouped on
#   firstCols     - First set of columns needed
#   secondCols    - Second set of columns needed

# Version    Author       Date      Affiliation
# 1.00       J K Summers  10/01/24  Wellington Lab - School of Life Sciences - University of Warwick
subData <- function(subsetData, groupLength, groupTypes, campaignType, dataType,
                    firstCols, secondCols) {

  sumData <- NULL

  for (iRow in (1:(nrow(subsetData) / groupLength)))
  {

    if (length(groupTypes) > 0)
    {

      for (iGroup in 1:length(groupTypes))
      {

        if (campaignType == "Env")
        {
          dataPoints <- subsetData[(subsetData$SamplingSite ==
                                      subsetData$SamplingSite[iRow * groupLength]), ]
        }else
        {
          dataPoints <- subsetData[(subsetData$SampleID ==
                                      subsetData$SampleID[iRow * groupLength]), ]
        }

        if (dataType == "SampleType")
        {

          if (campaignType == "Env")
          {
            dataPoints <- dataPoints[dataPoints$SampleType == groupTypes[iGroup], ]
          }

        }else
        {

          if (dataType == "Bacteria")
          {
            dataPoints <- dataPoints[dataPoints$Bacteria == groupTypes[iGroup], ]
          }else
          {

            if (dataType == "Media")
            {
              dataPoints <- dataPoints[dataPoints$Media == dataPoints$Media[iRow * groupLength]]
            }

          }

        }

        meanVal <- mean(dataPoints$MeanCfu, na.rm = TRUE)
        sdVal <- sd(dataPoints$MeanCfu, na.rm = TRUE)

        if (dataType == "Media")
        {
          sumDataRow <- cbind(dataPoints[1, firstCols], meanVal, sdVal,
                              dataPoints[1, secondCols])
        }else
        {
          sumDataRow <- cbind(dataPoints[1, firstCols], groupTypes[iGroup],
                              meanVal, sdVal, dataPoints[1, secondCols])
          colnames(sumDataRow)[length(firstCols) + 1] <- dataType
        }

        sumData <- rbind(sumData, sumDataRow)
      }

    }else
    {

      if (campaignType == "Env")
      {
        dataPoints <- subsetData[(subsetData$SamplingSite ==
                                    subsetData$SamplingSite[iRow * groupLength]), ]
      }else
      {
        dataPoints <- subsetData[(subsetData$SampleID ==
                                    subsetData$SampleID[iRow * groupLength]), ]
      }

      dataPoint <- dataPoints[dataPoints$Media == dataPoints$Media[iRow * groupLength]]

      meanVal <- mean(dataPoints$MeanCfu, na.rm = TRUE)
      sdVal <- sd(dataPoints$MeanCfu, na.rm = TRUE)

      sumDataRow <- cbind(dataPoints[1, firstCols], meanVal, sdVal,
                          dataPoints[1, secondCols])

      sumData <- rbind(sumData, sumDataRow)
    }

  }

  sumData$meanVal[sumData$meanVal == 0] <- NaN
  return(sumData)
}
