# averageCounts
#
# Average counts accross a dilution series
#
# file averageCounts
#
# inputs
# 	countData     - Initial Data to have counts averaged
#   numDilutions  - Number of replicates of data
#   rawData       - True if data is number of colonies, false for cfu / ml
#   convertCols   - The columns with counts to be converted
#   bacteriaTypes - Names of bacteria for which we have counts
#   campaignType  - "Env" for environmental, "HH" for household

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

averageCounts <- function(countData, numDilutions, rawData, convertCols,
                          bacteriaTypes, campaignType) {

  if (!rawData)
  {
    countData <- convertCfuCounts(countData, numDilutions, convertCols)
  }

  firstCols <- c(which(colnames(countData) == "SamplingCode"),
                 which(colnames(countData) == "SamplingSite"),
                 which(colnames(countData) == "SampleType"),
                 which(colnames(countData) == "Replicate"))

  for (iRow in (1:(nrow(countData) / numDilutions)))
  {
    iName <- 0

    for (iBacteria in convertCols)
    {
      iName <- iName + 1
      newRow <- countData[((iRow - 1) * numDilutions + 1), firstCols]

      totalCfu <- sum(countData[((iRow - 1) * numDilutions + 1):(iRow * numDilutions),
                                iBacteria], na.rm = TRUE)
      totalVol <- 0

      for (iRep in 1:numDilutions)
      {

        if (!is.na(countData[(iRow - 1) * numDilutions + 1, iBacteria]))
        {
          totalVol <- totalVol + 10^(-iRep)
        }

      }

      cfuPerml <- totalCfu / totalVol

      newRow <- cbind(newRow, bacteriaTypes[iName], totalCfu, totalVol,
                      cfuPerml)

      if (campaignType == "Env")
      {
          newRow <- cbind(newRow, countData$River[((iRow - 1) * numDilutions + 1)],
                          countData$Location[((iRow - 1) * numDilutions + 1)])
          colnames(newRow)[5:10] <- c("Bacteria", "Total colony count",
                                      "Total volume count", "MeanCfu", "River",
                                      "Location")
      }else
      {
        newRow <- cbind(newRow, countData$Household[((iRow - 1) * numDilutions + 1)])
        colnames(newRow)[5:9] <- c("Bacteria", "Total colony count",
                                    "Total volume count", "MeanCfu", "Household")
      }

      if (iRow == 1 & iName == 1)
      {
        subPlateData <- newRow
      }else
      {
        subPlateData <- rbind(subPlateData, newRow)
      }

    }

  }

  subPlateData$MeanCfu <- as.numeric(subPlateData$MeanCfu)

  return(subPlateData)
}