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
#   sampleType    - Type of sample
#   campaignType  - "Env" for environmental, "HH" for household
#.  incMedia      - Should a media column be included (true for comparison datasets)

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

averageCounts <- function(countData, numDilutions, rawData, convertCols,
                          bacteriaTypes, campaignType, incMedia = FALSE) {

  if (!rawData)
  {
    countData <- convertCfuCounts(countData, numDilutions, convertCols)
  }

  if (campaignType == "Env")
  {
    firstCols <- c(which(colnames(countData) == "SamplingCode"))
  }else
  {
    firstCols <- c(which(colnames(countData) == "SampleID"))
  }

  firstCols <- c(firstCols,
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

        if (!is.na(countData[((iRow - 1) * numDilutions + iRep), iBacteria]))
        {
          totalVol <- totalVol + 10^(-iRep)
        }

      }

      cfuPerml <- totalCfu / totalVol

      newRow <- cbind(newRow, bacteriaTypes[iName], totalCfu,
                      totalVol, cfuPerml)

      if (campaignType == "Env")
      {
        newRow <- cbind(newRow, countData$River[((iRow - 1) * numDilutions + 1)],
                        countData$RiverLocation[((iRow - 1) * numDilutions + 1)])
        colnames(newRow)[(length(firstCols) + 1):ncol(newRow)] <-
          c("Bacteria", "Total colony count", "Total volume count", "MeanCfu",
             "River", "Location")
      }else
      {
        newRow <- cbind(newRow, countData$Household[((iRow - 1) * numDilutions + 1)])
        colnames(newRow)[(length(firstCols) + 1):ncol(newRow)] <-
          c("Bacteria", "Total colony count", "Total volume count", "MeanCfu",
             "Household")
      }

      if (incMedia)
      {
        newRow <- cbind(newRow, countData$Media[((iRow - 1) * numDilutions + 1)])
        colnames(newRow)[ncol(newRow)] <- "Media"
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
