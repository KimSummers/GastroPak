# addRivers
#
# Add river data to the plate counts
#
# file addRivers
#
# inputs
# 	countData     - Count data to have river data added
#   metaData      - River data to add

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

addRivers <- function(countData, metaData) {

  for (iRow in 1:nrow(countData))
  {
    river <- metaData$River[metaData$Sampling.Site == countData$`Sampling Site`[iRow]]
    sampleLoc <- metaData$Location[metaData$Sampling.Site == countData$`Sampling Site`[iRow]]

    if (iRow == 1)
    {
      rivers <- cbind(river, sampleLoc)
    }else
    {
      rivers <- rbind(rivers, cbind(river, sampleLoc))
    }

  }

  countData <- cbind(countData, rivers)
  colnames(countData)[ncol(countData) - 1] <- "River"
  colnames(countData)[ncol(countData)] <- "Location"

  return(countData)
}
