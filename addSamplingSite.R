# addSamplingSite
#
# Add sampling site and season to the data with sample id
#
# file addSamplingSite
#
# inputs
# 	baseData      - Data to have sampling site and season added
#   metaData      - Sampling site and season data to add

# Version    Author       Date      Affiliation
# 1.00       J K Summers  22/04/24  Wellington Lab - School of Life Sciences - University of Warwick

addSamplingSite <- function(baseData, metaData) {

  samplingData <- NULL

  for (iRow in 1:nrow(baseData))
  {
    samplingSite <- metaData$Site_ID[metaData$Barcode == baseData$SampleID[iRow]]
    seasonHC <- metaData$Season[metaData$Barcode == baseData$SampleID[iRow]]
    seasonWD <- metaData$Season_WD[metaData$Barcode == baseData$SampleID[iRow]]
    season <- paste(seasonHC, seasonWD)

    samplingData <- rbind(samplingData, cbind(samplingSite, season))
  }

  baseData <- cbind(baseData, samplingData)
  colnames(baseData)[ncol(baseData) - 1] <- "SamplingSite"
  colnames(baseData)[ncol(baseData)] <- "Season"

  return(baseData)
}
