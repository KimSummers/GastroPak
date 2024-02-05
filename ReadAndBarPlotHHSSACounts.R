# ReadAndBarPlotHHSSACounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Salmonella, one for Shigella
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHSSACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   stoolReps         - Number of stool dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/05/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSSACounts <- function(plateCountHHSSAFile, plotsDir, metaDataFile,
                                      rawData, stoolReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  plateSSACountData <- numericCounts(plateSSACountData)

  waterData <- plateSSACountData[plateSSACountData$`Sample Type` == "Water", ]
  stoolData <- plateSSACountData[plateSSACountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateSSACountData <- rbind(waterData, stoolData)
  plateSSACountData <- nameCols(plateSSACountData)

  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- addHouseholds(plateSSACountData, metaData)

  plateSSACountData <- nameCols(plateSSACountData)

  plateSSACountData$SampleID <-
    as.numeric(substr(plateSSACountData$SampleID, 5, 8))

  plateSSACountData <- plateSSACountData %>% relocate(Shigella,
                                                    .after = Salmonella)
  plateSSACountData <- plateSSACountData %>% relocate(E.coli,
                                                      .after = Shigella)

  bacteriaTypes <- c("Salmonella", "Shigella", "E.coli")
  convertCols <- which(colnames(plateSSACountData) %in% bacteriaTypes)

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  stoolData <- plateSSACountData[plateSSACountData$SampleType == "Stool", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes, "HH"))

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household,
                                                       11, length(subPlateData$Household))),
                                     subPlateData$SampleType,
                                     subPlateData$SampleID), ]
  sampleIDs <- unique(subPlateData$SampleID)
  subPlateData$SampleID <- factor(subPlateData$SampleID,
                                     levels = sampleIDs)

  sampleSites <- unique(subPlateData$SamplingSite)
  subPlateData$SamplingSite <- factor(subPlateData$SamplingSite,
                                      levels = sampleSites)
  sampleType <- unique(subPlateData$SampleType)
  subPlateData$SampleType <- factor(subPlateData$SampleType,
                                    levels = sampleType)

  households <- unique(subPlateData$Household)
  subPlateData$Household <- factor(subPlateData$Household,
                                   levels = households)

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subBacteriaHHData(subPlateData, bacteriaTypes[iBacteria],
                                     c("Stool", "Water"))

    BarPlotGastroPak(bactSubData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "HH", "Salmonella-Shigella agar plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), "B", 4, plotsDir,
                     paste("Household", bacteriaTypes[iBacteria]))
  }

}
