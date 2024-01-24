# ReadAndBarPlotHHChromoSelectCounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by household and order by household, sample type and sample id
#
# file ReadAndBarPlotHHChromoSelectCounts
#
# inputs
# 	plateCountHHCSFile  - file containing plate count data
#   plotsDir            - directory to store plots in
#   metaDataFile        - file containing river meta data
#   rawData             - True for absolute counts, false for cfu / ml
#   stoolReps           - Number of stool dilutions
#   waterReps           - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/12/22  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHChromoSelectCounts <- function(plateCountHHCSFile, plotsDir,
                                               metaDataFile, rawData, stoolReps,
                                               waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountHHCSFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  plateCSCountData <- numericCounts(plateCSCountData)

  waterData <- plateCSCountData[plateCSCountData$`Sample Type` == "Water", ]
  stoolData <- plateCSCountData[plateCSCountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateCSCountData <- rbind(waterData, stoolData)

  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- addHouseholds(plateCSCountData, metaData)

  plateCSCountData <- nameCols(plateCSCountData)

  plateCSCountData$SampleID <-
    as.numeric(substr(plateCSCountData$SampleID, 5, 8))

  plateCSCountData <- plateCSCountData %>% relocate(Klebsiella,
                                                    .after = coliforms)

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  stoolData <- plateCSCountData[plateCSCountData$SampleType == "Stool", ]

  convertCols <- c(which(colnames(plateCSCountData) == "E.coli"),
                   which(colnames(plateCSCountData) == "coliforms"),
                   which(colnames(plateCSCountData) == "Klebsiella"))

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps, rawData,
                                                    convertCols, bacteriaTypes,
                                                    "HH"))

  colnames(subPlateData)[which(colnames(subPlateData) == "Sample-ID")] <- "SampleID"

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11, length(subPlateData$Household))),
                                     subPlateData$SampleType,
                                     subPlateData$SampleID), ]
  sampleNumbers <- unique(subPlateData$SampleID)
  subPlateData$SampleID <- factor(subPlateData$SampleID,
                                     levels = sampleNumbers)

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
                                     c("Water", "Stool"))

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "HH", "ChromoSelect plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), "B", 4, plotsDir,
                     bacteriaTypes[iBacteria])
  }

}

