# ReadAndBarPlotEnvChromoSelectCounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvChromoSelectCounts
#
# inputs
# 	plateCountCSFile  - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/12/22  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvChromoSelectCounts <- function(plateCountEnvCSFile, plotsDir,
                                                metaDataFile, rawData,
                                                sedimentReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountEnvCSFile)
  metaData <- read.csv(metaDataFile)

  plateCSCountData <- numericCounts(plateCSCountData)

  waterData <- plateCSCountData[plateCSCountData$`Sample Type` == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateCSCountData <- rbind(waterData, sedimentData)

  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- addRivers(plateCSCountData, metaData)

  plateCSCountData <- nameCols(plateCSCountData)

  plateCSCountData <- plateCSCountData %>% relocate(Klebsiella,
                                                    .after = coliforms)

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")
  convertCols <- which(colnames(plateCSCountData) %in% bacteriaTypes)

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$SampleType == "Sediment", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "Env")
  subPlateData <- rbind(subPlateData, averageCounts(sedimentData, sedimentReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes, "Env"))

  sampleCodes <- unique(subPlateData$SamplingCode)
  subPlateData$SamplingCode <- factor(subPlateData$SamplingCode,
                                      levels = sampleCodes)

  sampleSites <- unique(subPlateData$SamplingSite)
  subPlateData$SamplingSite <- factor(subPlateData$SamplingSite,
                                      levels = sampleSites)
  sampleType <- unique(subPlateData$SampleType)
  subPlateData$SampleType <- factor(subPlateData$SampleType,
                                    levels = sampleType)

  subPlateData <- subPlateData[order(subPlateData$SamplingSite), ]

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subBacteriaEnvData(subPlateData, bacteriaTypes[iBacteria],
                                      c("Sediment", "Water"))

    BarPlotGastroPak(bactSubData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "Chromocult plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), plotsDir, bacteriaTypes[iBacteria])
  }

}
