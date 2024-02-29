# ReadAndBarPlotEnvBGACounts
#
# Read in plate count data and barchart plot it
# One bar char for the Salmonella counts and one for E. coli
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvBGACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions
#   season            - The season the data is from

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/01/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvBGACounts <- function(plateCountEnvBGAFile, plotsDir, metaDataFile,
                                       rawData, sedimentReps, waterReps,
                                       season) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountEnvBGAFile)
  metaData <- read.csv(metaDataFile)

  plateBGACountData <- numericCounts(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$`Sample Type` == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateBGACountData <- rbind(waterData, sedimentData)

  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- addRivers(plateBGACountData, metaData)

  plateBGACountData <- nameCols(plateBGACountData)

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)
  bacteriaTypes <- c("Salmonella", "E.coli")
  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

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
                     "Env", "Brilliant Green Agar Plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), "B", 3, plotsDir,
                     paste(season, "BBA", bacteriaTypes[iBacteria]))
  }

}
