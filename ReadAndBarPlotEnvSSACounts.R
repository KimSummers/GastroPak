# ReadAndBarPlotEnvSSACounts
#
# Read in plate count data and barchart plot it
# One bar chart for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvSSACounts
#
# inputs
# 	plateCountFile    - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions
#   season            - The season the data is from

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/01/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvSSACounts <- function(plateCountSSAFile, plotsDir, metaDataFile,
                                       rawData, sedimentReps, waterReps, season) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountSSAFile)
  metaData <- read_csv(metaDataFile)

  plateSSACountData <- numericCounts(plateSSACountData)
  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- nameCols(plateSSACountData)

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  sedimentData <- plateSSACountData[plateSSACountData$SampleType == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateSSACountData <- rbind(waterData, sedimentData)

  plateSSACountData <- addRivers(plateSSACountData, metaData)

  bacteriaTypes <- c("Shigella", "Salmonella", "E.coli")

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  sedimentData <- plateSSACountData[plateSSACountData$SampleType == "Sediment", ]

  convertCols <- which(colnames(plateSSACountData) %in% bacteriaTypes)

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
    bactSubData <- subBacteriaEnvData(countData = subPlateData,
                                      bacteriaType = bacteriaTypes[iBacteria],
                                      sampleTypes = c("Sediment", "Water"))

    BarPlotGastroPak(graphData = bactSubData,
                     xAxisLabels = c("Upstream", "Midstream", "Downstream"),
                     plotType = "Sample", dataType = "Env",
                     mainData = "Salmonella-Shigella agar plates", bacteriaTypes[iBacteria],
                     fillColours = c('chocolate4', 'skyblue'), weightOrVol = "B",
                     rowGraphs = 3, plotsDir = plotsDir,
                     plotTitle = paste("Env", season, "SSA", bacteriaTypes[iBacteria]))
  }

}
