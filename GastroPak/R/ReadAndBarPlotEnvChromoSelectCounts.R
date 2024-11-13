# ReadAndBarPlotEnvChromoSelectCounts
#
# Read in plate count data and bar chart plot it
# One bar chart for the E. coli counts, one for coliforms and one for Klebsiella
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
#   season            - The season the data is from

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/12/22  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvChromoSelectCounts <- function(plateCountCSFile, plotsDir,
                                                metaDataFile, rawData,
                                                sedimentReps, waterReps,
                                                season) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountCSFile)
  metaData <- read_csv(metaDataFile)

  # Ensure count data is numeric, collate by species and name columns correctly
  plateCSCountData <- numericCounts(plateCSCountData)
  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- nameCols(plateCSCountData)

  # Fix gaps in data
  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$SampleType == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateCSCountData <- rbind(waterData, sedimentData)

  # Add meta data on rivers to table
  plateCSCountData <- addRivers(plateCSCountData, metaData)

  plateCSCountData <- plateCSCountData %>% relocate(Klebsiella,
                                                    .after = coliforms)

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")
  convertCols <- which(colnames(plateCSCountData) %in% bacteriaTypes)

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$SampleType == "Sediment", ]

  # Average plate counts per species and per sample across all dilutions
  subPlateData <- averageCounts(countData = waterData, numDilutions = waterReps, 
                                rawData = rawData, convertCols = convertCols,
                                bacteriaTypes = bacteriaTypes, campaignType = "Env")
  subPlateData <- rbind(subPlateData, averageCounts(countData = sedimentData, 
                                                    numDilutions = sedimentReps,
                                                    rawData = rawData, 
                                                    convertCols = convertCols,
                                                    bacteriaTypes = bacteriaTypes, 
                                                    campaignType = "Env"))

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
                     plotType = "Sample", dataType = "Env", mainData = "ChromoSelect plates",
                     subData = bacteriaTypes[iBacteria],
                     fillColours = c('chocolate4', 'skyblue'), weightOrVol = "B",
                     rowGraphs = 3, plotsDir = plotsDir,
                     plotTitle = paste("Env", season, "CS", bacteriaTypes[iBacteria]))
  }

}
