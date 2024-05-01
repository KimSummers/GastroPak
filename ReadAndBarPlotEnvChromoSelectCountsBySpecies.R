# ReadAndBarPlotEnvChromoSelectCountsBySpecies
#
# Read in plate count data and bar chart plot it
# One bar chart for the E. coli counts, one for coliforms and one for Klebsiella
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvChromoSelectCountsBySpecies
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
# 1.00       J K Summers  15/07/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvChromoSelectCountsBySpecies <- function(plateCountEnvCSFile,
                                                         plotsDir, metaDataFile,
                                                         rawData, sedimentReps,
                                                         waterReps, season) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountEnvCSFile)
  metaData <- read_csv(metaDataFile)

  plateCSCountData <- numericCounts(plateCSCountData)
  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- nameCols(plateCSCountData)

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$SampleType == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateCSCountData <- rbind(waterData, sedimentData)

  plateCSCountData <- addRivers(plateCSCountData, metaData)

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

  subPlateData$Bacteria <- factor(subPlateData$Bacteria,
                                  levels = bacteriaTypes)

  subPlateData <- subPlateData[order(subPlateData$SamplingCode, subPlateData$Bacteria), ]

  for (iSampleType in 1:length(sampleType))
  {
    sumData <- subSampleEnvData(subPlateData, sampleType[iSampleType], bacteriaTypes)

    if (sampleType[iSampleType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(graphData = sumData,
                     xAxisLabels = c("Upstream", "Midstream", "Downstream"),
                     plotType = "Species", dataType = "Env",
                     mainData = "ChromoSelect plates", subData = sampleType[iSampleType],
                     fillColours = c('darkblue','tomato2', 'pink'),
                     weightOrVol = measureType, rowGraphs = 3, plotsDir = plotsDir,
                     plotTitle = paste("Env", season, "CS", sampleType[iSampleType]))
  }

}



