# ReadAndBarPlotEnvBGACountsBySpecies
#
# Read in plate count data and bar chart plot it
# One bar chart for the Salmonella counts and one for E. coli
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvBGACountsBySpecies
#
# inputs
# 	plateCountBGAFile - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions
#   season            - which season campaign the data was collected in

# Version    Author       Date      Affiliation
# 1.00       J K Summers  17/07/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvBGACountsBySpecies <- function(plateCountBGAFile, plotsDir,
                                                metaDataFile, rawData,
                                                sedimentReps, waterReps,
                                                season) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountBGAFile)
  plateBGACountData <- plateBGACountData[!is.na(plateBGACountData$`Sample Type`), ]
  metaData <- read_csv(metaDataFile)

  # Ensure count data is numeric, collate by species and name columns correctly
  plateBGACountData <- numericCounts(plateBGACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- nameCols(plateBGACountData)

  # Fix gaps in data
  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateBGACountData <- rbind(waterData, sedimentData)

  # Add meta data on rivers to table
  plateBGACountData <- addRivers(plateBGACountData, metaData)

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)
  bacteriaTypes <- c("Salmonella", "E.coli")

  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

  # Average plate counts per species and per sample across all dilutions
  subPlateData <- averageCounts(countData = waterData, numDilutions = waterReps, 
                                rawData = rawData, convertCols = convertCols, 
                                bacteriaTypes =  bacteriaTypes, campaignType = "Env")
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

  subPlateData <- subPlateData[order(subPlateData$SamplingCode, subPlateData$Bacteria), ]

  for (iSampleType in 1:length(sampleType))
  {
    sumData <- subSampleEnvData(countData = subPlateData, 
                                sampleType = sampleType[iSampleType], 
                                bacteriaTypes = bacteriaTypes)

    if (sampleType[iSampleType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(graphData = sumData, xAxisLabels = c("Upstream", "Midstream", "Downstream"), 
                     plotType =  "Species", dataType = "Env", mainData = "Brilliant green agar", 
                     subData = sampleType[iSampleType], 
                     fillColours = c('palevioletred', 'yellowgreen'), 
                     weightOrVol = measureType, rowGraphs = 3, plotsDir = plotsDir,
                     plotTitle = paste("Env", season, "BGA", sampleType[iSampleType]))
  }

}

