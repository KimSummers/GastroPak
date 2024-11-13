# ReadAndBarPlotHHBGACountsBySpecies
#
# Read in plate count data and bar chart plot it
# One bar chart for the water samples and one for faecal
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHBGACountsBySpecies
#
# inputs
# 	plateCountBGAFile - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   faecalReps        - Number of faecal dilutions
#   waterReps         - Number of water dilutions
#   stream            - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  26/04/24  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHBGACountsBySpecies <- function(plateCountBGAFile, plotsDir,
                                               metaDataFile, rawData,
                                               faecalReps, waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountBGAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  plateBGACountData <- numericCounts(plateBGACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- nameCols(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  faecalData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateBGACountData <- rbind(waterData, faecalData)

  plateBGACountData <- addHouseholds(plateBGACountData, metaData)

  plateBGACountData$SampleID <-
    as.numeric(substr(plateBGACountData$SampleID, 5, 8))

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Salmonella", "E.coli")

  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  faecalData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(faecalData, faecalReps,
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

  locations <- unique(subPlateData$SamplingSite)

  for (ilocation in 1:length(locations))
  {
    subPlateLocData <- subPlateData[subPlateData$SamplingSite == locations[ilocation], ]

    for (iSampleType in 1:length(sampleType))
    {
      bactSubData <- subSampleHHData(subPlateLocData, sampleType[iSampleType],
                                     bacteriaTypes)

      if (sampleType[iSampleType] == "Water")
      {
        measureType = "V"
      }else
      {
        measureType = "W"
      }

      BarPlotGastroPak(graphData = bactSubData, xAxisLabels = NULL, plotType = "Species",
                       dataType = "HH", mainData = "Brilliant green agar plates",
                       subData = sampleType[iSampleType],
                       fillColours = c('palevioletred', 'yellowgreen'),
                       weightOrVol = measureType, rowGraphs = 4, plotsDir = plotsDir,
                       plotTitle = paste("Household", stream, locations[ilocation],
                                         "BGA", sampleType[iSampleType]))
    }

  }

}
