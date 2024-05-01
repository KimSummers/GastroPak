# ReadAndBarPlotHHChromoSelectCountsBySpecies
#
# Read in plate count data and barchart plot it
# One bar chart for the water samples and one for faecal
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHChromoSelectCountsBySpecies
#
# inputs
# 	plateCountCSFile  - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   faecalReps        - Number of faecal dilutions
#   waterReps         - Number of water dilutions
#   stream            - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  28/04/24  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHChromoSelectCountsBySpecies <- function(plateCountCSFile,
                                                        plotsDir, metaDataFile,
                                                        rawData, faecalReps,
                                                        waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountCSFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  plateCSCountData <- numericCounts(plateCSCountData)
  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- nameCols(plateCSCountData)

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  faecalData <- plateCSCountData[plateCSCountData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateCSCountData <- rbind(waterData, faecalData)

  plateCSCountData <- addHouseholds(plateCSCountData, metaData)

  plateCSCountData$SampleID <-
    as.numeric(substr(plateCSCountData$SampleID, 5, 8))

  plateCSCountData <- plateCSCountData %>% relocate(Klebsiella, .after = coliforms)

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")

  convertCols <- which(colnames(plateCSCountData) %in% bacteriaTypes)

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  faecalData <- plateCSCountData[plateCSCountData$SampleType == "Stool", ]

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

  for (iSampleType in 1:length(sampleType))
  {
    bactSubData <- subSampleHHData(subPlateData, sampleType[iSampleType], bacteriaTypes)

    if (sampleType[iSampleType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(graphData = bactSubData, xAxisLabels = NULL, plotType = "Species",
                     dataType = "HH", mainData = "ChromoSelect plates",
                     subData = sampleType[iSampleType],
                     fillColours = c('darkblue','tomato2', 'pink'),
                     weightOrVol = measureType, rowGraphs = 4, plotsDir = plotsDir,
                     plotTitle = paste("Household", stream, "CS", sampleType[iSampleType]))
  }

}
