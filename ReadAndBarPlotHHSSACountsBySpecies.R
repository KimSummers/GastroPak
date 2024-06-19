# ReadAndBarPlotHHSSACountsBySpecies
#
# Read in plate count data and barchart plot it
# One bar chart for the water samples and one for faecal
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHSSACountsBySpecies
#
# inputs
# 	plateCountSSAFile - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   faecalReps        - Number of faecal dilutions
#   waterReps         - Number of water dilutions
#   stream            - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  26/02/24  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSSACountsBySpecies <- function(plateCountSSAFile, plotsDir,
                                                metaDataFile, rawData,
                                                faecalReps, waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountSSAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  plateSSACountData <- numericCounts(plateSSACountData)
  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- nameCols(plateSSACountData)

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  faecalData <- plateSSACountData[plateSSACountData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateSSACountData <- rbind(waterData, faecalData)

  plateSSACountData <- addHouseholds(plateSSACountData, metaData)

  plateSSACountData$SampleID <-
    as.numeric(substr(plateSSACountData$SampleID, 5, 8))

  plateSSACountData <- plateSSACountData %>% relocate(Salmonella,
                                                      .after = Shigella)
  plateSSACountData <- plateSSACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Shigella", "Salmonella", "E.coli")

  convertCols <- which(colnames(plateSSACountData) %in% bacteriaTypes)

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  faecalData <- plateSSACountData[plateSSACountData$SampleType == "Stool", ]

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
                     dataType = "HH", mainData = "Salmonella-Shigella agar plates",
                     subData = sampleType[iSampleType],
                     fillColours = c('grey', 'black', 'pink'),
                     weightOrVol = measureType, rowGraphs = 4, plotsDir = plotsDir,
                     plotTitle = paste("Household", stream, "SSA", sampleType[iSampleType]))
  }

}
