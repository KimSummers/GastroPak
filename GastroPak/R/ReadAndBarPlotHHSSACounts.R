# ReadAndBarPlotHHSSACounts
#
# Read in plate count data and bar chart plot it
# One bar chart for the water counts and one for faecal
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHSSACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile    - file containing river meta data
#   rawData         - True for absolute counts, false for cfu / ml
#   faecalReps      - Number of faecal dilutions
#   waterReps       - Number of water dilutions
#   stream          - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/05/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSSACounts <- function(plateCountHHSSAFile, plotsDir, metaDataFile,
                                      rawData, faecalReps, waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  # Ensure count data is numeric, collate by species and name columns correctly
  plateSSACountData <- numericCounts(plateSSACountData)
  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- nameCols(plateSSACountData)

  # Fix gaps in data
  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  faecalData <- plateSSACountData[plateSSACountData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateSSACountData <- rbind(waterData, faecalData)

  # Add meta data on households to table
  plateSSACountData <- addHouseholds(plateSSACountData, metaData)

  plateSSACountData$SampleID <-
    as.numeric(substr(plateSSACountData$SampleID, 5, 8))

  plateSSACountData <- plateSSACountData %>% relocate(Shigella,
                                                    .after = Salmonella)
  plateSSACountData <- plateSSACountData %>% relocate(E.coli,
                                                      .after = Shigella)

  bacteriaTypes <- c("Salmonella", "Shigella", "E.coli")
  convertCols <- which(colnames(plateSSACountData) %in% bacteriaTypes)

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  faecalData <- plateSSACountData[plateSSACountData$SampleType == "Stool", ]

  # Average plate counts per species and per sample accross all dilutions
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

  subPlateData$SampleType[subPlateData$SampleType == "Stool"] <- "Faecal"
  sampleType <- unique(subPlateData$SampleType)
  subPlateData$SampleType <- factor(subPlateData$SampleType,
                                    levels = sampleType)

  households <- unique(subPlateData$Household)
  subPlateData$Household <- factor(subPlateData$Household,
                                   levels = households)

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subBacteriaHHData(countData = subPlateData,
                                     bacteriaType = bacteriaTypes[iBacteria],
                                     sampleTypes = c("Faecal", "Water"))

    BarPlotGastroPak(graphData = bactSubData, xAxisLabels = NULL, plotType = "Sample",
                     dataType = "HH", mainData = "Salmonella-Shigella agar plates",
                     subData = bacteriaTypes[iBacteria],
                     fillColours = c('chocolate4', 'skyblue'), weightOrVol = "B",
                     rowGraphs = 4, plotsDir = plotsDir,
                     plotTitle = paste("Household", stream, "SSA", bacteriaTypes[iBacteria]))
  }

}
