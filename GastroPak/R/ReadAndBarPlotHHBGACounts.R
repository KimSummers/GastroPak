# ReadAndBarPlotHHBGACounts
#
# Read in plate count data and bar chart plot it
# One bar chart for the E. coli counts and one for Salmonella
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHBGACounts
#
# inputs
# 	plateCountHHBGAFile     - file containing plate count data
#   plotsDir                - directory to store plots in
#   metaDataFile            - file containing river meta data
#   rawData                 - True for absolute counts, false for cfu / ml
#   faecalReps              - Number of faecal dilutions
#   waterReps               - Number of water dilutions
#   stream                  - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  16/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHBGACounts <- function(plateCountHHBGAFile, plotsDir, metaDataFile,
                                      rawData, faecalReps, waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountHHBGAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  # Ensure count data is numeric, collate by species and name columns correctly
  plateBGACountData <- numericCounts(countData = plateBGACountData)
  plateBGACountData <- coloursToSpeciesBGA(countData = plateBGACountData)
  plateBGACountData <- nameCols(countData = plateBGACountData)

  # Fix gaps in data
  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  faecalData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  waterData <- fixMissingData(fixData = waterData, fixReps = waterReps)
  faecalData <- fixMissingData(fixData = faecalData, fixReps = faecalReps)

  plateBGACountData <- rbind(waterData, faecalData)

  # Add meta data on households to table
  plateBGACountData <- addHouseholds(countData = plateBGACountData,
                                     metaData = metaData)

  plateBGACountData$SampleID <-
    as.numeric(substr(plateBGACountData$SampleID, 5, 8))

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Salmonella", "E.coli")
  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  faecalData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  # Average plate counts per species and per sample accross all dilutions
  subPlateData <- averageCounts(countData = waterData, numDilutions = waterReps,
                                rawData = rawData, convertCols = convertCols,
                                bacteriaTypes = bacteriaTypes, campaignType = "HH")
  subPlateData <- rbind(subPlateData, averageCounts(countData = faecalData,
                                                    numDilutions = faecalReps,
                                                    rawData = rawData,
                                                    convertCols = convertCols,
                                                    bacteriaTypes = bacteriaTypes,
                                                    campaignType = "HH"))

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11,
                                                       length(subPlateData$Household))),
                                     subPlateData$SampleType,
                                     subPlateData$SampleID), ]
  sampleNumbers <- unique(subPlateData$SampleID)
  subPlateData$SampleID <- factor(subPlateData$SampleID,
                                  levels = sampleNumbers)

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
                     dataType = "HH", mainData = "Brilliant green agar plates",
                     subData = bacteriaTypes[iBacteria],
                     fillColours = c('chocolate4', 'skyblue'), weightOrVol = "B",
                     rowGraphs = 4, plotsDir = plotsDir,
                     plotTitle = paste("Household", stream, "BGA", bacteriaTypes[iBacteria]))
  }

}

