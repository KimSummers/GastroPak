# ReadAndBarPlotHHBGACounts
#
# Read in plate count data and barchart plot it
# One bar chart for the E. coli counts and one for Salmonella
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHBGACounts
#
# inputs
# 	plateCountFile    - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   faecalReps        - Number of faecal dilutions
#   waterReps         - Number of water dilutions
#   stream            - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  16/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHBGACounts <- function(plateCountFile, plotsDir, metaDataFile,
                                      rawData, faecalReps, waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountHHBGAFile)
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

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  faecalData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(faecalData, faecalReps, rawData,
                                                    convertCols, bacteriaTypes, "HH"))

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

