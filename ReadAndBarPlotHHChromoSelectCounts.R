# ReadAndBarPlotHHChromoSelectCounts
#
# Read in plate count data and bar chart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by household and order by household, sample type and sample id
#
# file ReadAndBarPlotHHChromoSelectCounts
#
# inputs
# 	plateCountHHCSFile  - file containing plate count data
#   plotsDir            - directory to store plots in
#   metaDataFile        - file containing river meta data
#   rawData             - True for absolute counts, false for cfu / ml
#   faecalReps          - Number of faecal dilutions
#   waterReps           - Number of water dilutions
#   stream              - Up, mid or down stream

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/12/22  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHChromoSelectCounts <- function(plateCountHHCSFile, plotsDir,
                                               metaDataFile, rawData, faecalReps,
                                               waterReps, stream) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountHHCSFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  plateCSCountData <- numericCounts(plateCSCountData)

  waterData <- plateCSCountData[plateCSCountData$`Sample Type` == "Water", ]
  faecalData <- plateCSCountData[plateCSCountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateCSCountData <- rbind(waterData, faecalData)

  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- nameCols(plateCSCountData)
  plateCSCountData <- addHouseholds(plateCSCountData, metaData)

  plateCSCountData$SampleID <-
    as.numeric(substr(plateCSCountData$SampleID, 5, 8))

  plateCSCountData <- plateCSCountData %>% relocate(Klebsiella,
                                                    .after = coliforms)

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  faecalData <- plateCSCountData[plateCSCountData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateCSCountData) %in% bacteriaTypes)

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(faecalData, faecalReps, rawData,
                                                    convertCols, bacteriaTypes,
                                                    "HH"))

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
  subPlateData$Household <- factor(subPlateData$Household, levels = households)

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subBacteriaHHData(countData = subPlateData,
                                     bacteriaType = bacteriaTypes[iBacteria],
                                     sampleTypes = c("Faecal", "Water"))

    BarPlotGastroPak(graphData = bactSubData, xAxisLabels = NULL, plotType = "Sample",
                     dataType = "HH", mainData = "ChromoSelect plates",
                     subData = bacteriaTypes[iBacteria],
                     fillColours = c('chocolate4', 'skyblue'), weightOrVol = "B",
                     rowGraphs = 4, plotsDir = plotsDir,
                     plotTitle = paste("Household", stream, "CS", bacteriaTypes[iBacteria]))
  }

}

