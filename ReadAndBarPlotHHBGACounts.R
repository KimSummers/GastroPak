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
#   stoolReps         - Number of stool dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  16/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHBGACounts <- function(plateCountFile, plotsDir, metaDataFile,
                                      rawData, stoolReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountHHBGAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  plateBGACountData <- numericCounts(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$`Sample Type` == "Water", ]
  stoolData <- plateBGACountData[plateBGACountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateBGACountData <- rbind(waterData, stoolData)

  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- nameCols(plateBGACountData)
  plateBGACountData <- addHouseholds(plateBGACountData, metaData)

  plateBGACountData$SampleID <-
    as.numeric(substr(plateBGACountData$SampleID, 5, 8))

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Salmonella", "E.coli")

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  stoolData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps, rawData,
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
  sampleType <- unique(subPlateData$SampleType)
  subPlateData$SampleType <- factor(subPlateData$SampleType,
                                    levels = sampleType)

  households <- unique(subPlateData$Household)
  subPlateData$Household <- factor(subPlateData$Household,
                                   levels = households)

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subBacteriaHHData(subPlateData, bacteriaTypes[iBacteria],
                                     c("Stool", "Water"))

    BarPlotGastroPak(bactSubData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "HH", "Brilliant green agar plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), "B", 4, plotsDir,
                     paste("Household BGA", bacteriaTypes[iBacteria]))
  }

}

