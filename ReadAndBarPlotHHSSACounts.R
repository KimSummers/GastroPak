# ReadAndBarPlotHHSSACounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Salmonella, one for Shigella
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHSSACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   stoolReps         - Number of stool dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/05/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSSACounts <- function(plateCountHHSSAFile, plotsDir, metaDataFile,
                                      rawData, stoolReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  plateSSACountData <- numericCounts(plateSSACountData)

  waterData <- plateSSACountData[plateSSACountData$`Sample Type` == "Water", ]
  stoolData <- plateSSACountData[plateSSACountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateSSACountData <- rbind(waterData, stoolData)

  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- addHouseholds(plateSSACountData, metaData)

  plateSSACountData$`Sample-ID` <-
    as.numeric(substr(plateSSACountData$`Sample-ID`, 5, 8))

  plateSSACountData <- nameCols(plateSSACountData)

  plateSSACountData <- plateSSACountData %>% relocate(Shigella,
                                                    .after = Salmonella)
  plateSSACountData <- plateSSACountData %>% relocate(E.coli,
                                                      .after = Shigella)

  bacteriaTypes <- c("Salmonella", "Shigella", "E.coli")
  convertCols <- c(which(colnames(plateSSACountData) == "Salmonella"),
                   which(colnames(plateSSACountData) == "Shigella"),
                   which(colnames(plateSSACountData) == "E.coli"))

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  stoolData <- plateSSACountData[plateSSACountData$SampleType == "Stool", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes, "HH"))

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household,
                                                       11, length(subPlateData$Household))),
                                     subPlateData$SampleType,
                                     subPlateData$SampleID), ]
  sampleNumbers <- unique(subPlateData$SampleID)
  subPlateData$SampleID <- factor(subPlateData$SampleID,
                                  levels = sampleNumbers)

  sampleCodes <- unique(subPlateData$SamplingCode)
  subPlateData$SamplingCode <- factor(subPlateData$SamplingCode,
                                      levels = sampleCodes)

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
    bactSubData <- subBacteriaData(subPlateData, bacteriaTypes[iBacteria],
                                   c("Water", "Stool"))

    for (iRow in (1:(nrow(bactSubData) / 3)))
    {
      meanVal <- mean(bactSubData$MeanCfu[bactSubData$SampleID ==
                                            bactSubData$SampleID[iRow * 3]],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[bactSubData$SampleID ==
                                        bactSubData$SampleID[iRow * 3]],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 3, 1:8], meanVal, sdVal,
                          bactSubData[iRow * 3, 10])

      if (iRow == 1)
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    colnames(sumData)[9:11] <- c("MeanCfu", "StdDev", "Household")

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "HH", "Salmonella-Shigella agar plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), "B", 4, plotsDir,
                     bacteriaTypes[iBacteria])
  }

}
