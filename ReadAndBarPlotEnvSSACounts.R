# ReadAndBarPlotSSACounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotSSACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/01/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotSSACounts <- function(plateCountSSAFile, plotsDir, metaDataFile,
                                    rawCounts, sedimentReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountSSAFile)
  metaData <- read.csv(metaDataFile)

  plateSSACountData <- numericCounts(plateSSACountData)

  waterData <- plateSSACountData[plateSSACountData$`Sample Type` == "Water", ]
  sedimentData <- plateSSACountData[plateSSACountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateSSACountData <- rbind(waterData, sedimentData)

  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- addRivers(plateSSACountData, metaData)
  plateSSACountData <- nameCols(plateSSACountData)

  bacteriaTypes <- c("Shigella", "Salmonella", "E.coli")

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  sedimentData <- plateSSACountData[plateSSACountData$SampleType == "Sediment", ]

  convertCols <- c(which(colnames(plateSSACountData) == "Shigella"),
                   which(colnames(plateSSACountData) == "Salmonella"),
                   which(colnames(plateSSACountData) == "E.coli"))

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes)
  subPlateData <- rbind(subPlateData, averageCounts(sedimentData, sedimentReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes))

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

  subPlateData <- subPlateData[order(subPlateData$SamplingSite), ]

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subPlateData[subPlateData$Bacteria ==
                                  bacteriaTypes[iBacteria], ]

    for (iRow in (1:(nrow(bactSubData) / 6)))
    {
      meanVal <- mean(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                             bactSubData$SamplingSite[iRow * 6]) &
                                            (bactSubData$SampleType == "Sediment")],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                         bactSubData$SamplingSite[iRow * 6]) &
                                        (bactSubData$SampleType == "Sediment")],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 6, c(1, 3, 6)], "Sediment", meanVal,
                          sdVal, bactSubData[iRow * 6, 10:11])

      colnames(sumDataRow)[4] <- "SampleType"

      if (iRow == 1)
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

      meanVal <- mean(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                             bactSubData$SamplingSite[iRow * 6]) &
                                            (bactSubData$SampleType == "Water")],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                         bactSubData$SamplingSite[iRow * 6]) &
                                        (bactSubData$SampleType == "Water")],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 6, c(1, 3, 6)], "Water", meanVal,
                          sdVal, bactSubData[iRow * 6, 10:11])

      colnames(sumDataRow)[4] <- "SampleType"
      sumData <- rbind(sumData, sumDataRow)
    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    sumData <- sumData[order(sumData$SamplingSite), ]
    colnames(sumData)[5:8] <- c("MeanCfu", "StdDev", "River", "Location")

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "Env", "Salmonella-Shigella agar plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), "B", 3, plotsDir,
                     bacteriaTypes[iBacteria])
  }

}
