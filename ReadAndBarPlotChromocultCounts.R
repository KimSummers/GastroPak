# ReadAndBarPlotChromoSelectCounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotChromoSelectCounts
#
# inputs
# 	plateCountCSFile  - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/12/22  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotChromoSelectCounts <- function(plateCountEnvCSFile, plotsDir,
                                             metaDataFile, rawData,
                                             sedimentReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountEnvCSFile)
  metaData <- read.csv(metaDataFile)

  plateCSCountData <- numericCounts(plateCSCountData)

  waterData <- plateCSCountData[plateCSCountData$`Sample Type` == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateCSCountData <- rbind(waterData, sedimentData)

  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateCSCountData <- addRivers(plateCSCountData, metaData)

  plateCSCountData <- nameCols(plateCSCountData)

  plateCSCountData <- plateCSCountData %>% relocate(Klebsiella,
                                                    .after = coliforms)

  convertCols <- c(which(colnames(plateCSCountData) == "E.coli"),
                   which(colnames(plateCSCountData) == "coliforms"),
                   which(colnames(plateCSCountData) == "Klebsiella"))

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")

  waterData <- plateCSCountData[plateCSCountData$`Sample Type` == "Water"]
  sedimentData <- plateCSCountData[plateCSCountData$`Sample Type` == "Sediment"]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "Env")
  subPlateData <- rbind(subPlateData, averageCounts(sedimentData, sedimentReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes, "Env"))

  colnames(subPlateData)[2] <- "SampleID"

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

  subPlateData <- subPlateData[order(subPlateData$SampleID), ]

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
      sumDataRow <- cbind(bactSubData[iRow * 6, c(1,3,6)], "Sediment", meanVal,
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
      sumDataRow <- cbind(bactSubData[iRow * 6, c(1,3,6)], "Water", meanVal,
                          sdVal, bactSubData[iRow * 6, 10:11])

      colnames(sumDataRow)[4] <- "SampleType"
      sumData <- rbind(sumData, sumDataRow)
    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    sumData <- sumData[order(sumData$SamplingSite), ]
    colnames(sumData)[5:8] <- c("MeanCfu", "StdDev", "River", "Location")

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Sample",
                     "Chromocult plates", bacteriaTypes[iBacteria],
                     c('chocolate4', 'skyblue'), plotsDir, bacteriaTypes[iBacteria])
  }

}
