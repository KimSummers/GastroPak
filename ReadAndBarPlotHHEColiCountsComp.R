# ReadAndBarPlotHHEColiCountsComp
#
# Read in plate count and qPCR data and barchart plot it
# One bar char for the counts on ChromoSelect, one for SSA, one for BGA
# and one for the qPCR
# Group counts by household and order by household, sample type and sample id
#
# file ReadAndBarPlotHHEColiCountsComp
#
# inputs
# 	plateCountHHCSFile  - file containing plate count data on ChromoSelect
# 	plateCountHHSSAFile - file containing plate count data on SSA
# 	plateCountHHBGAFile - file containing plate count data on BGA
# 	qPCRFile            - file containing qPCR data
#   metaDataFile        - file containing meta data about the samples
#   rawData           - True for absolute counts, false for cfu / ml
#   stoolReps         - Number of stool dilutions
#   waterReps         - Number of water dilutions
#   plotsDir            - directory to store plots in

# Version    Author       Date      Affiliation
# 1.00       J K Summers  28/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHEColiCountsComp <- function(plateCountHHCSFile, plateCountHHSSAFile,
                                            plateCountHHBGAFile, qPCRFile,
                                            waterReps, stoolReps,
                                            metaDataFile, plotsDir) {

  # First load any libraries needed
  library(tidyverse)
  library(scales)
  library(ggplot2)

  # Read in data and put it in dataframes
  plateCSCountData <- read_csv(plateCountHHCSFile)
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  plateBGACountData <- read_csv(plateCountHHBGAFile)

  # For both the qPCR and meta data the first line has information not needed
  # by this code that should be removed
  qPCRData <- read_csv(qPCRFile)
  qPCRData <- qPCRData[2:nrow(qPCRData), ]
  qPCRData <- qPCRData[qPCRData$Sample %in% samplesIDs, ]
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  # Ensure all the counts data is numeric
  plateCSCountData <- numericCounts(plateCSCountData)
  plateCSCountData <- plateCSCountData[, 1:(ncol(plateCSCountData) - 3)]

  plateSSACountData <- numericCounts(plateSSACountData)
  plateBGACountData <- numericCounts(plateBGACountData)

  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)

  # Only want E. coli stool samples
  qPCRData <- qPCRData[(qPCRData$Target == "E. coli") & (qPCRData$Season == "Household"), ]

  # Add an extra column for the media type to each count table
  plateCCCountData <- cbind(plateCCCountData, rep("ChromoSelect", nrow(plateCCCountData)))
  colnames(plateCCCountData)[19] <- "Media"

  plateSSACountData <- cbind(plateSSACountData, rep("SSA", nrow(plateSSACountData)))
  colnames(plateSSACountData)[9] <- colnames(plateCCCountData)[18]
  colnames(plateSSACountData)[15] <- "Media"

  plateBGACountData <- cbind(plateBGACountData, rep("BGA", nrow(plateBGACountData)))
  colnames(plateBGACountData)[9] <- colnames(plateCCCountData)[18]
  colnames(plateBGACountData)[12] <- "Media"

  # qPCRData <- cbind(qPCRData, paste(qPCRData$Country, "qPCR"))
  qPCRData <- cbind(qPCRData, "qPCR")
  colnames(qPCRData)[ncol(qPCRData)] <- "Media"

  # Now combine the sample data and e. coli counts across the plates into one new table
  plateCountCombData <- rbind(plateCCCountData[c(2:6, 18:19)],
                              plateSSACountData[c(2:6, 9, 15)],
                              plateBGACountData[c(2:6, 9, 12)])

  waterData <- plateCountCombData[plateCountCombData$`Sample Type` == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateCountCombData <- rbind(waterData, stoolData)

  plateCountCombData <- addHousehold(plateCountCombData, metaData)

  colnames(qPCRData)[which(colnames(qPCRData) == "Sample")] <- "Sample-ID"
  qPCRData <- addHousehold(qPCRData, metaData)

  # Convert the sample ids into the numeric part for sorting
  plateCountCombData$`Sample-ID` <-
    as.numeric(substr(plateCountCombData$`Sample-ID`, 5, 8))

  # Amend the column names
  colnames(plateCountCombData)[2] <- "SamplingSite"
  colnames(plateCountCombData)[3] <- "SampleType"
  colnames(plateCountCombData)[6] <- "CfuCount"
  colnames(plateCountCombData)[8] <- "Household"

  mediaTypes <- c("ChromoSelect", "SSA", "BGA", "UK qPCR", "Pak qPCR")

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, 6, "E.coli")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps,
                                                    rawData, 6, "E.coli"))

  colnames(qPCRData)[1] <- "SampleID"
  colnames(qPCRData)[2] <- "SampleType"
  colnames(qPCRData)[5] <- "SamplingSite"
  colnames(qPCRData)[20] <- "MeanCfu"
  colnames(qPCRData)[22] <- "Household"

  qPCRData$SampleID <- as.numeric(substr(qPCRData$SampleID, 5, length(qPCRData$SampleID)))
  qPCRSubData <- cbind(qPCRData[, c(1, 5, 3, 7)], rep(0, nrow(qPCRData)),
                       rep(0, nrow(qPCRData)), qPCRData[, 20:22])

  colnames(qPCRSubData) <- colnames(subPlateData)

  subPlateData <- rbind(subPlateData, qPCRSubData)

  subPlateData$MeanCfu <- as.numeric(subPlateData$MeanCfu)
  colnames(subPlateData)[1] <- "SampleID"

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11, length(subPlateData$Household))),
                                     subPlateData$SampleType,
                                     subPlateData$SampleID), ]
  sampleNumbers <- unique(subPlateData$SampleID)
  subPlateData$SampleID <- factor(subPlateData$SampleID,
                                     levels = sampleNumbers)

  mediaTypes <- unique(subPlateData$Media)
  subPlateData$Media <- factor(subPlateData$Media, levels = mediaTypes)

  sampleSites <- unique(subPlateData$SamplingSite)
  subPlateData$SamplingSite <- factor(subPlateData$SamplingSite,
                                      levels = sampleSites)
  sampleType <- unique(subPlateData$SampleType)
  subPlateData$SampleType <- factor(subPlateData$SampleType,
                                    levels = sampleType)

  households <- unique(subPlateData$Household)
  subPlateData$Household <- factor(subPlateData$Household,
                                   levels = households)

  sampleTypes <- c("Stool", "Water")

  for (iType in 1:length(sampleTypes))
  {
    bactSubData <- subPlateData[subPlateData$SampleType ==
                                  sampleTypes[iType], ]

    for (iRow in (1:(nrow(bactSubData) / 3)))
    {
      meanVal <- mean(bactSubData$MeanCfu[(bactSubData$SampleID ==
                                            bactSubData$SampleID[iRow * 3]) &
                                            (bactSubData$Media ==
                                               bactSubData$Media[iRow * 3])],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[(bactSubData$SampleID ==
                                        bactSubData$SampleID[iRow * 3]) &
                                        (bactSubData$Media ==
                                           bactSubData$Media[iRow * 3])],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 3, 1:7], meanVal, sdVal,
                          bactSubData[iRow * 3, 8:9])

      if (iRow == 1)
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    sumData$sdVal[sumData$sdVal == 0] <- NaN
    sumData$sdVal[(sumData$meanVal - sumData$sdVal) < 1] <- NaN

    BarPlotGastroPak(sumData, "", "Media", "HH", "E coli counts", sampleTypes[iType],
                     c('darkblue', 'pink',  'lightgreen', 'gold', 'lightgrey'),
                     "B", 5, plotsDir, sampleTypes[iType])
  }

}

