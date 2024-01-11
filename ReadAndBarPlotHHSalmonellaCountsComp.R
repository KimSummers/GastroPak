# ReadAndBarPlotHHSalmonellaCountsComp
#
# Read in plate counts and qPCR data and barchart plot it
# One bar char for the counts on SSA, one for BGA and one for the qPCR
# Group counts by household and order by household, sample type and sample id
#
# file ReadAndBarPlotHHSalmonellaCountsComp
#
# inputs
# 	plateCountHHSSAFile - file containing plate count data on SSA
# 	plateCountHHBGAFile - file containing plate count data on BGA
# 	qPCRFile            - file containing qPCR data
#   metaDataFile        - file containing meta data about the samples
#   stoolReps           - Number of stool dilutions
#   waterReps           - Number of water dilutions
#   plotsDir            - directory to store plots in

# Version    Author       Date      Affiliation
# 1.00       J K Summers  03/08/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSalmonellaCountsComp <- function(plateCountHHSSAFile,
                                            plateCountHHBGAFile, qPCRFile,
                                            metaDataFile, waterReps, stoolReps,
                                            plotsDir) {

  # First load any libraries needed
  library(tidyverse)

  # Read in data and put it in dataframes
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  plateBGACountData <- read_csv(plateCountHHBGAFile)

  # For both the qPCR and meta data the first line has information not needed
  # by this code that should be removed
  qPCRData <- read_csv(qPCRFile)
  qPCRData <- qPCRData[2:nrow(qPCRData), ]
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  # Only want Salmonella stool samples
  qPCRData <- qPCRData[(qPCRData$Target == "Salmonella") & (qPCRData$Season == "Household"), ]

  # Ensure all the counts data is numeric
  plateSSACountData <- numericCounts(plateSSACountData)
  plateBGACountData <- numericCounts(plateBGACountData)

  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)

  # Add an extra column for the media type to each count table
  plateSSACountData <- cbind(plateSSACountData, rep("SSA", nrow(plateSSACountData)))
  colnames(plateSSACountData)[15] <- "Media"

  plateBGACountData <- cbind(plateBGACountData, rep("BGA", nrow(plateBGACountData)))
  colnames(plateBGACountData)[12] <- "Media"

  qPCRData <- cbind(qPCRData, paste(qPCRData$Country, "qPCR"))
  colnames(qPCRData)[ncol(qPCRData)] <- "Media"

  # Now combine the sample data and salmonella counts across the plates into one new table
  plateCountCombData <- rbind(plateSSACountData[c(2:6, 9, 15)],
                              plateBGACountData[c(2:6, 9, 12)])

  waterData <- plateCountCombData[plateCountCombData$`Sample Type` == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateCountCombData <- addHouseholds(plateCountCombData, metaData)

  colnames(qPCRData)[which(colnames(qPCRData) == "Sample")] <- "Sample-ID"
  qPCRData <- addHouseholds(qPCRData, metaData)

  # Convert the sample ids into the numeric part for sorting
  plateCountCombData$`Sample-ID` <-
    as.numeric(substr(plateCountCombData$`Sample-ID`, 5, 8))

  # Amend the column names
  colnames(plateCountCombData)[2] <- "SamplingSite"
  colnames(plateCountCombData)[3] <- "SampleType"
  colnames(plateCountCombData)[6] <- "CfuCount"
  colnames(plateCountCombData)[8] <- "Household"

  mediaTypes <- c("SSA", "BGA", "UK qPCR", "Pak qPCR")

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  subPlateData <- averageCounts(waterData)
  subPlateData <- rbind(subPlateData, averageCounts(stoolData))

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

    if (iType == 1)
    {
      BarPlotGastroPak(sumData, "", "Media", "HH", "Salmonella counts", "Stool",
                       fillColours = c('palevioletred', 'black',  'gold', 'lightgrey'),
                       "W", 4, plotsDir, sampleTypes[iType])
    }else
    {
      BarPlotGastroPak(sumData, "", "Media", "HH", "E coli counts", "Water",
                       fillColours = c('palevioletred', 'black',  'gold', 'lightgrey'),
                       "V", 4, plotsDir, sampleTypes[iType])
    }

  }

}

