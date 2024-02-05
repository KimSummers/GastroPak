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
                                                 metaDataFile, waterReps,
                                                 stoolReps, plotsDir) {

  # First load any libraries needed
  library(tidyverse)

  # Read in data and put it in dataframes
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  plateBGACountData <- read_csv(plateCountHHBGAFile)

  # For both the qPCR and meta data the first line has information not needed
  # by this code that should be removed
  qPCRData <- read_csv(qPCRFile)
  qPCRData <- qPCRData[!is.na(qPCRData$Target), ]
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

  # Only want Salmonella stool samples
  qPCRData <- qPCRData[(qPCRData$Target == "Salmonella") & (qPCRData$Season == "Household"), ]

  # Ensure all the counts data is numeric
  plateSSACountData <- numericCounts(plateSSACountData)
  plateBGACountData <- numericCounts(plateBGACountData)

  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)

  # Add an extra column for the media type to each count table
  plateSSACountData <- cbind(plateSSACountData, rep("SSA", nrow(plateSSACountData)))
  colnames(plateSSACountData)[ncol(plateSSACountData)] <- "Media"

  plateBGACountData <- cbind(plateBGACountData, rep("BGA", nrow(plateBGACountData)))
  colnames(plateBGACountData)[ncol(plateBGACountData)] <- "Media"

  qPCRData <- nameCols(qPCRData)

  if (nrow(qPCRData) > 0)
  {
    qPCRData <- cbind(qPCRData, paste(qPCRData$TestCountry, "qPCR"))
    colnames(qPCRData)[ncol(qPCRData)] <- "Media"
  }

  plateSSACountData <- nameCols(plateSSACountData)
  plateBGACountData <- nameCols(plateBGACountData)

  selSSACols <- c(which(colnames(plateSSACountData) == "SampleID"),
                  which(colnames(plateSSACountData) == "SamplingSite"),
                  which(colnames(plateSSACountData) == "SampleType"),
                  which(colnames(plateSSACountData) == "Replicate"),
                  which(colnames(plateSSACountData) == "Dilution"),
                  which(colnames(plateSSACountData) == "Salmonella"),
                  which(colnames(plateSSACountData) == "Media"))
  selBGACols <- c(which(colnames(plateBGACountData) == "SampleID"),
                  which(colnames(plateBGACountData) == "SamplingSite"),
                  which(colnames(plateBGACountData) == "SampleType"),
                  which(colnames(plateBGACountData) == "Replicate"),
                  which(colnames(plateBGACountData) == "Dilution"),
                  which(colnames(plateBGACountData) == "Salmonella"),
                  which(colnames(plateBGACountData) == "Media"))

  # Now combine the sample data and salmonella counts across the plates into one new table
  plateCountCombData <- rbind(plateSSACountData[selSSACols],
                              plateBGACountData[selBGACols])

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateCountCombData <- rbind(waterData, stoolData)
  plateCountCombData <- addHouseholds(plateCountCombData, metaData)

  qPCRData <- addHouseholds(qPCRData, metaData)

  # Convert the sample ids into the numeric part for sorting
  plateCountCombData$SampleID <-
    as.numeric(substr(plateCountCombData$SampleID, 5, 8))

  # Amend the column names
  plateCountCombData <- nameCols(plateCountCombData)
  colnames(plateCountCombData)[which(colnames(plateCountCombData) == "Salmonella")] <- "CfuCount"

  mediaTypes <- c("SSA", "BGA", "UK qPCR", "Pak qPCR")

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateCountCombData) == "CfuCount")

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                "CfuCount", "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps,
                                                    rawData, convertCols,
                                                    "CfuCount", "HH"))

  qPCRData <- nameCols(qPCRData)

  qPCRData$SampleID <- as.numeric(substr(qPCRData$SampleID, 5, length(qPCRData$SampleID)))
  qPCRSubData <- cbind(qPCRData[, c(1, 5, 3, 7)], rep(0, nrow(qPCRData)),
                       rep(0, nrow(qPCRData)), qPCRData[, 20:22])

  colnames(qPCRSubData) <- colnames(subPlateData)

  subPlateData <- rbind(subPlateData, qPCRSubData)

  subPlateData$MeanCfu <- as.numeric(subPlateData$MeanCfu)

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
    bactSubData <- subSampleHHData(subPlateData, sampleType[iSampleType], mediaTypes)

    if (sampleType[iType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    if (iType == 1)
    {
      BarPlotGastroPak(bactSubData, "", "Media", "HH", "Salmonella counts", sampleType[iType],
                       fillColours = c('palevioletred', 'black',  'gold', 'lightgrey'),
                       "W", 4, plotsDir, sampleTypes[iType])
    }else
    {
      BarPlotGastroPak(sumData, "", "Media", "HH", "E coli counts", "Water",
                       fillColours = c('palevioletred', 'black',  'gold', 'lightgrey'),
                       "V", 4, plotsDir, paste("Household", sampleTypes[iType]))
    }

  }

}

