# ReadAndBarPlotHHSalmonellaCountsComp
#
# Read in plate counts and qPCR data and bar chart plot it
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
#   waterReps           - Number of water dilutions
#   faecalReps          - Number of faecal dilutions
#   plotsDir            - directory to store plots in

# Version    Author       Date      Affiliation
# 1.00       J K Summers  03/08/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSalmonellaCountsComp <- function(plateCountHHSSAFile,
                                                 plateCountHHBGAFile, qPCRFile,
                                                 metaDataFile, waterReps,
                                                 faecalReps, plotsDir) {

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

  qPCRData <- nameCols(qPCRData)
  
  # Only want Salmonella faecal samples
  qPCRData <- qPCRData[(qPCRData$Target == "Salmonella") & (qPCRData$Season == "Household"), ]

  # Ensure all the counts data is numeric
  plateSSACountData <- numericCounts(plateSSACountData)
  plateBGACountData <- numericCounts(plateBGACountData)

  # collate count data by species
  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)

  # Ensure columns are correctly named
  plateBGACountData <- nameCols(plateBGACountData)
  plateSSACountData <- nameCols(plateSSACountData)
  qPCRData <- nameCols(qPCRData)

  # Add an extra column for the media type to each count table
  plateSSACountData <- cbind(plateSSACountData, rep("SSA", nrow(plateSSACountData)))
  colnames(plateSSACountData)[ncol(plateSSACountData)] <- "Media"

  plateBGACountData <- cbind(plateBGACountData, rep("BGA", nrow(plateBGACountData)))
  colnames(plateBGACountData)[ncol(plateBGACountData)] <- "Media"

  if (nrow(qPCRData) > 0)
  {
    qPCRData <- cbind(qPCRData, paste(qPCRData$TestCountry, "qPCR"))
    colnames(qPCRData)[ncol(qPCRData)] <- "Media"
  }

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

  # Fill in missing data
  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  faecalData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateCountCombData <- rbind(waterData, faecalData)

  # add metadata on households to plate count and qPCR data
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
  faecalData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateCountCombData) == "CfuCount")

  # Average colony counts across dilutions
  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                "Salmonella", "HH", TRUE)
  subPlateData <- rbind(subPlateData, averageCounts(faecalData, faecalReps,
                                                    rawData, convertCols,
                                                    "Salmonella", "HH", TRUE))

  # Select the columns we need from the qPCR data
  qPCRData <- nameCols(qPCRData)
  qPCRFirstSelCols <- c(which(colnames(qPCRData) == "SampleID"),
                        which(colnames(qPCRData) == "SamplingSite"),
                        which(colnames(qPCRData) == "SampleType"),
                        which(colnames(qPCRData) == "Replicate"))

  qPCRSecondSelCols <- c(which(colnames(qPCRData) == "MeanCfu"),
                         which(colnames(qPCRData) == "Household"),
                         which(colnames(qPCRData) == "Media"))

  qPCRData$SampleID <- as.numeric(substr(qPCRData$SampleID, 5, nchar(qPCRData$SampleID)))
  qPCRSubData <- cbind(qPCRData[, qPCRFirstSelCols], rep("Salmonella", nrow(qPCRData)),
                       rep(0, nrow(qPCRData)), rep(0, nrow(qPCRData)),
                       qPCRData[, qPCRSecondSelCols])

  colnames(qPCRSubData) <- colnames(subPlateData)

  # combine colony counts and qPCR data
  subPlateData <- rbind(subPlateData, qPCRSubData)

  subPlateData$MeanCfu <- as.numeric(subPlateData$MeanCfu)

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11, 
                                                       length(subPlateData$Household))),
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

  sampleTypes <- c("Faecal", "Water")

  for (iType in 1:length(sampleTypes))
  {
    bactSubData <- subCompHHData(subPlateData, sampleType[iSampleType])

    if (sampleType[iType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(bactSubData, "", "Media", "HH", "Salmonella counts", sampleType[iType],
                     fillColours = c('palevioletred', 'black',  'gold', 'lightgrey'),
                     measureType, 4, plotsDir, paste("Household Salmonella",
                                               sampleTypes[iType]))
  }

}
