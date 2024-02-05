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

  # Read in data and put it in dataframes
  plateCSCountData <- read_csv(plateCountHHCSFile)
  plateSSACountData <- read_csv(plateCountHHSSAFile)
  plateBGACountData <- read_csv(plateCountHHBGAFile)

  # For both the qPCR and meta data the first line has information not needed
  # by this code that should be removed
  qPCRData <- read_csv(qPCRFile)
  qPCRData <- qPCRData[!is.na(qPCRData$Target), ]
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[!is.na(metaData$`Sample code`), ]

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
  plateCSCountData <- nameCols(plateCSCountData)
  plateBGACountData <- nameCols(plateBGACountData)
  plateSSACountData <- nameCols(plateSSACountData)

  # Add an extra column for the media type to each count table
  plateCSCountData <- cbind(plateCSCountData, rep("ChromoSelect", nrow(plateCSCountData)))
  colnames(plateCSCountData)[ncol(plateCSCountData)] <- "Media"

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

  selCSCols <- c(which(colnames(plateCSCountData) == "SampleID"),
                 which(colnames(plateCSCountData) == "SamplingSite"),
                 which(colnames(plateCSCountData) == "SampleType"),
                 which(colnames(plateCSCountData) == "Replicate"),
                 which(colnames(plateCSCountData) == "Dilution"),
                 which(colnames(plateCSCountData) == "E.coli"),
                 which(colnames(plateCSCountData) == "Media"))
  selBGACols <- c(which(colnames(plateBGACountData) == "SampleID"),
                  which(colnames(plateBGACountData) == "SamplingSite"),
                  which(colnames(plateBGACountData) == "SampleType"),
                  which(colnames(plateBGACountData) == "Replicate"),
                  which(colnames(plateBGACountData) == "Dilution"),
                  which(colnames(plateBGACountData) == "E.coli"),
                  which(colnames(plateBGACountData) == "Media"))
  selSSACols <- c(which(colnames(plateSSACountData) == "SampleID"),
                  which(colnames(plateSSACountData) == "SamplingSite"),
                  which(colnames(plateSSACountData) == "SampleType"),
                  which(colnames(plateSSACountData) == "Replicate"),
                  which(colnames(plateSSACountData) == "Dilution"),
                  which(colnames(plateSSACountData) == "E.coli"),
                  which(colnames(plateSSACountData) == "Media"))
  # Now combine the sample data and e. coli counts across the plates into one new table
  plateCountCombData <- rbind(plateCSCountData[selCSCols],
                              plateSSACountData[selSSACols],
                              plateBGACountData[selBGACols])

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateCountCombData <- rbind(waterData, stoolData)

  plateCountCombData <- addHouseholds(plateCountCombData, metaData)

  colnames(qPCRData)[which(colnames(qPCRData) == "Sample")] <- "SampleID"
  qPCRData <- addHouseholds(qPCRData, metaData)

  # Convert the sample ids into the numeric part for sorting
  plateCountCombData$SampleID <-
    as.numeric(substr(plateCountCombData$SampleID, 5, 8))

  mediaTypes <- c("ChromoSelect", "SSA", "BGA", "UK qPCR", "Pak qPCR")

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  stoolData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateCountCombData) == "E.coli")

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                "E.coli", "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps,
                                                    rawData, convertCols,
                                                    "E.coli", "HH"))

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
    bactSubData <- subSampleHHData(subPlateData, sampleType[iSampleType], mediaTypes)

    if (sampleType[iType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(sumData, "", "Media", "HH", "E coli counts", sampleTypes[iType],
                     c('darkblue', 'pink',  'lightgreen', 'gold', 'lightgrey'),
                     "B", 5, plotsDir, sampleTypes[iType])
  }

}

