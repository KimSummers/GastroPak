# ReadAndBarPlotHHEColiCountsComp
#
# Read in plate count and qPCR data and barchart plot it
# One bar chart for the counts on ChromoSelect, one for SSA, one for BGA
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
#   waterReps           - Number of water dilutions
#   feacalReps          - Number of faecal dilutions
#   plotsDir            - directory to store plots in

# Version    Author       Date      Affiliation
# 1.00       J K Summers  28/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHEColiCountsComp <- function(plateCountHHCSFile, plateCountHHSSAFile,
                                            plateCountHHBGAFile, qPCRFile,
                                            metaDataFile, waterReps, feacalReps,
                                            plotsDir) {

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

  qPCRData <- nameCols(qPCRData)

  # Only want E. coli faecal samples
  qPCRData <- qPCRData[(qPCRData$Target == "E. coli") & (qPCRData$Season == "Household"), ]

  # Ensure all the counts data is numeric
  plateCSCountData <- numericCounts(plateCSCountData)
  plateSSACountData <- numericCounts(plateSSACountData)
  plateBGACountData <- numericCounts(plateBGACountData)

  # collate count data by species
  plateCSCountData <- coloursToSpeciesCS(plateCSCountData)
  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)

  # Ensure columns are correctly named
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

  # Fill in missing data
  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  faecalData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  faecalData <- fixMissingData(faecalData, faecalReps)

  plateCountCombData <- rbind(waterData, faecalData)

  # add metadata on households to plate count and qPCR data
  plateCountCombData <- addHouseholds(plateCountCombData, metaData)

  if (nrow(qPCRData) > 0)
  {
    qPCRData <- addHouseholds(qPCRData, metaData)
  }

  # Convert the sample ids into the numeric part for sorting
  plateCountCombData$SampleID <-
    as.numeric(substr(plateCountCombData$SampleID, 5, 8))

  # Amend the column names
  plateCountCombData <- nameCols(plateCountCombData)
  colnames(plateCountCombData)[which(colnames(plateCountCombData) == "E.coli")] <- "CfuCount"

  mediaTypes <- c("ChromoSelect", "SSA", "BGA", "UK qPCR", "Pak qPCR")

  waterData <- plateCountCombData[plateCountCombData$SampleType == "Water", ]
  faecalData <- plateCountCombData[plateCountCombData$SampleType == "Stool", ]

  convertCols <- which(colnames(plateCountCombData) == "CfuCount")

  # Average colony counts across dilutions
  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                "E.coli", "HH", TRUE)
  subPlateData <- rbind(subPlateData, averageCounts(faecalData, faecalReps,
                                                    rawData, convertCols,
                                                    "E.coli", "HH", TRUE))

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
  qPCRSubData <- cbind(qPCRData[, qPCRFirstSelCols], rep("E.coli", nrow(qPCRData)),
                       rep(0, nrow(qPCRData)), rep(0, nrow(qPCRData)),
                       qPCRData[, qPCRSecondSelCols])

  # combine colony counts and qPCR data
  if (nrow(qPCRData) > 0)
  {
    colnames(qPCRSubData) <- colnames(subPlateData)
    subPlateData <- rbind(subPlateData, qPCRSubData)
  }
  
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
    bactSubData <- subCompHHData(countData = subPlateData, 
                                 sampleType = sampleType[iSampleType])

    if (sampleType[iType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(graphData = bactSubData, xAxisLabels = "", plotType = "Media", 
                     dataType = "HH", mainData = "E coli counts", 
                     subData = sampleTypes[iType], 
                     fillColours = c('darkblue', 'pink',  'lightgreen', 'gold', 'lightgrey'),
                     weightOrVol = measureType, rowGraphs = 4, plotsDir = plotsDir, 
                     plotTitle = paste("Household E. coli", sampleTypes[iType]))
  }

}

