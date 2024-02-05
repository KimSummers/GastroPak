# ReadAndBarPlotEnvEColiCountsComp
#
# Read in plate count and qPCR data and barchart plot it
# One bar char for the counts on ChromoSelect, one for SSA, one for BGA
# and one for the qPCR
# Group counts by river and order by sample id
#
# file ReadAndBarPlotEnvEColiCountsComp
#
# inputs
# 	plateCountEnvCSFile   - file containing plate count data on ChromoSelect
# 	plateCountEnvSSAFile  - file containing plate count data on SSA
# 	plateCountEnvBGAFile  - file containing plate count data on BGA
# 	qPCRFile              - file containing qPCR data
#   metaDataFile          - file containing meta data about the samples
#   rawData               - True for absolute counts, false for cfu / ml
#   sedimentReps          - Number of sediment dilutions
#   waterReps             - Number of water dilutions
#   plotsDir              - directory to store plots in

# Version    Author       Date      Affiliation
# 1.00       J K Summers  11/12/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvEColiCountsComp <- function(plateCountEnvCSFile, plateCountEnvSSAFile,
                                             plateCountEnvBGAFile, qPCRFile,
                                             metaDataFile, rawData, sedimentReps,
                                             waterReps, plotsDir) {

  # First load any libraries needed
  library(tidyverse)

  # Read in data and put it in dataframes
  plateCSCountData <- read_csv(plateCountEnvCSFile)
  plateSSACountData <- read_csv(plateCountEnvSSAFile)
  plateBGACountData <- read_csv(plateCountEnvBGAFile)
  metaData <- read.csv(metaDataFile)

  # For both the qPCR and meta data the first line has information not needed
  # by this code that should be removed
  qPCRData <- read_csv(qPCRFile)
  # Only want E. coli stool samples
  qPCRData <- qPCRData[!is.na(qPCRData$Target), ]

  qPCRData <- nameCols(qPCRData)
  qPCRData <- qPCRData[(qPCRData$Target == "E. coli") & (qPCRData$Season == "Wet 22")
                       & (qPCRData$SampleType == "Sediment") &
                         (qPCRData$TestCountry == "UK"), ]

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

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")

  convertCols <- which(colnames(plateCSCountData) %in% bacteriaTypes)

  waterData <- plateCSCountData[plateCSCountData$SampleType == "Water", ]
  sedimentData <- plateCSCountData[plateCSCountData$SampleType == "Sediment", ]

  subPlateCSData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                  bacteriaTypes, "Env")
  subPlateCSData <- rbind(subPlateCSData, averageCounts(sedimentData, sedimentReps,
                                                        rawData, convertCols,
                                                        bacteriaTypes, "Env"))

  subPlateCSData <- subPlateCSData[subPlateCSData$Bacteria == "E.coli", ]

  # Add an extra column for the media type to each count table
  subPlateCSData <- cbind(subPlateCSData, rep("ChromoSelect", nrow(subPlateCSData)))
  colnames(subPlateCSData)[ncol(subPlateCSData)] <- "Media"

  # SSA data
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

  for (iBacteria in 1:length(bacteriaTypes))
  {

    if (iBacteria == 1)
    {
      convertCols <- which(colnames(plateSSACountData) == bacteriaTypes[iBacteria])
    }else
    {
      convertCols <- c(convertCols,
                       which(colnames(plateSSACountData) == bacteriaTypes[iBacteria]))
    }

  }

  subPlateSSAData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                   bacteriaTypes)
  subPlateSSAData <- rbind(subPlateSSAData, averageCounts(sedimentData, sedimentReps,
                                                        rawData, convertCols,
                                                        bacteriaTypes))

  subPlateSSAData <- subPlateSSAData[subPlateSSAData$Bacteria == "E.coli", ]

  # Add an extra column for the media type to each count table
  subPlateSSAData <- cbind(subPlateSSAData, rep("SSA", nrow(subPlateSSAData)))
  colnames(subPlateSSAData)[ncol(subPlateSSAData)] <- "Media"

  # BGA
  plateBGACountData <- numericCounts(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$`Sample Type` == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateBGACountData <- rbind(waterData, sedimentData)

  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- addRivers(plateBGACountData, metaData)

  colnames(plateBGACountData)[1] <- "SamplingCode"
  colnames(plateBGACountData)[3] <- "SamplingSite"
  colnames(plateBGACountData)[4] <- "SampleType"
  colnames(plateBGACountData)[12] <- "E.coli"
  colnames(plateBGACountData)[13] <- "Salmonella"

  bacteriaTypes <- c("E.coli", "Salmonella")

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

  subPlateBGAData <- averageCounts(waterData, waterReps, rawData, 12:13, bacteriaTypes)
  subPlateBGAData <- rbind(subPlateBGAData, averageCounts(sedimentData, sedimentReps,
                                                          rawData, 12:13, bacteriaTypes))

  subPlateBGAData <- subPlateBGAData[subPlateBGAData$Bacteria == "E.coli", ]

  # Add an extra column for the media type to each count table
  subPlateBGAData <- cbind(subPlateBGAData, rep("BGA", nrow(subPlateBGAData)))
  colnames(subPlateBGAData)[ncol(subPlateBGAData)] <- "Media"

  # qPCR
  colnames(qPCRData)[2] <- "Sampling Site"
  qPCRData <- addRivers(qPCRData, metaData)
  qPCRData <- cbind(qPCRData, rep("qPCR"))
  colnames(qPCRData)[ncol(qPCRData)] <- "Media"

  # Now combine the e. coli counts across the plates into one new table
  plateCountCombData <- rbind(subPlateCSData, subPlateSSAData, subPlateBGAData)
  plateCountCombData <- plateCountCombData[plateCountCombData$SampleType == "Sediment", ]

  # Convert the sample ids into the numeric part for sorting
  plateCountCombData$SamplingCode <-
    as.numeric(substr(plateCountCombData$SamplingCode, 5, 8))

  # Amend the column names
  colnames(plateCountCombData)[2] <- "SamplingSite"
  colnames(plateCountCombData)[3] <- "SampleType"
  colnames(plateCountCombData)[6] <- "CfuCount"

  mediaTypes <- c("ChromoSelect", "SSA", "BGA", "qPCR")

  colnames(qPCRData)[1] <- "SamplingCode"
  colnames(qPCRData)[2] <- "SamplingSite"
  colnames(qPCRData)[3] <- "SampleType"
  colnames(qPCRData)[7] <- "MeanCfu"

  qPCRData$SamplingCode <- as.numeric(substr(qPCRData$SamplingCode, 5,
                                             length(qPCRData$SamplingCode)))
  qPCRSubData <- qPCRData[, c(1, 2, 3, 5, 7:10)]

  subPlateComb <- plateCountCombData[, c(1, 2, 3, 5, 8:11)]
  colnames(qPCRSubData) <- colnames(subPlateComb)

  subPlateComb <- rbind(subPlateComb, qPCRSubData)

  subPlateComb$MeanCfu <- as.numeric(subPlateComb$MeanCfu)
  colnames(subPlateComb)[1] <- "SampleID"

  subPlateComb <- subPlateComb[order(as.numeric(substr(subPlateComb$Household, 11, length(subPlateData$Household))),
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

  sampleTypes <- c("Sediment", "Water")

  for (iType in 1:length(sampleTypes))
  {
    bactSubData <- subPlateComb[subPlateComb$SampleType ==
                                  sampleTypes[iType], ]

    for (iRow in (1:(nrow(bactSubData) / 3)))
    {
      meanVal <- mean(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                            bactSubData$SamplingSite[iRow * 3]) &
                                            (bactSubData$Media ==
                                               bactSubData$Media[iRow * 3])],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                        bactSubData$SamplingSite[iRow * 3]) &
                                        (bactSubData$Media ==
                                           bactSubData$Media[iRow * 3])],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 3, 1:4], meanVal, sdVal,
                          bactSubData[iRow * 3, 6:8])

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
    sumData$Media <- factor(sumData$Media, levels = mediaTypes)

    if (SampleType == 1)
    {
      sType = "W"
    }else
    {
      sType = "V"
    }

    BarPlotGastroPak(sumData, c("DS", "MS", "US"), "Species", "Env",
                     "E coli counts", "",
                     c('darkblue', 'pink',  'lightgreen', 'gold', "lightgrey"),
                     sType, 5, plotsDir, paste("Env Comp EC ", sampleTypes[iType]))
  }

}

