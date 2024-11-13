# FormatPlateCounts
#
# Read in plate count data and format it to mean cfu / ml across dilutions
# Also totals counts for species depending on media type and only keeps relevant data
#
# file FormatPlateCounts
#
# inputs
# 	plateCountFile      - Plate count data file
# 	mediaType           - "ChromoSelect", "BGA" or "SSA"
#   formattedFile       - file name for formatted file
#   metaDataFile        - file containing river or household meta data
#   rawData             - True for absolute counts, false for cfu / ml
#   otherReps           - Number of sediment or faecal dilutions
#   waterReps           - Number of water dilutions
#   campaignType        - "Env" for environmental and "HH" for household

# Version    Author       Date      Affiliation
# 1.00       J K Summers  02/08/23  Wellington Lab - School of Life Sciences - University of Warwick

FormatPlateCounts <- function(plateCountFile, mediaType, formattedFile,
                              metaDataFile, rawData, otherReps, waterReps,
                              campaignType)
{
  library(tidyverse)

  plateCountData <- read_csv(plateCountFile)
  metaData <- read_csv(metaDataFile)

  # Ensure count data is numeric, collate by species and name columns correctly
  plateCountData <- numericCounts(plateCountData)

  if (mediaType == "ChromoSelect")
  {
    plateCountData <- coloursToSpeciesCS(plateCountData)
  }else
  {

    if (mediaType == "SSA")
    {
      plateCountData <- coloursToSpeciesSSA(plateCountData)
    }else
    {
      plateCountData <- coloursToSpeciesBGA(plateCountData)
    }

  }

  plateCountData <- nameCols(plateCountData)

  # Fix gaps in data
  waterData <- plateCountData[plateCountData$SampleType == "Water", ]

  if (campaignType == "Env")
  {
    otherData <- plateCountData[plateCountData$SampleType == "Sediment", ]
  }else
  {
    otherData <- plateCountData[plateCountData$SampleType == "Stool", ]
  }

  waterData <- fixMissingData(waterData, waterReps)
  otherData <- fixMissingData(otherData, otherReps)

  plateCountData <- rbind(waterData, otherData)

  if (campaignType == "Env")
  {
    # Add meta data on rivers to table
    plateCountData <- addRivers(plateCountData, metaData)
  }else
  {
    # Add meta data on households to table
    plateCountData <- addHouseholds(plateCountData, metaData)
  }

  if (mediaType == "ChromoSelect")
  {
    plateCountData <- plateCountData %>% relocate(Klebsiella, .after = coliforms)

    bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")
    convertCols <- which(colnames(plateCountData) %in% bacteriaTypes)

  }else
  {

    if (mediaType == "SSA")
    {
      plateCountData <- plateCountData %>% relocate(Salmonella, .after = Shigella)
      plateCountData <- plateCountData %>% relocate(E.coli, .after = Salmonella)

      bacteriaTypes <- c("Shigella", "Salmonella", "E.coli")

      convertCols <- which(colnames(plateCountData) %in% bacteriaTypes)
    }else
    {

      if (mediaType == "BGA")
      {
        plateCountData <- plateCountData %>% relocate(E.coli, .after = Salmonella)
        bacteriaTypes <- c("Salmonella", "E.coli")
        convertCols <- which(colnames(plateCountData) %in% bacteriaTypes)
      }

    }

  }

  waterData <- plateCountData[plateCountData$SampleType == "Water", ]

  if (campaignType == "Env")
  {
    otherData <- plateCountData[plateCountData$SampleType == "Sediment", ]
  }else
  {
    otherData <- plateCountData[plateCountData$SampleType == "Stool", ]
  }

  # Average plate counts per species and per sample across all dilutions
  subPlateData <- averageCounts(countData = waterData, numDilutions = waterReps,
                                rawData = rawData, convertCols = convertCols,
                                bacteriaTypes =  bacteriaTypes, campaignType = campaignType)
  subPlateData <- rbind(subPlateData, averageCounts(countData = otherData,
                                                    numDilutions = otherReps,
                                                    rawData = rawData,
                                                    convertCols = convertCols,
                                                    bacteriaTypes = bacteriaTypes,
                                                    campaignType = campaignType))

  write_csv(subPlateData, formattedFile)
}
