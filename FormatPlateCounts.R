# FormatPlateCounts
#
# Read in plate count data and format it to mean cfu / ml accross dilutions
# Also totals counts for species depending on media type and only keeps relevant data
#
# file FormatPlateCounts
#
# inputs
# 	plateCountData    - Plate count data
# 	mediaType         - "ChromoSelect", "BGA" or "SSA"
#   formattedDir      - directory to store formatted files in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  02/08/23  Wellington Lab - School of Life Sciences - University of Warwick

FormatPlateCounts <- function(plateCountData, mediaType, plotsDir, formattedDir,
                              rawData, sedimentReps, waterReps)
{
  library(tidyverse)

  plateCountData <- read_csv(plateCountFile)

  plateCountData <- numericCounts(plateCountData)

  waterData <- plateCountData[plateCountData$`Sample Type` == "Water", ]
  sedimentData <- plateCountData[plateCountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData)
  sedimentData <- fixMissingData(sedimentData)
  plateCountData <- rbind(waterData, sedimentData)

  if (mediaType == "ChromoSelect")
  {
    plateCountData <- coloursToSpeciesCS(plateCountData)
    plateCountData <- nameCols(plateCountData)
    convertCols <- c(which(colnames(plateCountData) == "E.coli"),
                     which(colnames(plateCountData) == "Coliforms"),
                     which(colnames(plateCountData) == "Klebsiella"))

    plateCountData <- plateCountData %>% relocate(E.coli,
                                                      .after = Dilution)
    plateCountData <- plateCountData %>% relocate(Coliforms,
                                                      .after = E.coli)
    plateCountData <- plateCountData %>% relocate(Klebsiella,
                                                      .after = Coliforms)
  }else
  {

    if (mediaType == "SSA")
    {
      plateCountData <- coloursToSpeciesSSA(plateCountData)
      plateCountData <- nameCols(plateCountData)

      convertCols <- c(which(colnames(plateCountData) == "Salmonella"),
                       which(colnames(plateCountData) == "E.coli"),
                       which(colnames(plateCountData) == "Shigella"))
      plateCountData <- plateCountData %>% relocate(Salmonella,
                                                        .after = Dilution)
      plateCountData <- plateCountData %>% relocate(Shigella,
                                                        .after = Salmonella)
      plateCountData <- plateCountData %>% relocate(E.coli,
                                                        .after = Shigella)
    }else
    {

      if (mediaType == "BGA")
      {
        plateCountData <- coloursToSpeciesBGA(plateCountData)
        plateCountData <- nameCols(plateCountData)

        convertCols <- c(which(colnames(plateCountData) == "Salmonella"),
                         which(colnames(plateCountData) == "E.coli"))
        plateCountData <- plateCountData %>% relocate(Salmonella,
                                                          .after = Dilution)
        plateCountData <- plateCountData %>% relocate(E.coli,
                                                          .after = Salmonella)

      }

    }

  }

  waterData <- plateCountData[plateCountData$SampleType == "Water", ]
  sedimentData <- plateCountData[plateCountData$SampleType == "Sediment", ]

  if (!rawData)
  {
    plateCountData <- convertCfuCounts(waterData, waterReps, convertCols)
    plateCountData <- rbind(plateCountData, convertCfuCounts(sedimentData, sedimentReps,
                                                             convertCols))
  }

  plateCountData$`Sample-ID` <-
    as.numeric(substr(plateCountData$`Sample-ID`, 5, 8))
}
