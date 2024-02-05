# ReadAndBarPlotEnvSSACountsBySpecies
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Samonella, and one for Shigella
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvSSACountsBySpecies
#
# inputs
# 	plateCountSSAFile- file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions
#   season            - The season the data is from

# Version    Author       Date      Affiliation
# 1.00       J K Summers  16/07/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvSSACountsBySpecies <- function(plateCountSSAFile, plotsDir,
                                                metaDataFile, rawData,
                                                sedimentReps, waterReps,
                                                season) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateSSACountData <- read_csv(plateCountSSAFile)
  metaData <- read.csv(metaDataFile)

  plateSSACountData <- numericCounts(plateSSACountData)

  waterData <- plateSSACountData[plateSSACountData$`Sample Type` == "Water", ]
  sedimentData <- plateSSACountData[plateSSACountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateSSACountData <- rbind(waterData, sedimentData)

  plateSSACountData <- coloursToSpeciesSSA(plateSSACountData)
  plateSSACountData <- addRivers(plateSSACountData, metaData)

  plateSSACountData <- nameCols(plateSSACountData)

  plateSSACountData <- plateSSACountData %>% relocate(Salmonella,
                                                      .after = Shigella)
  plateSSACountData <- plateSSACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Shigella", "Salmonella", "E.coli")

  convertCols <- which(colnames(plateSSACountData) %in% bacteriaTypes)

  waterData <- plateSSACountData[plateSSACountData$SampleType == "Water", ]
  sedimentData <- plateSSACountData[plateSSACountData$SampleType == "Sediment", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "Env")
  subPlateData <- rbind(subPlateData, averageCounts(sedimentData, sedimentReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes, "Env"))

  sampleCodes <- unique(subPlateData$SamplingCode)
  subPlateData$SamplingCode <- factor(subPlateData$SamplingCode,
                                      levels = sampleCodes)

  sampleSites <- unique(subPlateData$SamplingSite)
  subPlateData$SamplingSite <- factor(subPlateData$SamplingSite,
                                      levels = sampleSites)

  sampleType <- unique(subPlateData$SampleType)
  subPlateData$SampleType <- factor(subPlateData$SampleType,
                                    levels = sampleType)

  subPlateData$Bacteria <- factor(subPlateData$Bacteria,
                                  levels = bacteriaTypes)

  subPlateData <- subPlateData[order(subPlateData$SamplingSite, subPlateData$Bacteria), ]

  for (iSampleType in 1:length(sampleType))
  {
    sumData <- subSampleEnvData(subPlateData, sampleType[iSampleType], bacteriaTypes)

    if (sampleType[iSampleType] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Species",
                     "Env", "Salmonella-Shigella agar plates", sampleType[iSampleType],
                     c('grey', 'black', 'pink'), measureType, 3, plotsDir,
                     paste("SSA Env", season, sampleType[iSampleType]))
  }

}
