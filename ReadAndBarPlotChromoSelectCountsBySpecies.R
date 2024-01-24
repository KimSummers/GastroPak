# ReadAndBarPlotChromoSelectEnvCountsBySpecies
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotChromoSelectEnvCountsBySpecies
#
# inputs
# 	plateCountCSFile  - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/07/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotChromoSelectEnvCountsBySpecies <- function(plateCountEnvCSFile,
                                                         plotsDir, metaDataFile,
                                                         rawData, sedimentReps,
                                                         waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCSCountData <- read_csv(plateCountEnvCSFile)
  metaData <- read.csv(metaDataFile)

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

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "Env")
  subPlateData <- rbind(subPlateData, averageCounts(sedimentData, sedimentReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes, "Env"))

  sampleNumbers <- unique(subPlateData$SampleID)
  subPlateData$SampleID <- factor(subPlateData$SampleID,
                                     levels = sampleNumbers)

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

  subPlateData <- subPlateData[order(subPlateData$SampleID, subPlateData$Bacteria), ]

  for (iSampleType in 1:length(sampleType))
  {
    sampleSubData <- subPlateData[subPlateData$SampleType ==
                                    sampleType[iSampleType], ]
    groupLength <- 3 * length(bacteriaTypes)

    for (iRow in (1:(nrow(sampleSubData) / groupLength)))
    {

      for (iBacteria in 1:length(bacteriaTypes))
      {
        meanVal <- mean(sampleSubData$MeanCfu[(sampleSubData$SamplingSite ==
                                                 sampleSubData$SamplingSite[iRow * groupLength]) &
                                                (sampleSubData$Bacteria == bacteriaTypes[iBacteria])],
                        na.rm = TRUE)
        sdVal <- sd(sampleSubData$MeanCfu[(sampleSubData$SamplingSite ==
                                             sampleSubData$SamplingSite[iRow * groupLength]) &
                                            (sampleSubData$Bacteria == bacteriaTypes[iBacteria])],
                    na.rm = TRUE)
        sumDataRow <- cbind(sampleSubData[iRow * groupLength, c(1, 3, 4)],
                            bacteriaTypes[iBacteria], meanVal, sdVal,
                            sampleSubData[iRow * groupLength, 10:11])

        colnames(sumDataRow)[4] <- "Bacteria"

        if ((iRow == 1) & (iBacteria == 1))
        {
          sumData <- sumDataRow
        }else
        {
          sumData <- rbind(sumData, sumDataRow)
        }

      }

    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    sumData$Bacteria <- factor(sumData$Bacteria, levels = unique(sumData$Bacteria))
    sumData <- sumData[order(sumData$SamplingSite), ]
    colnames(sumData)[5:8] <- c("MeanCfu", "StdDev", "River", "Location")

    if (SampleType[iSample] == "Water")
    {
      measureType = "V"
    }else
    {
      measureType = "W"
    }

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Species",
                     "Env", "ChromoSelect plates", sampleType[iSampleType],
                     c('darkblue','tomato2', 'pink'), measureType, 3, plotsDir,
                     sampleType[iSampleType])
  }

}



