# ReadAndBarPlotBGACountsBySpecies
#
# Read in plate count data and barchart plot it
# One bar char for the Salmonella counts and one for E. coli
# Group counts by river and order by site
#
# file ReadAndBarPlotBGACountsBySpecies
#
# inputs
# 	plateCountBGAFile - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  17/07/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotBGACounts <- function(plateCountBGAFile, plotsDir,
                                    metaDataFile, rawData,
                                    sedimentReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountEnvBGAFile)
  plateBGACountData <- plateBGACountData[1:279, ]
  metaData <- read.csv(metaDataFile)

  plateBGACountData <- numericCounts(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$`Sample Type` == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateBGACountData <- rbind(waterData, sedimentData)

  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- addRivers(plateBGACountData, metaData)

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)
  convertCols <- c(which(colnames(plateCountData) == "Salmonella"),
                   which(colnames(plateCountData) == "E.coli"))

  bacteriaTypes <- c("Salmonella", "E.coli")

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes)
  subPlateData <- rbind(subPlateData, averageCounts(sedimentData, sedimentReps,
                                                    rawData, convertCols,
                                                    bacteriaTypes))

  colnames(subPlateData)[2] <- "SampleID"

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

    BarPlotGastroPak(sumData, c("Upstream", "Midstream", "Downstream"), "Species",
                     "Env", "Brilliant green agar", sampleType[iSampleType],
                     c('palevioletred', 'yellowgreen'), 3, plotDir,
                     sampleType[iSampleType])
  }

}

