# ReadAndBarPlotEnvBGACounts
#
# Read in plate count data and barchart plot it
# One bar char for the Salmonella counts and one for E. coli
# and one for the coliforms
# Group counts by river and order by site
#
# file ReadAndBarPlotEnvBGACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  05/01/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotEnvBGACounts <- function(plateCountEnvBGAFile, plotsDir, metaDataFile,
                                       rawData, sedimentReps, waterReps) {

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

  plateBGACountData <- nameCols(plateBGACountData)

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)
  bacteriaTypes <- c("Salmonella", "E.coli")
  convertCols <- which(colnames(plateBGACountData) %in% bacteriaTypes)

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

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

  subPlateData <- subPlateData[order(subPlateData$SamplingSite), ]

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subPlateData[subPlateData$Bacteria ==
                                  bacteriaTypes[iBacteria], ]

    for (iRow in (1:(nrow(bactSubData) / 6)))
    {
      meanVal <- mean(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                             bactSubData$SamplingSite[iRow * 6]) &
                                            (bactSubData$SampleType == "Sediment")],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                         bactSubData$SamplingSite[iRow * 6]) &
                                        (bactSubData$SampleType == "Sediment")],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 6, c(1,3,6)], "Sediment", meanVal,
                          sdVal, bactSubData[iRow * 6, 10:11])

      colnames(sumDataRow)[4] <- "SampleType"

      if (iRow == 1)
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

      meanVal <- mean(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                             bactSubData$SamplingSite[iRow * 6]) &
                                            (bactSubData$SampleType == "Water")],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[(bactSubData$SamplingSite ==
                                         bactSubData$SamplingSite[iRow * 6]) &
                                        (bactSubData$SampleType == "Water")],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 6, c(1,3,6)], "Water", meanVal,
                          sdVal, bactSubData[iRow * 6, 10:11])

      colnames(sumDataRow)[4] <- "SampleType"
      sumData <- rbind(sumData, sumDataRow)
    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    sumData <- sumData[order(sumData$SamplingSite), ]
    colnames(sumData)[5:8] <- c("MeanCfu", "StdDev", "River", "Location")

    bacteriaPlot <- ggplot(data = sumData, aes(x = SamplingSite, y = MeanCfu,
                                               fill = SampleType))
    bacteriaPlot <- bacteriaPlot + geom_bar(stat = "identity",
                                            position = "dodge")
    bacteriaPlot <- bacteriaPlot +
      geom_errorbar(aes(ymin = MeanCfu - StdDev, ymax = MeanCfu + StdDev),
                    width = 0.2, position = position_dodge(.9))

    bacteriaPlot <- bacteriaPlot + theme(axis.text.x = element_text(size = 8))

    bacteriaPlot <- bacteriaPlot +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                       size = 12))
    bacteriaPlot <- bacteriaPlot + labs(y = bquote('Mean cfu '~ml^-1),
                                        title = paste("Brilliant green plates ",
                                                      bacteriaTypes[iBacteria],
                                                      sep = ""),
                                        fill = "Sample Type")
    bacteriaPlot <- bacteriaPlot + facet_wrap( ~ River, ncol = 3,
                                               scales = "free_x")
    bacteriaPlot <- bacteriaPlot +
      theme(strip.background = element_rect(fill = "white", colour = "grey"),
            strip.text = element_text(size = 8, colour = "black"))

    bacteriaPlot <- bacteriaPlot +
      theme(panel.background = element_rect(fill = "white",
                                            colour = "grey",
                                            size = 0.5, linetype = "solid"))
    bacteriaPlot <- bacteriaPlot +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(size = .5, color = "grey"),
            panel.grid.minor.y = element_line(size = .2, color = "grey"))

    bacteriaPlot <- bacteriaPlot + theme(axis.title.y = element_text(size = 12))
    bacteriaPlot <- bacteriaPlot + theme(axis.title.x = element_blank())
    bacteriaPlot <- bacteriaPlot + theme(plot.title = element_text(size = 14,
                                                                   hjust = 0.5))
    bacteriaPlot <- bacteriaPlot + theme(legend.title = element_text(size = 10),
                                         legend.text = element_text(size = 8))

    bacteriaPlot <- bacteriaPlot + theme(axis.text.y = element_text(size = 8))

    bacteriaPlot <- bacteriaPlot +
      scale_fill_manual(values = c('chocolate4', 'skyblue'))
    bacteriaPlot <- bacteriaPlot + scale_y_continuous(trans = 'log10')
    bacteriaPlot <- bacteriaPlot + scale_x_discrete(labels = c("Upstream",
                                                               "Midstream",
                                                               "Downstream"))

    bacteriaPlot
    ggsave(paste(plotsDir, bacteriaTypes[iBacteria], ".pdf", sep = ""),
           bacteriaPlot, width = 5, height = 6)

  }

}
