# ReadAndBarPlotChromocultCounts
#
# Read in plate count data and barchart plot it
# One bar char for the E. coli counts, one for Klebsiella, one for Salmonella / Shigella
# and one for the coliforms
# Group counts by household and order by household, sample type and sample id
#
# file ReadAndBarPlotChromocultCounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile    - file containing river meta data
#   rawData         - True for absolute counts, false for cfu / ml
#   stoolReps       - Number of stool dilutions
#   waterReps       - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  15/12/22  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotChromocultCounts <- function(plateCountHHCCFile, plotsDir, metaDataFile,
                                           rawData, stoolReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateCCCountData <- read_csv(plateCountHHCCFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  plateCCCountData <- numericCounts(plateCCCountData)

  waterData <- plateCCCountData[plateCCCountData$`Sample Type` == "Water", ]
  stoolData <- plateCCCountData[plateCCCountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateCCCountData <- rbind(waterData, stoolData)

  plateCCCountData <- coloursToSpeciesCS(plateCCCountData)
  plateCCCountData <- addHousehold(plateCCCountData, metaData)

  plateCCCountData$`Sample-ID` <-
    as.numeric(substr(plateCCCountData$`Sample-ID`, 5, 8))

  plateCCCountData <- nameCols(plateCCCountData)

  plateCCCountData <- plateCCCountData %>% relocate(Klebsiella,
                                                    .after = coliformCounts)

  bacteriaTypes <- c("E.coli", "coliforms", "Klebsiella")

  waterData <- plateCCCountData[plateCCCountData$SampleType == "Water", ]
  stoolData <- plateCCCountData[plateCCCountData$SampleType == "Stool", ]

  subPlateData <- averageCounts(waterData, waterReps, rawData, 9:11, bacteriaTypes)
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps, rawData,
                                                    9:11, bacteriaTypes))

  colnames(subPlateData)[2] <- "SampleID"

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11, length(subPlateData$Household))),
                                     subPlateData$SampleType,
                                     subPlateData$SampleID), ]
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

  households <- unique(subPlateData$Household)
  subPlateData$Household <- factor(subPlateData$Household,
                                   levels = households)

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subPlateData[subPlateData$Bacteria ==
                                  bacteriaTypes[iBacteria], ]

    for (iRow in (1:(nrow(bactSubData) / 3)))
    {
      meanVal <- mean(bactSubData$MeanCfu[bactSubData$SampleID ==
                                            bactSubData$SampleID[iRow * 3]],
                      na.rm = TRUE)
      sdVal <- sd(bactSubData$MeanCfu[bactSubData$SampleID ==
                                        bactSubData$SampleID[iRow * 3]],
                  na.rm = TRUE)
      sumDataRow <- cbind(bactSubData[iRow * 3, 1:8], meanVal, sdVal,
                          bactSubData[iRow * 3, 10])

      if (iRow == 1)
      {
        sumData <- sumDataRow
      }else
      {
        sumData <- rbind(sumData, sumDataRow)
      }

    }

    sumData$meanVal[sumData$meanVal == 0] <- NaN
    colnames(sumData)[9:11] <- c("MeanCfu", "StdDev", "Household")

    bacteriaPlot <- ggplot(data = sumData, aes(x = SampleID, y = MeanCfu,
                                               fill = SampleType))
    bacteriaPlot <- bacteriaPlot + geom_bar(stat = "identity", position = "dodge")
    bacteriaPlot <- bacteriaPlot +
      geom_errorbar(aes(ymin = MeanCfu - sdVal, ymax = MeanCfu + sdVal), width = 0.2)

    bacteriaPlot <- bacteriaPlot + theme(axis.text.x = element_text(size = 8))
  #  bacteriaPlot <- bacteriaPlot + labs(y = bquote('Mean cfu '~ml^-1),
  #                                      title = paste("Chromocult plates ",
  #                                                    bacteriaTypes[iBacteria],
  #                                                    sep = ""),
  #                                      fill = "Sample Type")
    bacteriaPlot <- bacteriaPlot + labs(y = "Mean cfu / ml (water) or cfu / 0.1 g (faeces)",
                                        title = paste("Chromocult plates ",
                                                      bacteriaTypes[iBacteria],
                                                      sep = ""),
                                        fill = "Sample Type")
    bacteriaPlot <- bacteriaPlot + facet_wrap( ~ Household, ncol = 4,
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

    bacteriaPlot
    ggsave(paste(plotsDir, bacteriaTypes[iBacteria], ".pdf", sep = ""),
           bacteriaPlot, width = 5, height = 6)

  }

}

