# ReadAndBarPlotHHBGACounts
#
# Read in plate count data and barchart plot it
# One bar chart for the E. coli counts and one for Salmonella
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHBGACounts
#
# inputs
# 	plateCountFile  - file containing plate count data
#   plotsDir        - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   sedimentReps      - Number of sediment dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  16/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHSSACounts <- function(plateCountFile, plotsDir, metaDataFile,
                                      rawData, sedimentReps, waterReps) {

  library(tidyverse)
  library(scales)
  library(ggplot2)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountHHBGAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  plateBGACountData <- numericCounts(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$`Sample Type` == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$`Sample Type` == "Sediment", ]

  waterData <- fixMissingData(waterData, waterReps)
  sedimentData <- fixMissingData(sedimentData, sedimentReps)

  plateBGACountData <- rbind(waterData, sedimentData)

  plateBGACountData <- coloursToSpeciesCS(plateBGACountData)
  plateBGACountData <- addHousehold(plateBGACountData, metaData)

  plateBGACountData$`Sample-ID` <-
    as.numeric(substr(plateBGACountData$`Sample-ID`, 5, 8))

  plateBGACountData <- nameCols(plateBGACountData)
  colnames(plateBGACountData)[1] <- "SamplingCode"
  colnames(plateBGACountData)[3] <- "SamplingSite"
  colnames(plateBGACountData)[4] <- "SampleType"
  colnames(plateBGACountData)[9] <- "E.coli"
  colnames(plateBGACountData)[12] <- "Salmonella"
  colnames(plateBGACountData)[13] <- "Household"

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Salmonella", "E.coli")

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  sedimentData <- plateBGACountData[plateBGACountData$SampleType == "Sediment", ]

  subPlateData <- averageCounts(waterData)
  subPlateData <- rbind(subPlateData, sedimentData)

  colnames(subPlateData)[2] <- "SampleID"

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11,
                                                       length(subPlateData$Household))),
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
    #                                      title = paste("Brilliant green plates ",
    #                                                    bacteriaTypes[iBacteria],
    #                                                    sep = ""),
    #                                      fill = "Sample Type")
    bacteriaPlot <- bacteriaPlot + labs(y = "Mean cfu / ml (water) or cfu / 0.1 g (faeces)",
                                        title = paste("Brilliant green plates ",
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

