# ReadAndBarPlotHHBGACounts
#
# Read in plate count data and barchart plot it
# One bar chart for the E. coli counts and one for Salmonella
# Group counts by house and order by sample type and number
#
# file ReadAndBarPlotHHBGACounts
#
# inputs
# 	plateCountFile    - file containing plate count data
#   plotsDir          - directory to store plots in
#   metaDataFile      - file containing river meta data
#   rawData           - True for absolute counts, false for cfu / ml
#   stoolReps         - Number of stool dilutions
#   waterReps         - Number of water dilutions

# Version    Author       Date      Affiliation
# 1.00       J K Summers  16/06/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotHHBGACounts <- function(plateCountFile, plotsDir, metaDataFile,
                                      rawData, stoolReps, waterReps) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  plateBGACountData <- read_csv(plateCountHHBGAFile)
  metaData <- read_csv(metaDataFile)
  metaData <- metaData[2:nrow(metaData), ]

  plateBGACountData <- numericCounts(plateBGACountData)

  waterData <- plateBGACountData[plateBGACountData$`Sample Type` == "Water", ]
  stoolData <- plateBGACountData[plateBGACountData$`Sample Type` == "Stool", ]

  waterData <- fixMissingData(waterData, waterReps)
  stoolData <- fixMissingData(stoolData, stoolReps)

  plateBGACountData <- rbind(waterData, stoolData)

  plateBGACountData <- coloursToSpeciesBGA(plateBGACountData)
  plateBGACountData <- addHouseholds(plateBGACountData, metaData)

  plateBGACountData <- nameCols(plateBGACountData)

  plateBGACountData$SampleID <-
    as.numeric(substr(plateBGACountData$SampleID, 5, 8))

  plateBGACountData <- plateBGACountData %>% relocate(E.coli,
                                                      .after = Salmonella)

  bacteriaTypes <- c("Salmonella", "E.coli")

  waterData <- plateBGACountData[plateBGACountData$SampleType == "Water", ]
  stoolData <- plateBGACountData[plateBGACountData$SampleType == "Stool", ]

  convertCols <- c(which(colnames(plateBGACountData) == "Salmonella"),
                   which(colnames(plateBGACountData) == "E.coli"))

  subPlateData <- averageCounts(waterData, waterReps, rawData, convertCols,
                                bacteriaTypes, "HH")
  subPlateData <- rbind(subPlateData, averageCounts(stoolData, stoolReps, rawData,
                                                    convertCols, bacteriaTypes, "HH"))

  colnames(subPlateData)[which(colnames(subPlateData) == "Sample-ID")] <- "SampleID"

  subPlateData <- subPlateData[order(as.numeric(substr(subPlateData$Household, 11,
                                                       length(subPlateData$Household))),
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

  households <- unique(subPlateData$Household)
  subPlateData$Household <- factor(subPlateData$Household,
                                   levels = households)

  for (iBacteria in 1:length(bacteriaTypes))
  {
    bactSubData <- subBacteriaHHData(subPlateData, bacteriaTypes[iBacteria],
                                     c("Water", "Stool"))

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

