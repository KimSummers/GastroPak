# ReadAndBarPlotMetagenomics
#
# Read in metagenomic data and bar chart plot it
# One bar chart for the ResFinder combined abundance data, one for India only
# one for combined percentage data, one for India only and one for Hubac data
#
# file ReadAndBarPlotMetagenomics
#
# inputs
# 	metagenomicDataFile   - file containing plate count data
#   plotsDir              - directory to store plots in
# Version    Author       Date      Affiliation
# 1.00       J K Summers  21/11/23  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndBarPlotMetagenomics <- function(metagenomicDataFile, plotsDir) {

  library(tidyverse)
  library(scales)
  library(ggplot2)

  # Read in data and put it in a dataframe
  metagenomicsData <- read_csv(metagenomicDataFilePol)

  metagenomicsData <- metagenomicsData[c(13:15, 10:12, 6:9, 1:5, 26:27, 25, 24,
                                         22:23, 20:21, 18:19, 16:17, 38:39, 37,
                                         36, 34:35, 32:33, 30:31, 28:29), ]
  metagenomicsData$Season_WD[metagenomicsData$Season_WD == "Dry_Hot"] <- "Dry"
  metagenomicsData$Season_WD[metagenomicsData$Season_WD == "Wet_Hot"] <- "Wet"

  xLabsIndia <- metagenomicsData$site_ID[substr(metagenomicsData$Sample_description, 1, 1) == "I"]
  xLabsUKSed <- metagenomicsData$Season_WD[(substr(metagenomicsData$Sample_description, 1, 1) == "S") &
                                             (metagenomicsData$Sample_type == "sediment")]
  xLabsUKWater <- metagenomicsData$Season_WD[(substr(metagenomicsData$Sample_description, 1, 1) == "S") &
                                             (metagenomicsData$Sample_type == "water")]
  xLabsComb <- c(xLabsIndia, xLabsUKSed, xLabsUKWater)

  xVals <- paste(metagenomicsData$Season_WD, metagenomicsData$site_ID)

  metagenomicsData <- cbind(metagenomicsData, as.data.frame(xVals), as.data.frame(xLabsComb))
  metagenomicsDataLong <- pivot_longer(metagenomicsData, cols = c(9:(ncol(metagenomicsData) - 4)),
                                       names_to = "ARG", values_to = "Abundance")
  metagenomicsDataLongSed <- metagenomicsDataLong[metagenomicsDataLong$Sample_type != "water", ]
  metagenomicsDataSed <- metagenomicsData[metagenomicsData$Sample_type != "water", ]

  BarPlotMetagenomics(metagenomicsDataLongSed,
                      xLabsComb, TRUE, FALSE, plotsDir, "Combined Abundances")

  metLongIndia <-
    metagenomicsDataLong[substr(metagenomicsDataLong$Sample_description, 1, 1) == "I", ]
  metIndia <- metagenomicsData[substr(metagenomicsDataLong$Sample_description, 1, 1) == "I", ]

  BarPlotMetagenomics(metLongIndia, xLabsIndia, TRUE, FALSE, plotsDir, "India Abundances")

  metLongUKSed <-
    metagenomicsDataLong[(substr(metagenomicsDataLong$Sample_description, 1, 1) == "S") &
                           (metagenomicsDataLong$Sample_type == "sediment"), ]
  metUKSed <-
    metagenomicsData[(substr(metagenomicsData$Sample_description, 1, 1) == "S") &
                       (metagenomicsData$Sample_type == "sediment"), ]

  BarPlotMetagenomics(metLongUKSed, xLabsUK, TRUE, FALSE, plotsDir, "UK Sediment Abundances")

  metLongUKWater <-
    metagenomicsDataLong[substr(metagenomicsDataLong$Sample_description, 1, 1) == "S" &
                           metagenomicsDataLong$Sample_type == "water", ]
  metUKWater <-
    metagenomicsData[substr(metagenomicsData$Sample_description, 1, 1) == "S" &
                       metagenomicsData$Sample_type == "water", ]

  BarPlotMetagenomics(metLongUKWater, xLabsUK, TRUE, FALSE, plotsDir, "UK Water Abundances")

  BarPlotMetagenomics(metagenomicsDataLongSed, xLabsComb, FALSE, FALSE, plotsDir, "Combined Percent Abundances")
  BarPlotMetagenomics(metLongIndia, xLabsIndia, FALSE, FALSE, plotsDir, "India Percent Abundances")
  BarPlotMetagenomics(metLongUKSed, xLabsUK, FALSE, FALSE, plotsDir, "UK Sediment Percent Abundances")
  BarPlotMetagenomics(metLongUKWater, xLabsUK, FALSE, FALSE, plotsDir, "UK Water Percent Abundances")

  BarPlotMetagenomics(metagenomicsDataSed, xLabsComb, FALSE, TRUE, plotsDir, "Combined Hubac")
  BarPlotMetagenomics(metIndia, xLabsIndia, FALSE, TRUE, plotsDir, "India Hubac")
  BarPlotMetagenomics(metUKSed, xLabsUK, FALSE, TRUE, plotsDir, "UK Sediment Hubac")
  BarPlotMetagenomics(metUKWater, xLabsUK, FALSE, TRUE, plotsDir, "UK Water Hubac")

}



