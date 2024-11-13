# ReadAndPlotGPHuBac
#
# Read in metagenomic data and bar chart plot it
#
# file ReadAndPlotGPHuBac
#
# inputs
# 	metagenomicDataFile   - file containing plate count data
# 	abClassesFile         - file containing antibiotic class data
# 	normalisationFile     - file containing normalisation data
#   huBacCovFile          - file containing HuBac and RuBac coverage data
# 	metaDataFile          - file containing metadata
#   campaignType          - "Env" or "HH" for environmental or household
#   normType              - type of normalisation to use
#   multiplier            - value to multiple data by
#   average               - TRUE = average data across samples, FALSE do not
#   absLab                - Label for the y-axis for absolute data plots
#   yAxisLims             - Limits for the y-axis
#   plotsDir              - directory to store plots in
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  12/08/24  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndPlotGPHuBacOnly <- function(normalisationFile, huBacCovFile, metaDataFile, 
                                   campaignType, normType, multiplier, average,
                                   sites, absLab, yAxisLims, plotsDir) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  normData <- read_tsv(normalisationFile)
  metaData <- read_csv(metaDataFile)
  huBacData <- read_tsv(huBacCovFile)

  colnames(huBacData)[1] <- "Bac"
  
  if (campaignType == "Env")
  {
    huBacData <- huBacData[, c(1, grep("PAK0", colnames(huBacData)))]
    colnames(huBacData) <- sub("coassemblies_Pak_sediments_PAK0", "PAK0", 
                               colnames(huBacData))
  }else
  {
    huBacData <- huBacData[, grep("PHS0", colnames(huBacData))]
    colnames(huBacData) <- sub("coassemblies_Pak_faecal_PHS0", "PHS0", 
                               colnames(huBacData))
  }
  
  normHuBac <- huBacData
  
  for (iSample in 2:(ncol(normHuBac) - 1))
  {
    
    for (iData in 1:nrow(normHuBac))
    {
      normHuBac[iData, iSample] <- normHuBac[iData, iSample] /
        normData[which(normData$Normalisation == normType),
                 which(colnames(normData) == colnames(normHuBac)[iSample])]
    }
    
  }
  
  
  # Transpose data for ease of working - this then requires further adjustments
  normHuBac <- as.data.frame(t(normHuBac))

  colnames(normHuBac) <- normHuBac[1, ]
  normHuBac <- normHuBac[2:nrow(normHuBac), ]
  normHuBac <- cbind(rownames(normHuBac), normHuBac)

  # Ensure sample naming is consistent across data sources
  colnames(normHuBac)[1] <- "SampleID"
  rownames(normHuBac) <- 1:nrow(normHuBac)
  normHuBac$SampleID <- sub("PAK0", "PAK_0", normHuBac$SampleID)
  normHuBac$SampleID <- sub("PHS0", "PHS_0", normHuBac$SampleID)
  
  if (campaignType == "Env")
  {
    # Add sampling site and season
    normHuBac <- addSamplingSite(normHuBac, metaData)
    normHuBac <- normHuBac %>% relocate(SamplingSite, .after = SampleID)
    normHuBac <- normHuBac %>% relocate(Season, .after = SamplingSite)

    startCol <- 4
  }else
  {
    startCol <- 2
  }

  for (iCol in startCol:ncol(normHuBac))
  {
    normHuBac[, iCol] <- as.numeric(normHuBac[, iCol])
  }
  
  if (campaignType == "Env")
  {
    # Split data into seasons
    normHuBacCamp <- normHuBac[normHuBac$Season == "Summer Dry", ]
    normHuBacCamp <- normHuBacCamp[normHuBacCamp$SamplingSite %in% sites, ]
  }else
  {
    normHuBacCamp <- normHuBac
  }

  # aggregate  samples then pivot data
  if (average)
  {
    normBacAverage <- aggregate(normHuBacCamp[ , startCol:ncol(normHuBacCamp)],
                                by = list(normHuBacCamp$SamplingSite),
                                FUN = mean)
    colnames(normBacAverage)[1] <- "SamplingSite"
    metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(2:ncol(normBacAverage)),
                                        names_to = "Bac", values_to = "Abundance")
    xAxisLabels <- normBacAverage$SamplingSite
  }else
  {
    metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(startCol:ncol(normBacAverage)),
                                        names_to = "Bac", values_to = "Abundance")

    if (campaignType == "Env")
    {
      xAxisLabels <- normHuBacCamp$SamplingSite
    }else
    {
      xAxisLabels <- normHuBacCamp$SampleID
    }

  }

  metagenomicsBacLong <- cbind(metagenomicsBacLong[, 1], metagenomicsBacLong)
  colnames(metagenomicsBacLong)[1] <- "xVals"
  
  if (campaignType == "Env")
  {
    plotTitleBase <- "Dry Summer 2022"
  }else
  {
    plotTitleBase <- "Household"
  }

  # Plot dry season absolute data
  BarPlotMetagenomics(graphData = metagenomicsBacLong, xAxisLabels = xAxisLabels, 
                      absolute = TRUE, absLab = absLab, graphType =  "HuBac", 
                      yAxisVals = yAxisLims, plotsDir = plotsDir,
                      plotTitle = paste(plotTitleBase, "HuBac", average, normType,
                                        "Normalised Absolute Abundances"))

  if (campaignType == "Env")
  {
    normHuBacWet <- normHuBac[normHuBac$Season == "Summer Wet", ]
    normHuBacWet <- normHuBacWet[normHuBacWet$SamplingSite %in% sites, ]
    
    # aggregate  samples then pivot data
    if (average)
    {
      normBacAverage <- aggregate(normHuBacWet[ , 4:ncol(normHuBacWet)],
                                  by = list(normHuBacWet$SamplingSite),
                                  FUN = mean)
      colnames(normBacAverage)[1] <- "SamplingSite"
      
      metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(2:ncol(normBacAverage)),
                                          names_to = "Bac", values_to = "Abundance")
      xAxisLabels <- normBacAverage$SamplingSite
    }else
    {
      normBacAverage <- normHuBacWet
      metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(4:ncol(normBacAverage)),
                                          names_to = "Bac", values_to = "Abundance")
      if (campaignType == "Env")
      {
        xAxisLabels <- normHuBacCamp$SamplingSite
      }else
      {
        xAxisLabels <- normHuBacCamp$SampleID
      }
      
    }

    metagenomicsBacLong <- cbind(metagenomicsBacLong[, 1], metagenomicsBacLong)
    colnames(metagenomicsBacLong)[1] <- "xVals"
    
    # Plot wet season absolute data
    BarPlotMetagenomics(graphData = metagenomicsBacLong, xAxisLabels = xAxisLabels,
                        absolute =  TRUE, absLab = absLab, graphType = "HuBac",
                        yAxisVals = yAxisLims, plotsDir = plotsDir,
                        plotTitle = paste("Wet 2022 HuBac", average, normType,
                                          "Normalised Absolute Abundances"))
  }

}
