# ReadAndPlotGPMetagenomics
#
# Read in metagenomic data and bar chart plot it
#
# file ReadAndPlotGPMetagenomics
#
# inputs
# 	metagenomicDataFile   - file containing plate count data
# 	abClassesFile         - file containing antibiotic class data
# 	normalisationFile     - file containing normalisation data
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
# 1.00       J K Summers  26/03/24  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndPlotGPMetagenomics <- function(metagenomicDataFile, abClassesFile,
                                      normalisationFile, metaDataFile, campaignType,
                                      normType, multiplier, average, absLab,
                                      yAxisLims, plotsDir) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  metagenomicsData <- read_tsv(metagenomicDataFile)
  classData <- read_tsv(abClassesFile)
  normData <- read_tsv(normalisationFile)
  metaData <- read_csv(metaDataFile)

  classCol <- NULL

  metagenomicsData <- metagenomicsData[, c(TRUE, colSums(metagenomicsData[, -1]) > 0)]

  # Add antibiotic class to metagenomics table
  for (iData in 1:nrow(metagenomicsData))
  {
    abClass <- unique(classData$Drug_class[classData$Subject ==
                                      metagenomicsData$Annotation[iData]])
    classCol <- rbind(classCol, abClass)
  }

  rownames(classCol) <- 1:nrow(classCol)

  # Amalgamate antibiotic classes
  classCol[classCol == "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)"] <- "Rifamycin"
  classCol[classCol == "major facilitator superfamily (MFS) antibiotic efflux pump" ] <- "MDR"
  classCol[classCol == "resistance-nodulation-cell division (RND) antibiotic efflux pump"] <- "MDR"
  classCol[classCol == "aminocoumarin self resistant parY,aminocoumarin resistant parY"] <- "Aminocoumarin"
  classCol[classCol == "ampC-type beta-lactamase"] <- "Beta-lactam"
  classCol[classCol == "ADC beta-lactamase"] <- "Beta-lactam"
  classCol[classCol == "antibiotic resistant isoleucyl-tRNA synthetase (ileS)"] <- "Mupirocin"
  classCol[classCol == "triclosan"] <- "Triclosan"
  classCol[classCol == "Ampenicol"] <- "Amphenicol"
  classCol[classCol == "ATP-binding cassette (ABC) antibiotic efflux pump"] <- "MDR"

  drugClasses <- c("Aminocoumarin", "Aminoglycoside", "Amphenicol", "Beta-lactam",
                   "Fluoroquinolone", "Fosfomycin", "Glycopeptide", "MDR", "MLSB",
                   "Mupirocin", "Nucleoside", "Peptide antibiotic", "Rifamycin", "Sulfonamide",
                   "Tetracycline", "Triclosan", "Trimethoprim")
  metagenomicsData <- cbind(metagenomicsData, classCol)
  colnames(metagenomicsData)[ncol(metagenomicsData)] <- "DrugClass"

  drugClasses <- unique(c(drugClasses, metagenomicsData$DrugClass))
  metagenomicsData$DrugClass <- factor(metagenomicsData$DrugClass, levels = drugClasses)

  # For each sample add all the resistances in a particular class together
  metagenomicsAMR <- aggregate(metagenomicsData[ , 2:(ncol(metagenomicsData) - 1)],
                               by = list(metagenomicsData$DrugClass),
                               FUN = sum)

  colnames(metagenomicsAMR)[1] <- "ARG"

  # Normalise data by normType (passed in)
  normMGD <- metagenomicsAMR

  for (iSample in 2:(ncol(normMGD) - 1))
  {

    for (iData in 1:nrow(normMGD))
    {
      normMGD[iData, iSample] <- normMGD[iData, iSample] /
        normData[which(normData$Normalisation == normType),
                 which(colnames(normData) == colnames(normMGD)[iSample])]
    }

  }

  # Multiple data by multiplier (passed in)
  normMGD[, 2:ncol(normMGD)] <- multiplier * normMGD[, 2:(ncol(normMGD) - 1)]

  # Transpose data for ease of working - this then requires further adjustments
  normMGD <- t(normMGD)
  normMGD <- as.data.frame(normMGD)

  colnames(normMGD) <- normMGD[1, ]
  normMGD <- normMGD[2:nrow(normMGD), ]
  normMGD <- cbind(rownames(normMGD), normMGD)

  colnames(normMGD)[1] <- "SampleID"
  rownames(normMGD) <- 1:nrow(normMGD)
  # Ensure sample naming is consistent across data sources
  normMGD$SampleID <- sub("PAK0", "PAK_0", normMGD$SampleID)
  normMGD$SampleID <- sub("PHS0", "PHS_0", normMGD$SampleID)

  if (campaignType == "Env")
  {
    # Add sampling site and season
    normMGD <- addSamplingSite(normMGD, metaData)
    normMGD <- normMGD %>% relocate(SamplingSite, .after = SampleID)
    normMGD <- normMGD %>% relocate(Season, .after = SamplingSite)
    startCol <- 4
  }else
  {
    startCol <- 2
  }

  for (iCol in startCol:ncol(normMGD))
  {
    normMGD[, iCol] <- as.numeric(normMGD[, iCol])
  }

  if (campaignType == "Env")
  {
    # Split data into seasons
    normCampaign <- normMGD[normMGD$Season == "Summer Dry", ]
  }else
  {
    normCampaign <- normMGD
  }

  # aggregate  samples then pivot data
  if (average)
  {
    normMGDAverage <- aggregate(normCampaign[ , startCol:ncol(normCampaign)],
                                by = list(normCampaign$SamplingSite),
                                FUN = mean)
    colnames(normMGDAverage)[1] <- "SamplingSite"

    metagenomicsAMRLong <- pivot_longer(normMGDAverage, cols = c(2:ncol(normMGDAverage)),
                                        names_to = "ARG", values_to = "Abundance")
    xAxisLabels <- normMGDAverage$SamplingSite
  }else
  {
    metagenomicsAMRLong <- pivot_longer(normCampaign, cols = c(startCol:ncol(normCampaign)),
                                        names_to = "ARG", values_to = "Abundance")

    if (campaignType == "Env")
    {
      xAxisLabels <- normCampaign$SamplingSite
    }else
    {
      xAxisLabels <- normCampaign$SampleID
    }

  }

  metagenomicsAMRLong <- cbind(metagenomicsAMRLong[, 1], metagenomicsAMRLong)
  metagenomicsAMRLong$ARG <- factor(metagenomicsAMRLong$ARG, levels = drugClasses)
  colnames(metagenomicsAMRLong)[1] <- "xVals"

  if (campaignType == "Env")
  {
    plotTitleBase <- "Dry season"
  }else
  {
    plotTitleBase <- "Household"
  }

  # Plot dry season absolute data
  BarPlotMetagenomics(graphData = metagenomicsAMRLong, xAxisLabels = xAxisLabels,
                      absolute = TRUE, absLab = absLab, graphType = "ARGs", 
                      yAxisVals = yAxisLims, plotsDir = plotsDir,
                      plotTitle = paste(plotTitleBase, average, normType,
                                        "Normalised Absolute Abundances"))

  # Plot dry season relative data
  BarPlotMetagenomics(graphData = metagenomicsAMRLong, xAxisLabels = xAxisLabels,
                      absolute =  FALSE, absLab = NULL, graphType = "ARGs", 
                      plotsDir = plotsDir,
                      plotTitle = paste(plotTitleBase, average, normType,
                                        "Normalised Relative Abundances"))

  if (campaignType == "Env")
  {
    normWet <- normMGD[normMGD$Season == "Summer Wet", ]

    # aggregate  samples then pivot data
    if (average)
    {
      normMGDAverage <- aggregate(normWet[ , 4:ncol(normWet)],
                                  by = list(normWet$SamplingSite),
                                  FUN = mean)
      colnames(normMGDAverage)[1] <- "SamplingSite"

      metagenomicsAMRLong <- pivot_longer(normMGDAverage, cols = c(2:ncol(normMGDAverage)),
                                          names_to = "ARG", values_to = "Abundance")
    }else
    {
      normMGDAverage <- normWet
      metagenomicsAMRLong <- pivot_longer(normMGDAverage, cols = c(4:ncol(normMGDAverage)),
                                          names_to = "ARG", values_to = "Abundance")
    }

    metagenomicsAMRLong <- cbind(metagenomicsAMRLong[, 1], metagenomicsAMRLong)
    colnames(metagenomicsAMRLong)[1] <- "xVals"

    # Plot wet season absolute data
    BarPlotMetagenomics(graphData = metagenomicsAMRLong, xAxisLabels = normMGDAverage$SamplingSite,
                        absolute =  TRUE, absLab = absLab,
                        graphType = "ARGs", yAxisVals = yAxisLims, plotsDir = plotsDir,
                        plotTitle = paste("Wet season", average, normType,
                                          "Normalised Absolute Abundances"))

    # Plot wet season relative data
    BarPlotMetagenomics(graphData = metagenomicsAMRLong, xAxisLabels = normMGDAverage$SamplingSite,
                        absolute = FALSE, absLab = NULL,
                        graphType = "ARGs", plotsDir = plotsDir,
                        plotTitle = paste("Wet season", average, normType,
                                          "Normalised Relative Abundances"))
  }

}
