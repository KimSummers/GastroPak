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

ReadAndPlotGPHuBac <- function(metagenomicDataFile, abClassesFile, 
                               normalisationFile, huBacCovFile, metaDataFile, 
                               campaignType, normType, multiplier, average, 
                               absLab, yAxisLims, plotsDir) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  metagenomicsData <- read_tsv(metagenomicDataFile)
  classData <- read_tsv(abClassesFile)
  normData <- read_tsv(normalisationFile)
  metaData <- read_csv(metaDataFile)
  huBacData <- read_tsv(huBacCovFile)

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
  
  for (iSample in 2:(ncol(normMGD) - 1))
  {

    for (iData in 1:nrow(normMGD))
    {
      normMGD[iData, iSample] <- normMGD[iData, iSample] /
        normData[which(normData$Normalisation == normType),
                 which(colnames(normData) == colnames(normMGD)[iSample])]
    }

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
  
  
  # Multiply data by multiplier (passed in)
  normMGD[, 2:ncol(normMGD)] <- multiplier * normMGD[, 2:(ncol(normMGD) - 1)]

  # Transpose data for ease of working - this then requires further adjustments
  normMGD <- as.data.frame(t(normMGD))
  normHuBac <- as.data.frame(t(normHuBac))

  colnames(normMGD) <- normMGD[1, ]
  normMGD <- normMGD[2:nrow(normMGD), ]
  normMGD <- cbind(rownames(normMGD), normMGD)
  
  colnames(normHuBac) <- normHuBac[1, ]
  normHuBac <- normHuBac[2:nrow(normHuBac), ]
  normHuBac <- cbind(rownames(normHuBac), normHuBac)

  colnames(normMGD)[1] <- "SampleID"
  rownames(normMGD) <- 1:nrow(normMGD)
  # Ensure sample naming is consistent across data sources
  normMGD$SampleID <- sub("PAK0", "PAK_0", normMGD$SampleID)
  normMGD$SampleID <- sub("PHS0", "PHS_0", normMGD$SampleID)
  
  colnames(normHuBac)[1] <- "SampleID"
  rownames(normHuBac) <- 1:nrow(normHuBac)
  normHuBac$SampleID <- sub("PAK0", "PAK_0", normHuBac$SampleID)
  normHuBac$SampleID <- sub("PHS0", "PHS_0", normHuBac$SampleID)
  
  if (campaignType == "Env")
  {
    # Add sampling site and season
    normMGD <- addSamplingSite(normMGD, metaData)
    normMGD <- normMGD %>% relocate(SamplingSite, .after = SampleID)
    normMGD <- normMGD %>% relocate(Season, .after = SamplingSite)
    startCol <- 4
    
    normHuBac <- addSamplingSite(normHuBac, metaData)
    normHuBac <- normHuBac %>% relocate(SamplingSite, .after = SampleID)
    normHuBac <- normHuBac %>% relocate(Season, .after = SamplingSite)
  }else
  {
    startCol <- 2
  }

  for (iCol in startCol:ncol(normMGD))
  {
    normMGD[, iCol] <- as.numeric(normMGD[, iCol])
  }

  for (iCol in startCol:ncol(normHuBac))
  {
    normHuBac[, iCol] <- as.numeric(normHuBac[, iCol])
  }
  
  if (campaignType == "Env")
  {
    # Split data into seasons
    normCampaign <- normMGD[normMGD$Season == "Summer Dry", ]
    sites <- unique(normCampaign$SamplingSite)
    normHuBacCamp <- normHuBac[normHuBac$Season == "Summer Dry", ]
    normHuBacCamp <- normHuBacCamp[normHuBacCamp$SamplingSite %in% sites, ]
  }else
  {
    normCampaign <- normMGD
    normHuBacCamp <- normHuBac
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
    
    normBacAverage <- aggregate(normHuBacCamp[ , startCol:ncol(normHuBacCamp)],
                                by = list(normHuBacCamp$SamplingSite),
                                FUN = mean)
    colnames(normBacAverage)[1] <- "SamplingSite"
    metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(2:ncol(normBacAverage)),
                                        names_to = "Bac", values_to = "Abundance")
  }else
  {
    metagenomicsAMRLong <- pivot_longer(normCampaign, cols = c(startCol:ncol(normCampaign)),
                                        names_to = "ARG", values_to = "Abundance")

    metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(startCol:ncol(normBacAverage)),
                                        names_to = "Bac", values_to = "Abundance")

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
  BarPlotMetagenomics(graphData = metagenomicsAMRLong,  
                      graphData2 = metagenomicsBacLong, xAxisLabels = xAxisLabels, 
                      absolute = TRUE, absLab = absLab, graphType =  "Both", 
                      yAxisVals = yAxisLims, plotsDir = plotsDir,
                      plotTitle = paste(plotTitleBase, "HuBac", average, normType,
                                        "Normalised Absolute Abundances"))

  # Plot dry season relative data
  BarPlotMetagenomics(graphData = metagenomicsAMRLong, 
                      graphData2 = metagenomicsBacLong, xAxisLabels = xAxisLabels,
                      absolute =  FALSE, absLab = NULL, graphType = "Both", 
                      plotsDir = plotsDir,
                      plotTitle = paste(plotTitleBase, "HuBac", average, normType,
                                        "Normalised Relative Abundances"))

  if (campaignType == "Env")
  {
    normWet <- normMGD[normMGD$Season == "Summer Wet", ]
    normHuBacWet <- normHuBac[normHuBac$Season == "Summer Wet", ]
    sites <- unique(normWet$SamplingSite)
    normHuBacWet <- normHuBacWet[normHuBacWet$SamplingSite %in% sites, ]
    
    # aggregate  samples then pivot data
    if (average)
    {
      normMGDAverage <- aggregate(normWet[ , 4:ncol(normWet)],
                                  by = list(normWet$SamplingSite),
                                  FUN = mean)
      colnames(normMGDAverage)[1] <- "SamplingSite"

      metagenomicsAMRLong <- pivot_longer(normMGDAverage, cols = c(2:ncol(normMGDAverage)),
                                          names_to = "ARG", values_to = "Abundance")

      normBacAverage <- aggregate(normHuBacWet[ , 4:ncol(normHuBacWet)],
                                  by = list(normHuBacWet$SamplingSite),
                                  FUN = mean)
      colnames(normBacAverage)[1] <- "SamplingSite"
      
      metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(2:ncol(normBacAverage)),
                                          names_to = "Bac", values_to = "Abundance")
    }else
    {
      normMGDAverage <- normWet
      metagenomicsAMRLong <- pivot_longer(normMGDAverage, cols = c(4:ncol(normMGDAverage)),
                                          names_to = "ARG", values_to = "Abundance")
      normBacAverage <- normHuBacWet
      metagenomicsBacLong <- pivot_longer(normBacAverage, cols = c(4:ncol(normBacAverage)),
                                          names_to = "Bac", values_to = "Abundance")
    }

    metagenomicsAMRLong <- cbind(metagenomicsAMRLong[, 1], metagenomicsAMRLong)
    colnames(metagenomicsAMRLong)[1] <- "xVals"

    metagenomicsBacLong <- cbind(metagenomicsBacLong[, 1], metagenomicsBacLong)
    colnames(metagenomicsBacLong)[1] <- "xVals"
    
    # Plot wet season absolute data
    BarPlotMetagenomics(graphData = metagenomicsAMRLong, graphData2 = metagenomicsBacLong, 
                        xAxisLabels = normMGDAverage$SamplingSite,
                        absolute =  TRUE, absLab = absLab, graphType = "Both",
                        yAxisVals = yAxisLims, plotsDir = plotsDir,
                        plotTitle = paste("Wet 2022 HuBac", average, normType,
                                          "Normalised Absolute Abundances"))

    # Plot wet season relative data
    BarPlotMetagenomics(graphData = metagenomicsAMRLong, graphData2 = metagenomicsBacLong, 
                        xAxisLabels = normMGDAverage$SamplingSite,
                        absolute = FALSE, absLab = NULL, graphType = "Both", 
                        plotsDir = plotsDir,
                        plotTitle = paste("Wet 2022 HuBac", average, normType,
                                          "Normalised Relative Abundances"))
  }

}
