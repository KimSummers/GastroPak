# ReadAndPlotGPMetaPhlan
#
# Read in metagenomic metaPhlan data and bar chart plot it
#
# file ReadAndPlotGPMetaPhlan
#
# inputs
# 	speciesFile           - file containing species data
# 	phylaFile             - file containing phylogeny data
# 	metaDataFile          - file containing metadata
#   campaignType          - "Env" or "HH" for environmental or household
#   average               - TRUE = average data across samples, FALSE do not
#   otherPercent          - Species below this threshold are amalgamated
#   absLab                - Label for the y-axis for absolute data plots
#   yAxisLims             - Limits for the y-axis
#   plotsDir              - directory to store plots in
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  02/05/24  Wellington Lab - School of Life Sciences - University of Warwick

ReadAndPlotGPMetaPhlan <- function(speciesFile, phylaFile, metaDataFile,
                                   campaignType, average, otherPercent, absLab,
                                   yAxisLims, plotsDir) {

  library(tidyverse)

  # Read in data and put it in a dataframe
  species <- read.csv(speciesFile, sep = "\t")
  phylas <- read.csv(phylaFile, sep = "\t", header = FALSE)
  metaData <- read_csv(metaDataFile)

  colnames(species) <- species[1, ]
  species <- species[-1, ]

  # subset species data to only include relevant data
  if (campaignType == "Env")
  {
    species <- species[, c(1, which(substr(colnames(species), 1, 3) == "PAK"))]
  }else
  {
    species <- species[, c(1, which(substr(colnames(species), 1, 3) == "PHS"))]
  }

  species <- species[, order(colnames(species))]
  species <- as.data.frame(species)

  # ensure data is numeric
  for (iCol in 2:ncol(species))
  {
    species[, iCol] <- as.numeric(species[, iCol])
  }

  # remove samples for which there
  species <- species[, c(TRUE, colSums(species[, -1]) > 0)]

  familyCol <- NULL

  for (iData in 1:nrow(species))
  {
    bactFam <- phylas[phylas[, 1] == species$clade_name[iData], 6]
    familyCol <- rbind(familyCol, bactFam)
  }

  species <- cbind(species, familyCol)
  colnames(species)[ncol(species)] <- "Family"
  species$Familiy <- substr(species$Familiy, 4, nchar(species$Familiy))

  # For each sample add all the species in a particular families together
  species <- aggregate(species[ , 2:(ncol(species) - 1)],
                       by = list(species$Familiy), FUN = sum)

  # Transpose data for ease of working - this then requires further adjustments
  species <- t(species)
  species <- as.data.frame(species)

  colnames(species) <- species[1, ]
  species <- species[2:(nrow(species) - 1), ]

  species <- species %>% rownames_to_column("SampleID")

  colnames(species)[1] <- "SampleID"
  # Ensure sample naming is consistent across data sources
  species$SampleID <- sub("PAK0", "PAK_0", species$SampleID)
  species$SampleID <- sub("PHS0", "PHS_0", species$SampleID)

  if (campaignType == "Env")
  {
    # Add sampling site and season
    species <- addSamplingSite(species, metaData)
    species <- species %>% relocate(SamplingSite, .after = SampleID)
    species <- species %>% relocate(Season, .after = SamplingSite)
    startCol <- 4
  }else
  {
    startCol <- 2
  }

  for (iCol in startCol:ncol(species))
  {
    species[, iCol] <- as.numeric(species[, iCol])
  }

  keepCols <- NULL

  for (iCol in 4:ncol(species))
  {

    if (max(species[, iCol]) > otherPercent)
    {
      keepCols <- cbind(keepCols, TRUE)
    }else
    {
      keepCols <- cbind(keepCols, FALSE)
    }

  }

  # ncol(species[, cbind(TRUE, TRUE, TRUE, keepCols)])

  otherFamily <- NULL

  for (iRow in 1:nrow(species))
  {
    otherRow <- sum(species[iRow, cbind(FALSE, FALSE, FALSE, !keepCols)])
    otherFamily <- rbind(otherFamily, otherRow)
  }

  speciesSub <- cbind(species[, cbind(TRUE, TRUE, TRUE, keepCols)], otherFamily)
  colnames(speciesSub)[ncol(speciesSub)] <- "Other"

  if (campaignType == "Env")
  {
    # Split data into seasons
    speciesCampaign <- speciesSub[speciesSub$Season == "Summer Dry", ]
  }else
  {
    speciesCampaign <- speciesSub
  }

  # aggregate  samples then pivot data
  if (average)
  {
    speciesAverage <- aggregate(speciesCampaign[ , startCol:ncol(speciesCampaign)],
                                by = list(speciesCampaign$SamplingSite),
                                FUN = mean)
    colnames(speciesAverage)[1] <- "SamplingSite"

    metagenomicsSpeciesLong <- pivot_longer(speciesAverage, cols = c(2:ncol(speciesAverage)),
                                            names_to = "Family", values_to = "Abundance")
    xAxisLabels <- speciesAverage$SamplingSite
  }else
  {
    metagenomicsSpeciesLong <- pivot_longer(speciesCampaign, cols = c(startCol:ncol(speciesCampaign)),
                                            names_to = "Family", values_to = "Abundance")

    if (campaignType == "Env")
    {
      xAxisLabels <- speciesCampaign$SamplingSite
    }else
    {
      xAxisLabels <- speciesCampaign$SampleID
    }

  }

  metagenomicsSpeciesLong <- cbind(metagenomicsSpeciesLong[, 1], metagenomicsSpeciesLong)
  metagenomicsSpeciesLong$Family <- factor(metagenomicsSpeciesLong$Family,
                                           levels = unique(metagenomicsSpeciesLong$Family))
  colnames(metagenomicsSpeciesLong)[1] <- "xVals"

  if (campaignType == "Env")
  {
    plotTitleBase <- "Dry season"
  }else
  {
    plotTitleBase <- "Household"
  }

  # Plot dry season data
  BarPlotMetagenomicsSpecies(graphData = metagenomicsSpeciesLong,
                             xAxisLabels = xAxisLabels, absolute = FALSE,
                             absLab = absLab, huBac = FALSE,
                             yAxisVals = yAxisLims, plotsDir = plotsDir,
                             plotTitle = paste(plotTitleBase, average, "Relative Abundances"))

  if (campaignType == "Env")
  {
    # Split data into seasons
    speciesCampaign <- speciesSub[speciesSub$Season == "Summer Wet", ]

      # aggregate  samples then pivot data
    if (average)
    {
      speciesAverage <- aggregate(speciesCampaign[ , startCol:ncol(speciesCampaign)],
                                  by = list(speciesCampaign$SamplingSite),
                                  FUN = mean)
      colnames(speciesAverage)[1] <- "SamplingSite"

      metagenomicsSpeciesLong <- pivot_longer(speciesAverage, cols = c(2:ncol(speciesAverage)),
                                              names_to = "Family", values_to = "Abundance")
      xAxisLabels <- speciesAverage$SamplingSite
    }else
    {
      metagenomicsSpeciesLong <- pivot_longer(speciesCampaign, cols = c(startCol:ncol(speciesCampaign)),
                                              names_to = "Family", values_to = "Abundance")
      xAxisLabels <- speciesCampaign$SamplingSite
    }


    metagenomicsSpeciesLong <- cbind(metagenomicsSpeciesLong[, 1], metagenomicsSpeciesLong)
    metagenomicsSpeciesLong$Family <- factor(metagenomicsSpeciesLong$Family,
                                             levels = unique(metagenomicsSpeciesLong$Family))
    colnames(metagenomicsSpeciesLong)[1] <- "xVals"

    # Plot wet season data
    BarPlotMetagenomicsSpecies(graphData = metagenomicsSpeciesLong,
                               xAxisLabels = xAxisLabels, absolute =  FALSE,
                               absLab = absLab, huBac = FALSE,
                               yAxisVals = yAxisLims, plotsDir = plotsDir,
                               plotTitle = paste("Wet season", average, "Relative Abundances"))

  }

}
