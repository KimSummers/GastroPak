# GPMetagenomicsNMDSGenus
#
# NDMS analysis of Genus data for GastroPak - code adapted from code supplied by Emma Travis
#
# file GPMetagenomicsNMDSGenus
#
# inputs
# 	metagenomicDataFile   - file containing plate count data
# 	metaDataFile          - file containing metadata
#   riversFile            - file to conver site into river name and stream position
# 	phylaFile             - file containing phylogeny data
#.  dataDir               - directory to save data in
# 	categories            - which categories e.g. seasons should data be split into
# 	sampleOrder           - Desired order for samples
#   campaignType          - "Env" or "HH" for environmental or household
#   plotsDir              - directory to store plots in

# 	normalisationFile     - file containing normalisation data
#   normType              - type of normalisation to use
#   multiplier            - value to multiple data by
#   average               - TRUE = average data across samples, FALSE do not
#   absLab                - Label for the y-axis for absolute data plots
#   yAxisLims             - Limits for the y-axis
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  26/03/24  Wellington Lab - School of Life Sciences - University of Warwick

GPMetagenomicsNMDSGenus <- function(metagenomicDataFile, metaDataFile, riversFile,
                                    phylaFile, distMatrixFile, categories, sampleOrder,
                                    campaignType, plotsDir) {
  # Load necessary packages
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(pals)

  # Read in data
  mgGenusWide <- read.csv(metagenomicDataFile, sep = "\t")
  metaDataWide <- read.csv(metaDataFile)
  phylas <- read.csv(phylaFile, sep = "\t", header = FALSE)

  colnames(phylas) <- c("Species", "Kingdom", "Phyla", "Class", "Order",
                        "Family", "Genus", "species")
  phylas$Genus <- sub("g__", "", phylas$Genus)
  phylas$Class <- sub("c__", "", phylas$Class)

  # Define the desired order of levels and corresponding labels for treatment
  desiredOrder <- c("Dry", "Wet")
  newLabels <- c("Dry 2022", "Wet 2022")

  # Define the desired order of levels and corresponding labels for River_loc
  # desired_order1 <- c("Aire_up", "Aire_down", "Sowe_up", "Sowe_down")
  # new_labels1 <- c("River Aire upstream", "River Aire downstream", "River Sowe upstream", "River Sowe downstream")
  #

  metaDataWide <- metaDataWide[metaDataWide$Sample_ID != "", ]
  riverSites <- unique(metaDataWide$Site_ID)
  riversData <- read_csv(file = riversFile)
  colnames(metaDataWide)[which(colnames(metaDataWide) == "Site_ID")] <- "SamplingSite"

  metaDataWide <- addRivers(countData = metaDataWide, metaData = riversData)
  gpRivers <- unique(metaDataWide$River)
  locations <- unique(metaDataWide$RiverLocation)

  # Convert categorical variable to factor with desired order of levels and new labels
  metaDataWide$Campaign <- factor(metaDataWide$Season_WD, levels = desiredOrder,
                                  labels = newLabels)

  metaDataWide$SamplingSite <- factor(metaDataWide$SamplingSite, levels = riverSites,
                                      labels = riverSites)
  metaDataWide$River <- factor(metaDataWide$River, levels = gpRivers)
  metaDataWide$RiverLocation <- factor(metaDataWide$RiverLocation, levels = locations)
  metaDataWide$CampaignLoc <- paste(metaDataWide$Campaign, metaDataWide$RiverLocation)
  campaignLocs <- unique(metaDataWide$CampaignLoc)
  metaDataWide$CampaignLoc <- factor(metaDataWide$CampaignLoc, levels = campaignLocs)

  # preparing Genus data
  mgGenusWide[is.na(mgGenusWide)] <- 0
  # colnames(mgGenusWide) <- mgGenusWide[1, ]
  mgGenusWide <- mgGenusWide[2:nrow(mgGenusWide), ]
  rownames(mgGenusWide) <- 1:nrow(mgGenusWide)
  mgGenusWide <- mgGenusWide %>% column_to_rownames("clade_name")

  for (iCol in 1:ncol(mgGenusWide))
  {
    mgGenusWide[, iCol] <- as.numeric(mgGenusWide[, iCol])
  }

  mgGenusWide <- mgGenusWide[, colSums(mgGenusWide) > 0]

  # subset species data to only include relevant data
  if (campaignType == "Env")
  {
    mgGenusWide <- mgGenusWide[, which(substr(colnames(mgGenusWide), 1, 3) == "PAK")]
  }else
  {
    mgGenusWide <- mgGenusWide[, which(substr(colnames(mgGenusWide), 1, 3) == "PHS")]
  }

  mgGenusWide <- mgGenusWide[, order(colnames(mgGenusWide))]

  mgGenusWide <- as.data.frame(t(mgGenusWide))
  rownames(mgGenusWide) <- sub("PAK0", "PAK_0", rownames(mgGenusWide))

  metaDataWide <- metaDataWide[metaDataWide$Barcode %in% rownames(mgGenusWide), ]

  # transform to relative abundance
  # hellinger might be better- square root of total and accounts for more high and low numbers
  genus.rel <- decostand(mgGenusWide, method = "total")

  # Calculate distance matrix - make distance matrix data frame - write .csv
  genusDistMat <- vegdist(genus.rel, method = "bray")
  genusDistanceMatrix <- as.matrix(genusDistMat, labels = T)
  write.csv(genusDistanceMatrix, paste(plotsDir, "GenusDistMat.csv", sep = ""))

  # Running NMDS in vegan (metaMDS)
  genusNMS <- metaMDS(genusDistMat, distance = "bray", k = 2, maxit = 999,
                    trymax = 250, autotransform = FALSE)
  genusScores <- `sppscores<-`(genusNMS, mgGenusWide)
  genusCorrelation <- cor(mgGenusWide, genusNMS$points, use = "complete.obs",
                        method = "pearson")
  write.csv(genusCorrelation, file = paste(plotsDir, "genus_PearsonCor.csv", sep = ""))

  # Shepards test/goodness of fit
  fitGoodness <- goodness(genusNMS)
  stressplot(genusNMS)

  pdf(file = paste(plotsDir, "GenusStressPlot.pdf", sep = ""))
  stressplot(genusNMS)
  dev.off()

  # Extract NMDS coordinates
  nmdsCoordinates <- scores(genusNMS, choices = c(1, 2))

  NMDSGenus <- cbind(metaDataWide, nmdsCoordinates)

  genusSppFit <- envfit(genusNMS, genus.rel, permutations = 999)

  # obtain all info from NMDS as dataframe
  genusScores <- as.data.frame(scores(genusScores, display = "species"))
  genusFit <- as.data.frame(scores(genusSppFit, display = "vectors"))
  genusScores <- cbind(genusScores, Species = rownames(genusScores))

  # add species names to dataframe
  genuspVal <- as.data.frame(genusSppFit$vectors$pvals)

  # add pvalues to dataframe so you can select species which are significant
  genusScores <- cbind(genusScores, pval = genuspVal$`genusSppFit$vectors$pvals`)

  # add metadata information to the genus_count file
  genusScores <- inner_join(genusScores, phylas, by = "Species")

  #abbreviate species names
  #spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6))

  #subset data to show species significant at 0.01
  sig.genus.scrs <- subset(genusScores, pval <= 0.01)

  #calculate vector (distance from 0,0)

  sig.genus.scrs$distance <- sqrt(sig.genus.scrs$NMDS1 ^ 2 + sig.genus.scrs$NMDS2 ^ 2)

  top20Sig.genus.scrs<- head(sig.genus.scrs[order(-sig.genus.scrs$distance), ], 20)

  #sigVars <- genus_spp_fit$signif
  graphColours <- c("deepskyblue", "deepskyblue4", "coral",
                    "coral4", "red", "blue", "chartreuse3",
                    "lightblue", "cornflowerblue", "darkslategray",
                    "cadetblue4", "orange", "darkolivegreen", "grey",
                    "darkblue", "purple", "pink", "yellow", "cyan",
                    "aquamarine3", "peachpuff")
  graphColours2 <- c("cyan2", "blue", "lightcyan3",
                     "lightpink2", "burlywood1", "darkorange3", "chartreuse3",
                     "lightblue", "cornflowerblue","darkslategray",
                     "cadetblue4", "orange", "darkolivegreen", "grey",
                     "darkblue", "purple", "pink", "yellow", "cyan",
                     "aquamarine3")

  # nmds colored by river location and shape by season
  #sets up the plot
  nmds.plot.genus <- ggplot(NMDSGenus, aes(x = NMDS1, y = NMDS2))
  #adds site points to plot, shape determined by month, colour determined by site
  nmds.plot.genus <- nmds.plot.genus + geom_point(aes(NMDS1, NMDS2, colour = factor(River),
                                                  shape = factor(CampaignLoc)), size = 4)
  nmds.plot.genus <- nmds.plot.genus + stat_ellipse(geom = "polygon", type = "t", alpha = 0.1,
                                                mapping = aes(fill = Campaign, color = Campaign))
  # add legend labels for Location and Month
  nmds.plot.genus <- nmds.plot.genus + labs(color = "River", shape = "River Location")
  # add legend labels for Location and Month
  nmds.plot.genus <- nmds.plot.genus + ggtitle("Species NMDS - Grouping by season")
  nmds.plot.genus <- nmds.plot.genus + theme_bw() + guides(fill = FALSE)
  # Set colors
  # nmds.plot.genus <- nmds.plot.genus + scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral",
  #                                 "coral4", "red","blue", "chartreuse3",
  #                                 "lightblue","cornflowerblue","darkslategray",
  #                                 "cadetblue4","orange","darkolivegreen","grey",
  #                                 "darkblue"))
  nmds.plot.genus <- nmds.plot.genus + scale_color_manual(values = graphColours2)
  nmds.plot.genus <- nmds.plot.genus + scale_fill_manual(values = graphColours2[6:19])
  #  nmds.plot.genus <- nmds.plot.genus + scale_color_manual(values = glasbey())
  # Set shapes (19: circle, 17: triangle, 15: square)
  nmds.plot.genus <- nmds.plot.genus + scale_shape_manual(values = c(19, 17, 15, 1, 2, 0))
  nmds.plot.genus <- nmds.plot.genus + guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

  nmds.plot.genus

  ggsave(paste(plotsDir, "All Genus NMDS by season.pdf", sep = ""), nmds.plot.genus,
         height = 6, width = 8)

  nmds.plot.genus2 <- ggplot(NMDSGenus, aes(x = NMDS1, y = NMDS2))
  #adds site points to plot, shape determined by month, colour determined by site
  nmds.plot.genus2 <- nmds.plot.genus2 + geom_point(aes(NMDS1, NMDS2, colour = factor(River),
                                                    shape = factor(CampaignLoc)), size = 4)
  nmds.plot.genus2 <- nmds.plot.genus2 + stat_ellipse(geom = "polygon", type = "t", alpha = 0.1,
                                                  mapping = aes(fill = River, color = River))
  # add legend labels for Location and Month
  nmds.plot.genus2 <- nmds.plot.genus2 + labs(color = "River", shape = "Location")
  # add legend labels for Location and Month
  nmds.plot.genus2 <- nmds.plot.genus2 + ggtitle("Antibiotic resistance gene NMDS - Grouping by location")
  nmds.plot.genus2 <- nmds.plot.genus2 + theme_bw() + guides(fill = FALSE)
  # Set colors
  # nmds.plot.genus2 <- nmds.plot.genus2 + scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral",
  #                                   "coral4", "red","blue", "chartreuse3",
  #                                   "lightblue","cornflowerblue","darkslategray",
  #                                   "cadetblue4","orange","darkolivegreen","grey",
  #                                   "darkblue"))
  nmds.plot.genus2 <- nmds.plot.genus2 + scale_color_manual(values = graphColours2)
  nmds.plot.genus2 <- nmds.plot.genus2 + scale_fill_manual(values = graphColours2)
  # Set shapes (19: circle, 17: triangle, 15: square)
  nmds.plot.genus2 <- nmds.plot.genus2 + scale_shape_manual(values = c(19, 17, 15, 1, 2, 0))
  nmds.plot.genus2 <- nmds.plot.genus2 + guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

  nmds.plot.genus2

  ggsave(paste(plotsDir, "All Genus NMDS by river.pdf", sep = ""), nmds.plot.genus2,
         height = 6, width = 8)

  #ADD significant species (i.e. genus driving diversity)
  nmds.plot.genus3 <- nmds.plot.genus +
    geom_segment(data = sig.genus.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)
  #add vector arrows of significant species
  nmds.plot.genus3 <- nmds.plot.genus3 +
    ggrepel::geom_text_repel(data = sig.genus.scrs, aes(x = NMDS1, y = NMDS2,
                                                      label = Genus),
                             cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.genus3
  #too many genus driving diversity so how to reduce this? Reduced sig to <0.01 not <0.05 - maybe also use size of vector??! / combine CTX-M etc?

  ggsave(paste(plotsDir, "All Genus NMDS by season plus vectors.pdf", sep = ""),
         nmds.plot.genus3, height = 6, width = 8)

  #ADD significant species (i.e. genus driving diversity)
  nmds.plot.genus4 <- nmds.plot.genus +
    geom_segment(data = top20Sig.genus.scrs, aes(x = 0, xend = NMDS1, y = 0,
                                               yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)

  #add vector arrows of significant species
  nmds.plot.genus4 <- nmds.plot.genus4 +
    ggrepel::geom_text_repel(data = top20Sig.genus.scrs,
                             aes(x = NMDS1, y = NMDS2, label = Genus), cex = 3,
                             direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.genus4

  ggsave(paste(plotsDir, "All Genus NMDS by season plus top 20 vectors.pdf", sep = ""),
         nmds.plot.genus4, height = 6, width = 8)

  # link the name of genus to drug_class etc

  #ADD significant species (i.e. classes driving diversity)
  nmds.plot.genus5 <- nmds.plot.genus +
    geom_segment(data = top20Sig.genus.scrs,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, colour = Class),
                 arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3)
  # # add vector arrows of significant species
  # nmds.plot.genus5 <- nmds.plot.genus5 +
  #   ggrepel::geom_text_repel(data = sig.genus.scrs,
  #                            aes(x = NMDS1, y = NMDS2, label = Class, colour = Class),
  #                            cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.genus5

  ggsave(paste(plotsDir, "All Genus NMDS by season plus classes.pdf", sep = ""), nmds.plot.genus5,
         height = 6, width = 8)

  #ADD significant species (i.e. genus driving diversity)
  nmds.plot.genus6 <- nmds.plot.genus2 +
    geom_segment(data = sig.genus.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)
  #add vector arrows of significant species
  nmds.plot.genus6 <- nmds.plot.genus6 +
    ggrepel::geom_text_repel(data = sig.genus.scrs, aes(x = NMDS1, y = NMDS2,
                                                        label = Genus),
                             cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.genus6
  #too many genus driving diversity so how to reduce this? Reduced sig to <0.01 not <0.05 - maybe also use size of vector??! / combine CTX-M etc?

  ggsave(paste(plotsDir, "NMDS Genus All main by river plus vectors.pdf", sep = ""), nmds.plot.genus6)

  #ADD significant species (i.e. genus driving diversity)
  nmds.plot.genus7 <- nmds.plot.genus2 +
    geom_segment(data = top20Sig.genus.scrs, aes(x = 0, xend = NMDS1, y = 0,
                                                 yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)

  #add vector arrows of significant species
  nmds.plot.genus7 <- nmds.plot.genus7 +
    ggrepel::geom_text_repel(data = top20Sig.genus.scrs,
                             aes(x = NMDS1, y = NMDS2, label = Genus), cex = 3,
                             direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.genus7

  ggsave(paste(plotsDir, "NMDS Genus All main by river plus top 20 vectors.pdf", sep = ""), nmds.plot.genus7)

  # link the name of genus to drug_class etc

  #ADD significant species (i.e. classes driving diversity)
  nmds.plot.genus8 <- nmds.plot.genus2 +
    geom_segment(data = top20Sig.genus.scrs,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, colour = Class),
                 arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3)
  # # add vector arrows of significant species
  # nmds.plot.genus5 <- nmds.plot.genus5 +
  #   ggrepel::geom_text_repel(data = sig.genus.scrs,
  #                            aes(x = NMDS1, y = NMDS2, label = Class, colour = Class),
  #                            cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.genus8

  ggsave(paste(plotsDir, "NMDS Genus All main by river plus classes.pdf", sep = ""), nmds.plot.genus8)

  # # Filter out a particular class, for example, class 'A'
  # filt_sig.genus.scrs.meta <- subset(sig.genus.scrs.meta, drug_class != "MDR efflux pump")
  #
  # nmds.plot.genus9 <- nmds.plot.genus +
  #   geom_segment(data = filt_sig.genus.scrs.meta,
  #                aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, colour = Drug_class),
  #                arrow = arrow(length = unit(0.25, "cm")),
  #                lwd = 0.3)
  # # add vector arrows of significant species
  # nmds.plot.genus9 <- nmds.plot.genus9 +
  #   ggrepel::geom_text_repel(data = filt_sig.genus.scrs.meta,
  #                            aes(x = NMDS1, y = NMDS2, label = Species, colour = Drug_class),
  #                            cex = 3, direction = "both", segment.size = 0.25)
  # # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  #
  # nmds.plot.genus9
  #
  ## -------Perform PERMANOVA to check significance of clustering from ordination
  # Produce distance matrix
  # Load data
  # Assuming row names are sample IDs in both data.matrix and sample_metadata

  # Calculating distance matrix using vegan function vegdist
  bray_dist_mat_genus <- vegdist(mgGenusWide, method = "bray")
  genus_dist <- as.data.frame(genusDistanceMatrix)
  genus_dist2 <- cbind(barcode = rownames(genus_dist), genus_dist)
  colnames(metaDataWide)[which(colnames(metaDataWide) == "Barcode")] <- "barcode"

  meta_dist <- inner_join(metaDataWide, genus_dist2, by = "barcode")

  #permutations for tests
  perm <- 999
  test1 <- adonis2(genus_dist~River, data = meta_dist, permutations = perm,
                   method = euclidean)
  write.csv(test1, file = paste(plotsDir, " Genus river sig.csv", sep = ""))
  
  # sig p = 0.001
  test2 <- adonis2(genus_dist~River*Campaign, data = meta_dist, permutations = perm,
                   method=euclidean)
  write.csv(test2, file = paste(plotsDir, " Genus river and campaign sig.csv", sep = ""))

  #River p = 0.001, Campaign p = 0.001, River * Campaign p = 0.001
  test3 <- adonis2(genus_dist~Campaign, data = meta_dist, permutations=perm,
                   method = euclidean)
  write.csv(test3, file = paste(plotsDir, " Genus campaign sig.csv", sep = ""))
  # sig p=0.001

  test8 <- adonis2(genus_dist~RiverLocation*Campaign*River, data = meta_dist, permutations=perm,
                   method = euclidean)
  write.csv(test8, file = paste(plotsDir, " Genus location, campaign and river sig.csv", sep = ""))
  
  # Run pairwise adonis tests (change location/month to check for each)
  # pairwise_adonis_results <- pairwise.perm.manova(genus_dist, wimp_wide_env$River_loc, nperm = 999, p.method = "bonferroni")
  # pairwise perm manova not working! can't find function
  # View the results
  # print(pairwise_adonis_results)

}

