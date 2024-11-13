# GPMetagenomicsNMDSGenusHH
#
# NDMS analysis of Genus data for GastroPak - code adapted from code supplied by Emma Travis
#
# file GPMetagenomicsNMDSGenusHH
#
# inputs
# 	metagenomicDataFile   - file containing plate count data
# 	metaDataFile          - file containing metadata
# 	phylaFile             - file containing phylogeny data
# 	categories            - which categories e.g. seasons should data be split into
# 	sampleOrder           - Desired order for samples
#   plotsDir              - directory to store plots in
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  02/05/24  Wellington Lab - School of Life Sciences - University of Warwick

GPMetagenomicsNMDSGenusHH <- function(metagenomicDataFile, metaDataFile, phylaFile,
                                      categories, sampleOrder, plotsDir) {
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


  # name columns and remove prefixes from data
  colnames(phylas) <- c("Species", "Kingdom", "Phyla", "Class", "Order",
                        "Family", "Genus", "species")
  phylas$Genus <- sub("g__", "", phylas$Genus)
  phylas$Class <- sub("c__", "", phylas$Class)

  colnames(metaDataWide)[which(colnames(metaDataWide) == "Sample.ID")] <- "Sample_ID"
  metaDataWide <- metaDataWide[metaDataWide$Sample_ID != "", ]

  households <- unique(metaDataWide$Household)
  locations <- unique(metaDataWide$Location)

  metaDataWide$Household <- factor(metaDataWide$Household, levels = households)
  metaDataWide$Location <- factor(metaDataWide$Location, levels = locations)

  # preparing Genus data
  mgGenusWide[is.na(mgGenusWide)] <- 0
  colnames(mgGenusWide) <- mgGenusWide[1, ]
  mgGenusWide <- mgGenusWide[2:nrow(mgGenusWide), ]
  rownames(mgGenusWide) <- 1:nrow(mgGenusWide)
  mgGenusWide <- mgGenusWide %>% column_to_rownames("clade_name")

  # ensure data is numeric
  for (iCol in 1:ncol(mgGenusWide))
  {
    mgGenusWide[, iCol] <- as.numeric(mgGenusWide[, iCol])
  }

  # remove samples which have 0 for all genii
  mgGenusWide <- mgGenusWide[, colSums(mgGenusWide) > 0]
  # only want household samples
  mgGenusWide <- mgGenusWide[, which(substr(colnames(mgGenusWide), 1, 3) == "PHS")]

  # order the columns
  mgGenusWide <- mgGenusWide[, order(colnames(mgGenusWide))]

  # want to have samples in rows and genus in columns
  mgGenusWide <- as.data.frame(t(mgGenusWide))
  rownames(mgGenusWide) <- sub("PHS0", "PHS_0", rownames(mgGenusWide))

  # only want metadata for these samples
  metaDataWide <- metaDataWide[metaDataWide$Sample_ID %in% rownames(mgGenusWide), ]

  # transform to relative abundance
  # hellinger might be better- square root of total and accounts for more high and low numbers
  genus.rel <- decostand(mgGenusWide, method = "total")

  # Calculate distance matrix - make distance matrix data frame - write .csv
  genusDistMat <- vegdist(genus.rel, method = "bray")
  genusDistanceMatrix <- as.matrix(genusDistMat, labels = T)
  write.csv(genusDistanceMatrix, paste(dataDir, "HH GenusDistMat.csv", sep = ""))

  # Running NMDS in vegan (metaMDS)
  genusNMS <- metaMDS(genusDistMat, distance = "bray", k = 2, maxit = 999,
                      trymax = 250, autotransform = FALSE)
  genusScores <- `sppscores<-`(genusNMS, mgGenusWide)
  genusCorrelation <- cor(mgGenusWide, genusNMS$points, use = "complete.obs",
                          method = "pearson")
  write.csv(genusCorrelation, file = paste(dataDir, "HH genus_PearsonCor.csv", sep = ""))

  # Shepards test/goodness of fit
  fitGoodness <- goodness(genusNMS)
  stressplot(genusNMS)

  pdf(file = paste(plotsDir, "HH GenusStressPlot.pdf", sep = ""))
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

  # nmds colored by river location and shape by season
  #sets up the plot
  nmds.plot.genus <- ggplot(NMDSGenus, aes(x = NMDS1, y = NMDS2))
  #adds site points to plot, shape determined by month, colour determined by site
  nmds.plot.genus <- nmds.plot.genus + geom_point(aes(NMDS1, NMDS2, colour = factor(Household),
                                                      shape = factor(Location)), size = 4)
  nmds.plot.genus <- nmds.plot.genus + stat_ellipse(geom = "polygon", type = "t", alpha = 0.1,
                                                    mapping = aes(fill = Location, color = Location))
  # add legend labels for Location and Month
  nmds.plot.genus <- nmds.plot.genus + labs(color = "Location", shape = "Household")
  # add legend labels for Location and Month
  nmds.plot.genus <- nmds.plot.genus + ggtitle("Species NMDS - Grouping by location")
  nmds.plot.genus <- nmds.plot.genus + theme_bw() + guides(fill = FALSE)
  # Set colors
  # nmds.plot.genus <- nmds.plot.genus + scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral",
  #                                 "coral4", "red","blue", "chartreuse3",
  #                                 "lightblue","cornflowerblue","darkslategray",
  #                                 "cadetblue4","orange","darkolivegreen","grey",
  #                                 "darkblue"))
  nmds.plot.genus <- nmds.plot.genus + scale_color_manual(values = graphColours)
  nmds.plot.genus <- nmds.plot.genus + scale_fill_manual(values = graphColours)
  #  nmds.plot.genus <- nmds.plot.genus + scale_color_manual(values = glasbey())
  # Set shapes (19: circle, 17: triangle, 15: square)
  nmds.plot.genus <- nmds.plot.genus + scale_shape_manual(values = c(19, 17, 15, 1, 2, 0))
  nmds.plot.genus <- nmds.plot.genus + guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

  nmds.plot.genus

  ggsave(paste(plotsDir, "NMDS HH Genus All main.pdf", sep = ""), nmds.plot.genus)

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

  ggsave(paste(plotsDir, "NMDS HH Genus All main plus vectors.pdf", sep = ""), nmds.plot.genus3)

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

  ggsave(paste(plotsDir, "NMDS HH Genus All main plus top 20 vectors.pdf", sep = ""), nmds.plot.genus4)

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

  ggsave(paste(plotsDir, "NMDS HH Genus All main plus classes.pdf", sep = ""), nmds.plot.genus5)

  # Calculating distance matrix using vegan function vegdist
  bray_dist_mat_genus <- vegdist(mgGenusWide, method = "bray")
  genus_dist <- as.data.frame(genusDistanceMatrix)
  genus_dist2 <- cbind(barcode = rownames(genus_dist), genus_dist)
  metaDataWide$barcode <- metaDataWide$Sample_ID

  meta_dist <- inner_join(metaDataWide, genus_dist2, by = "barcode")

  #permutations for tests
  perm <- 999
  test1 <- adonis2(genus_dist~Household, data = meta_dist, permutations = perm,
                   method = euclidean)
  write.csv(test1, file = paste(dataDir, "Genus HH household sig.csv", sep = ""))

  # sig p = 0.001
  test2 <- adonis2(genus_dist~Household*Location, data = meta_dist, permutations = perm,
                   method=euclidean)
  write.csv(test2, file = paste(dataDir, " Genus HH household and location sig.csv", sep = ""))

  #River p = 0.001, Campaign p = 0.001, River * Campaign p = 0.001
  test3 <- adonis2(genus_dist~Location, data = meta_dist, permutations=perm,
                   method = euclidean)
  write.csv(test2, file = paste(dataDir, " Genus HH location sig.csv", sep = ""))
  
  # Run pairwise adonis tests (change location/month to check for each)
  # pairwise_adonis_results <- pairwise.permmanova(genus_dist, wimp_wide_env$River_loc, nperm = 999, p.method = "bonferroni")
  # pairwise permmanova not working! can't find function
  # View the results
  # print(pairwise_adonis_results)

}

