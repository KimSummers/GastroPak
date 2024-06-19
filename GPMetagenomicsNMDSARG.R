# GPMetagenomicsNMDSAMR
#
# NDMS analysis of AMR data for GastroPak - code adapted from code supplied by Emma Travis
#
# file GPMetagenomicsNMDSAMR
#
# inputs
# 	metagenomicDataFile   - file containing plate count data
# 	metaDataFile          - file containing metadata
#   riversFile            - file to conver site into river name and stream position
# 	abClassesFile         - file containing antibiotic class data
#.  dataDir               - directory to save data in
# 	categories            - which categories e.g. seasons should data be split into
# 	sampleOrder           - Desired order for samples
#   plotsDir              - directory to store plots in

# 	normalisationFile     - file containing normalisation data
#   campaignType          - "Env" or "HH" for environmental or household
#   normType              - type of normalisation to use
#   multiplier            - value to multiple data by
#   average               - TRUE = average data across samples, FALSE do not
#   absLab                - Label for the y-axis for absolute data plots
#   yAxisLims             - Limits for the y-axis
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  26/03/24  Wellington Lab - School of Life Sciences - University of Warwick

GPMetagenomicsNMDSAMR <- function(metagenomicDataFile, metaDataFile, riversFile,
                                  distMatrixFile, categories, sampleOrder,
                                  plotsDir) {
  # Load necessary packages
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(pals)

  # Read in data
  mgAMRWide <- read.csv(metagenomicDataFile, sep = "\t")
  metaDataWide <- read.csv(metaDataFile)
  abClasses <- read.csv(abClassesFile, sep = "\t")

  abClasses$Drug_class[abClasses$Drug_class == "Ampenicol"] <- "Amphenicol"
  abClasses$Drug_class[abClasses$Drug_class == "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)"] <- "Rifamycin"
  abClasses$Drug_class[abClasses$Drug_class == "triclosan"] <- "Triclosan"

  # Define the desired order of levels and corresponding labels for treatment
  desiredOrder <- c("Dry", "Wet")
  newLabels <- c("Dry 2022", "Wet 2022")

  # Define the desired order of levels and corresponding labels for River_loc
  # desired_order1 <- c("Aire_up", "Aire_down", "Sowe_up", "Sowe_down")
  # new_labels1 <- c("River Aire upstream", "River Aire downstream", "River Sowe upstream", "River Sowe downstream")
  #

  colnames(metaDataWide)[which(colnames(metaDataWide) == "Sample.ID")] <- "Sample_ID"
  metaDataWide <- metaDataWide[metaDataWide$Sample_ID != "", ]
  riverSites <- unique(metaDataWide$Site_ID)
  riversData <- read_csv(file = riversFile)
  colnames(metaDataWide)[which(colnames(metaDataWide) == "Site_ID")] <- "SamplingSite"
  riversData$River[riversData$River == "Jinnah "] <- "Jinnah"
  riversData$River[riversData$River == "Korang "] <- "Korang"

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

  # preparing ARG data
  mgAMRWide[is.na(mgAMRWide)] <- 0
  mgAMRWide <- mgAMRWide %>% column_to_rownames("Annotation")

  mgAMRWide <- mgAMRWide[, colSums(mgAMRWide) > 0]
  mgAMRWide <- as.data.frame(t(mgAMRWide))
  rownames(mgAMRWide) <- sub("PAK0", "PAK_0", rownames(mgAMRWide))

  metaDataWide <- metaDataWide[metaDataWide$Barcode %in% rownames(mgAMRWide), ]

  # transform to relative abundance
  # hellinger might be better- square root of total and accounts for more high and low numbers
  amr.rel <- decostand(mgAMRWide, method = "total")

  # Calculate distance matrix - make distance matrix data frame - write .csv
  amrDistMat <- vegdist(amr.rel, method = "bray")
  amrDistanceMatrix <- as.matrix(amrDistMat, labels = T)
  write.csv(amrDistanceMatrix, paste(dataDir, "AMRDistMat.csv", sep = ""))

  # Running NMDS in vegan (metaMDS)
  amrNMS <- metaMDS(amrDistMat, distance = "bray", k = 2, maxit = 999,
                    trymax = 250, autotransform = FALSE)
  amrScores <- `sppscores<-`(amrNMS, mgAMRWide)
  amrCorrelation <- cor(mgAMRWide, amrNMS$points, use = "complete.obs",
                        method = "pearson")
  write.csv(amrCorrelation, file = paste(dataDir, "amr_PearsonCor.csv", sep = ""))

  # Shepards test/goodness of fit
  fitGoodness <- goodness(amrNMS)
  stressplot(amrNMS)

  pdf(file = paste(plotsDir, "AMRStressPlot.pdf", sep = ""))
  stressplot(amrNMS)
  dev.off()

  # Extract NMDS coordinates
  nmdsCoordinates <- scores(amrNMS, choices = c(1, 2))

  NMDSArg <- cbind(metaDataWide, nmdsCoordinates)

  argSppFit <- envfit(amrNMS, amr.rel, permutations = 999)

  # obtain all info from NMDS as dataframe
  amrScores <- as.data.frame(scores(amrScores, display = "species"))
  amrFit <- as.data.frame(scores(argSppFit, display = "vectors"))
  amrScores <- cbind(amrScores, Species = rownames(amrScores))

  # add species names to dataframe
  amrpVal <- as.data.frame(argSppFit$vectors$pvals)

  # add pvalues to dataframe so you can select species which are significant
  amrScores <- cbind(amrScores, pval = amrpVal$`argSppFit$vectors$pvals`)

  #abbreviate species names
  #spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6))

  #subset data to show species significant at 0.01
  sig.arg.scrs <- subset(amrScores, pval <= 0.01)

  #calculate vector (distance from 0,0)

  sig.arg.scrs$distance <- sqrt(sig.arg.scrs$NMDS1 ^ 2 + sig.arg.scrs$NMDS2 ^ 2)

  top20Sig.arg.scrs<- head(sig.arg.scrs[order(-sig.arg.scrs$distance), ], 20)

  #sigVars <- arg_spp_fit$signif
  graphColours <- c("deepskyblue", "deepskyblue4", "coral",
                    "coral4", "red", "blue", "chartreuse3",
                    "lightblue", "cornflowerblue","darkslategray",
                    "cadetblue4", "orange", "darkolivegreen", "grey",
                    "darkblue", "purple", "pink", "yellow", "cyan",
                    "aquamarine3")
  graphColours2 <- c("cyan2", "blue", "lightcyan3",
                    "lightpink2", "burlywood1", "darkorange3", "chartreuse3",
                    "lightblue", "cornflowerblue","darkslategray",
                    "cadetblue4", "orange", "darkolivegreen", "grey",
                    "darkblue", "purple", "pink", "yellow", "cyan",
                    "aquamarine3")

  #---- additional objects
  annColors <- list(
    River = c("Korang" = "lightpink2", "Nurpur" = "cyan2", "Rawal Lake" = "burlywood1",
              "Jinnah" = "lightcyan3", "Shahdara" = "blue"),
    `Sample Type` = c("River sediment" = "chocolate4" , "River water" = "skyblue"),
    group = c(Integrons = "#3288bd", Aminoglycoside = "#1a1a1a", Sulfonamide ="#fdae61",
              MGE = "#66c2a5",  "Beta-Lactam" = "#d53e4f", Other = "#e6f598" ,
              Phenicol ="#ffffbf", Quinolone = "#fee08b", Trimethoprim ="#8c510a",
              Tetracycline  ="#f46d43", "16S rRNA" = "deeppink", MDR = "grey49",
              MLSB = "#abdda4", Vancomycin = "#9e0142"),
    `River Location` = c(Upstream = "skyblue1",  Midstream = "slateblue4", Downstream = "violetred3"),
    Season = c(Wet = 'lightblue2', Dry = 'darkorange3'))

  riverColors <- annColors$River[match(argOrdMd$River, names(annColors$River))]
  seasonColours <- annColors$Season[match(argOrdMd$Season_WD, names(annColors$Season))]

  # nmds colored by river location and shape by season
  #sets up the plot
  nmds.plot.arg <- ggplot(NMDSArg, aes(x = NMDS1, y = NMDS2))
  #adds site points to plot, shape determined by month, colour determined by site
  nmds.plot.arg <- nmds.plot.arg + geom_point(aes(NMDS1, NMDS2, colour = factor(River),
                                                  shape = factor(CampaignLoc)), size = 4)
  nmds.plot.arg <- nmds.plot.arg + stat_ellipse(geom = "polygon", type = "t", alpha = 0.1,
                                                mapping = aes(fill = Campaign,
                                                              color = Campaign),
                                                level = 0.95)
  # add legend labels for Location and Month
  nmds.plot.arg <- nmds.plot.arg + labs(color = "River", shape = "River Site")
      # add legend labels for Location and Month
  nmds.plot.arg <- nmds.plot.arg + ggtitle("Antibiotic resistance gene NMDS - Grouping by season")
  nmds.plot.arg <- nmds.plot.arg + theme_bw() + guides(fill = FALSE)
  # Set colors
  # nmds.plot.arg <- nmds.plot.arg + scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral",
  #                                 "coral4", "red","blue", "chartreuse3",
  #                                 "lightblue","cornflowerblue","darkslategray",
  #                                 "cadetblue4","orange","darkolivegreen","grey",
  #                                 "darkblue"))
  nmds.plot.arg <- nmds.plot.arg + scale_color_manual(values = graphColours2)
  nmds.plot.arg <- nmds.plot.arg + scale_fill_manual(values = graphColours2[6:19])
  #  nmds.plot.arg <- nmds.plot.arg + scale_color_manual(values = glasbey())
  # Set shapes (19: circle, 17: triangle, 15: square)
  nmds.plot.arg <- nmds.plot.arg + scale_shape_manual(values = c(19, 17, 15, 1, 2, 0))
  nmds.plot.arg <- nmds.plot.arg + guides(color = guide_legend(order = 1),
                                          shape = guide_legend(order = 2))

  nmds.plot.arg

  ggsave(paste(plotsDir, "All ARG NMDS by season.pdf", sep = ""), nmds.plot.arg,
         height = 6, width = 8)

  nmds2.plot.arg <- ggplot(NMDSArg, aes(x = NMDS1, y = NMDS2))
  #adds site points to plot, shape determined by month, colour determined by site
  nmds2.plot.arg <- nmds2.plot.arg + geom_point(aes(NMDS1, NMDS2, colour = factor(River),
                                                    shape = factor(CampaignLoc)), size = 4)
  nmds2.plot.arg <- nmds2.plot.arg + stat_ellipse(geom = "polygon", type = "t", alpha = 0.1,
                                                  mapping = aes(fill = River, color = River))
  # add legend labels for Location and Month
  nmds2.plot.arg <- nmds2.plot.arg + labs(color = "location", shape = "Campaign")
  # add legend labels for Location and Month
  nmds2.plot.arg <- nmds2.plot.arg + ggtitle("Antibiotic resistance gene NMDS - Grouping by location")
  nmds2.plot.arg <- nmds2.plot.arg + theme_bw() + guides(fill = FALSE)
  # Set colors
  # nmds2.plot.arg <- nmds2.plot.arg + scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral",
  #                                   "coral4", "red","blue", "chartreuse3",
  #                                   "lightblue","cornflowerblue","darkslategray",
  #                                   "cadetblue4","orange","darkolivegreen","grey",
  #                                   "darkblue"))
  nmds2.plot.arg <- nmds2.plot.arg + scale_color_manual(values = graphColours2)
  nmds2.plot.arg <- nmds2.plot.arg + scale_fill_manual(values = graphColours2)
  # Set shapes (19: circle, 17: triangle, 15: square)
  nmds2.plot.arg <- nmds2.plot.arg + scale_shape_manual(values = c(19, 17, 15, 1, 2, 0))
  nmds2.plot.arg <- nmds2.plot.arg + guides(color = guide_legend(order = 1),
                                            shape = guide_legend(order = 2))

  nmds2.plot.arg

  ggsave(paste(plotsDir, "All ARG NMDS by river.pdf", sep = ""), nmds2.plot.arg,
         height = 6, width = 8)

  #ADD significant species (i.e. arg driving diversity)
  nmds.plot.arg3 <- nmds.plot.arg +
    geom_segment(data = sig.arg.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)
  #add vector arrows of significant species
  nmds.plot.arg3 <- nmds.plot.arg3 +
    ggrepel::geom_text_repel(data = sig.arg.scrs, aes(x = NMDS1, y = NMDS2,
                                                      label = Species),
                             cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.arg3
  #too many arg driving diversity so how to reduce this? Reduced sig to <0.01 not <0.05 - maybe also use size of vector??! / combine CTX-M etc?

  ggsave(paste(plotsDir, "All ARG NMDS plus vectors.pdf", sep = ""), nmds.plot.arg3,
         height = 6, width = 8)

  #ADD significant species (i.e. arg driving diversity)
  nmds.plot.arg4 <- nmds.plot.arg +
    geom_segment(data = top20Sig.arg.scrs, aes(x = 0, xend = NMDS1, y = 0,
                                                yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)

  #add vector arrows of significant species
  nmds.plot.arg4 <- nmds.plot.arg4 +
    ggrepel::geom_text_repel(data = top20Sig.arg.scrs,
                             aes(x = NMDS1, y = NMDS2, label = Species), cex = 3,
                             direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.arg4

  ggsave(paste(plotsDir, "All ARG NMDS plus top 20 vectors.pdf", sep = ""), nmds.plot.arg4,
         height = 6, width = 8)

  # link the name of arg to drug_class etc

  # Create new column 'name' with the same data as 'Species'
  sig.arg.scrs$Subject <- sig.arg.scrs$Species

  # add metadata information to the amr_count file
  sig.arg.scrs.meta <- inner_join(sig.arg.scrs, abClasses, by = "Subject")

  #ADD significant species (i.e. arg driving diversity)
  nmds.plot.arg5 <- nmds.plot.arg +
    geom_segment(data = sig.arg.scrs.meta,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, colour = Drug_class),
                 arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3)
  # add vector arrows of significant species
  nmds.plot.arg5 <- nmds.plot.arg5 +
    ggrepel::geom_text_repel(data = sig.arg.scrs.meta,
                             aes(x = NMDS1, y = NMDS2, label = Species, colour = Drug_class),
                             cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.arg5

  ggsave(paste(plotsDir, "All ARG plus ARG classes.pdf", sep = ""),
         nmds.plot.arg5, height = 6, width = 8)

  #ADD significant species (i.e. arg driving diversity)
  nmds.plot.arg6 <- nmds2.plot.arg +
    geom_segment(data = sig.arg.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)
  #add vector arrows of significant species
  nmds.plot.arg6 <- nmds.plot.arg6 +
    ggrepel::geom_text_repel(data = sig.arg.scrs, aes(x = NMDS1, y = NMDS2,
                                                      label = Species),
                             cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.arg6
  #too many arg driving diversity so how to reduce this? Reduced sig to <0.01 not <0.05 - maybe also use size of vector??! / combine CTX-M etc?

  ggsave(paste(plotsDir, "All ARG NMDS by river plus vectors.pdf", sep = ""),
         nmds.plot.arg6, height = 6, width = 8)

  #ADD significant species (i.e. arg driving diversity)
  nmds.plot.arg7 <- nmds2.plot.arg +
    geom_segment(data = top20Sig.arg.scrs, aes(x = 0, xend = NMDS1, y = 0,
                                               yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd = 0.3)

  #add vector arrows of significant species
  nmds.plot.arg7 <- nmds.plot.arg7 +
    ggrepel::geom_text_repel(data = top20Sig.arg.scrs,
                             aes(x = NMDS1, y = NMDS2, label = Species), cex = 3,
                             direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.arg7

  ggsave(paste(plotsDir, "All ARG NMDS by river plus top 20 vectors.pdf", sep = ""),
         nmds.plot.arg7, height = 6, width = 8)

  #ADD significant species (i.e. arg driving diversity)
  nmds.plot.arg8 <- nmds2.plot.arg +
    geom_segment(data = sig.arg.scrs.meta,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, colour = Drug_class),
                 arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3)
  # add vector arrows of significant species
  nmds.plot.arg8 <- nmds.plot.arg8 +
    ggrepel::geom_text_repel(data = sig.arg.scrs.meta,
                             aes(x = NMDS1, y = NMDS2, label = Species, colour = Drug_class),
                             cex = 3, direction = "both", segment.size = 0.25)
  # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  nmds.plot.arg8

  ggsave(paste(plotsDir, "All ARG NMDS by river plus AMR classes.pdf", sep = ""),
         nmds.plot.arg8, height = 6, width = 8)

  # # Filter out a particular class, for example, class 'A'
  # filt_sig.arg.scrs.meta <- subset(sig.arg.scrs.meta, drug_class != "MDR efflux pump")
  #
  # nmds.plot.arg9 <- nmds.plot.arg +
  #   geom_segment(data = filt_sig.arg.scrs.meta,
  #                aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, colour = Drug_class),
  #                arrow = arrow(length = unit(0.25, "cm")),
  #                lwd = 0.3)
  # # add vector arrows of significant species
  # nmds.plot.arg9 <- nmds.plot.arg9 +
  #   ggrepel::geom_text_repel(data = filt_sig.arg.scrs.meta,
  #                            aes(x = NMDS1, y = NMDS2, label = Species, colour = Drug_class),
  #                            cex = 3, direction = "both", segment.size = 0.25)
  # # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  #
  # nmds.plot.arg9
  #
  ## -------Perform PERMANOVA to check significance of clustering from ordination
  # Produce distance matrix
  # Load data
  # Assuming row names are sample IDs in both data.matrix and sample_metadata

  # Calculating distance matrix using vegan function vegdist
  bray_dist_mat_arg <- vegdist(mgAMRWide, method = "bray")
  amr_dist <- as.data.frame(amrDistanceMatrix)
  amr_dist2 <- cbind(barcode = rownames(amr_dist), amr_dist)
  colnames(metaDataWide)[which(colnames(metaDataWide) == "Barcode")] <- "barcode"

  meta_dist <- inner_join(metaDataWide, amr_dist2, by = "barcode")

  #permutations for tests
  perm <- 999
  test1 <- adonis2(amr_dist~River, data = meta_dist, permutations = perm,
                   method = euclidean)
  # sig p = 0.001
  test2 <- adonis2(amr_dist~River*Campaign, data = meta_dist, permutations = perm,
                   method=euclidean)
  #River p = 0.001, Campaign p = 0.001, River * Campaign p = 0.001
  test3 <- adonis2(amr_dist~Campaign, data = meta_dist, permutations=perm,
                   method = euclidean)
  # sig p=0.001

  test8 <- adonis2(amr_dist~RiverLocation*Campaign*River, data = meta_dist, permutations=perm,
                   method = euclidean)

  # Run pairwise adonis tests (change location/month to check for each)
  # pairwise_adonis_results <- pairwise.perm.manova(amr_dist, wimp_wide_env$River_loc, nperm = 999, p.method = "bonferroni")
  # pairwise perm manova not working! can't find function
  # View the results
  # print(pairwise_adonis_results)

}

