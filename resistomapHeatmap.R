# resistomapHeatmap
#
# Heatmap ResistoMap data for GastroPak - code adapted from code supplied by Emma Travis
#
# file resistomapHeatmap
#
# inputs
# 	dataFile      - file containing plate count data
#   metaDataFile  - file containing metadata on samples
#   sampleNumber  - Number of samples
#   riversFile    - file with metadata of rivers
#   plotsDir      - directory to store plots in
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  31/05/24  Wellington Lab - School of Life Sciences - University of Warwick

resistomapHeatmap <- function(dataFile, metaDataFile, sampleNumber, plotsDir) {

  #---- RESISTOMAP Dry summer (May 2022)----#

  # libraries <-c("vegan", "pheatmap", "ggplot2", "dplyr", "sf", "tidyverse",
  #               "ggmap", "modeest", "ggpubr", "tidyheatmaps", "maps",
  #               "ggrepel", "cowplot", "gridExtra", "readxl", "pals", "dunn.test")
  # lapply(libraries, require, character.only = TRUE)
  #
  library(tidyverse)
  library(ggpubr)
  library(ggplot2)
  library(tidyheatmaps)
  library(vegan)
  library(ggmap)

  metaData <- read.delim(metaDataFile, sep = ",")
  metaData <- metaData[1:sampleNumber, ]
  colnames(metaData)[which(colnames(metaData) == "Site_ID")] <- "SamplingSite"

  riversData <- read_csv(file = riversFile)
  riversData$River[7:9] <- "Jinnah"
  riversData$River[10:12] <- "Korang"
  
  metaData <- addRivers(countData = metaData, metaData = riversData)
  
  colnames(metaData)[which(colnames(metaData) == "SamplingSite")] <- "Site_ID"
  metaData$Sample_type[metaData$Sample_type == "River_water"] <- "River water"
  metaData$Sample_type[metaData$Sample_type == "River_sediment"] <- "River sediment"
  colnames(metaData)[which(colnames(metaData) == "RiverLocation")] <- "River Location"
  colnames(metaData)[which(colnames(metaData) == "Sample_type")] <- "Sample Type"

  rmData <- read.csv(rmFile)
  
  # transform rmData into long format 
  rmDataLong <- pivot_longer(rmData, cols = starts_with("PAK"), names_to = "Sample_ID",
                             values_to = "relative_abundance")

  #statistical analysis
  # test for normality
  ggdensity(rmDataLong$relative_abundance, fill = "lightgray")
  ggqqplot(rmDataLong$relative_abundance)
  shapiro.test(rmDataLong$relative_abundance)

  # Shapiro-Wilk normality test
  # data:  rmDataLong$relative_abundance

  # use log10 and NA instead of 0 values and test again for normality
  rmDataLongLog <- mutate(rmDataLong, Log_RA = log10(rmDataLong$relative_abundance))

  ggdensity(rmDataLongLog$Log_RA, fill = "lightgray")
  ggqqplot(rmDataLongLog$Log_RA)
  shapiro.test(rmDataLongLog$Log_RA)

  # Shapiro-Wilk normality test
  # data:  rmDataLongLog$Log_RA

  colnames(metaData)[which(colnames(metaData) == "Sample_ID")] <- "Location_ID"
  colnames(metaData)[which(colnames(metaData) == "Barcode")] <- "Sample_ID"

  # add sample metadata using sample ID to merge files
  rmDataLongLog <- full_join(rmDataLongLog, metaData, by = "Sample_ID")

  # boxplot all targets
  plot <- ggboxplot(rmDataLongLog, x = "Sample_ID", y = "Log_RA", color = "Sample Type")
  plot <- facet(plot, facet.by = "group", ncol = 3) + rotate_x_text(angle = 90)
  plot <- plot + theme(axis.text = element_text(size = 8))
  plot <- plot + theme(axis.title = element_text(size = 10))
  plot

  ggsave(paste(plotsDir, "Log RA sample type comp.pdf", sep = ""), plot, height = 20,
         width = 20)

  plot <- ggboxplot(rmDataLongLog, x = "Sample Type", y = "Log_RA", color = "Site_ID")
  facet(plot, facet.by = "group", ncol = 3) + rotate_x_text(angle = 90)

  ggsave(paste(plotsDir, "Log RA sample ID comp.pdf", sep = ""), plot, height = 20,
         width = 20)

  plot <- ggboxplot(rmDataLongLog, x = "Sample Type", y = "Log_RA", color = "Location_ID")
  facet(plot, facet.by = "group", ncol = 3) + rotate_x_text(angle = 90)

  ggsave(paste(plotsDir, "Log RA location ID comp.pdf", sep = ""), plot, height = 20,
         width = 20)
  
  plot <- ggboxplot(rmDataLongLog, x = "Site_ID", y = "Log_RA", color = "Sample Type")
  facet(plot, facet.by = "group", ncol = 3) + rotate_x_text(angle = 90)

  ggsave(paste(plotsDir, "Log RA sample ID and type comp.pdf", sep = ""), plot, height = 20,
         width = 20)
  
  #remove genes that had no detection across all samples or were only detected in 1-2 samples
  filterTab <- rmData[rowSums(rmData[, 4:ncol(rmData)], na.rm = TRUE) > 0, ]

  filterTabLong <- pivot_longer(filterTab, cols = starts_with("PAK"),
                                names_to = "Sample_ID",
                                values_to = "relative_abundance")
  filterTabLongLog <- mutate(filterTabLong, Log_RA = log10(filterTabLong$relative_abundance))
  filterTabLongLog <- full_join(filterTabLongLog, metaData, by = "Sample_ID")

  # divide water and sediment
  waterfiltLongLog <- filter(filterTabLongLog, `Sample Type` != "River sediment")
  sedfiltLongLog <- filter(filterTabLongLog, `Sample Type` == "River sediment")

  aggregateGenes <- aggregate(waterfiltLongLog$relative_abundance,
                              by = list(waterfiltLongLog$gene), FUN = sum, na.rm = TRUE)
  waterfiltLongLog <- filter(waterfiltLongLog, gene %in% aggregateGenes$Group.1[aggregateGenes$x > 0])

  aggregateGenes <- aggregate(sedfiltLongLog$relative_abundance,
                              by = list(sedfiltLongLog$gene), FUN = sum, na.rm = TRUE)
  sedfiltLongLog <- filter(sedfiltLongLog, gene %in% aggregateGenes$Group.1[aggregateGenes$x > 0])

  # additional objects
  annColors <- list(
    River = c(Jinnah = "lightcyan3", Korang = "lightpink2", Nurpur = "cyan2",
              "Rawal Lake" = "burlywood1", Shahdara = "blue"),
    `Sample Type` = c("River sediment" = "chocolate4" , "River water" = "skyblue"),
    group = c(Integrons = "#3288bd", Aminoglycoside = "#1a1a1a", Sulfonamide ="#fdae61",
              MGE = "#66c2a5",  "Beta-Lactam" = "#d53e4f", Other = "#e6f598" ,
             Phenicol ="#ffffbf", Quinolone = "#fee08b", Trimethoprim ="#8c510a",
             Tetracycline  ="#f46d43", "16S rRNA" = "deeppink", MDR = "grey49",
             MLSB = "#abdda4", Vancomycin = "#9e0142"),
    `River Location` = c(Upstream = "skyblue1",  Midstream = "slateblue4", Downstream = "violetred3"))

  colorRampPalette(c('#2166ac', '#67a9cf','#d1e5f0','#fddbc7', '#ef8a62','#b2182b'))
  colorRampPalette(c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac'))

  filterTabLongLog %>%
    arrange(group) %>%
    tidyheatmap(
      filterTabLongLog,
      rows = gene,
      columns =  Sample_ID,
      values =  Log_RA,
      #scale = "row",
      color_legend_n = 6,
      color_legend_min = -5,
      color_legend_max = 1,
      annotation_col = c(River, `Sample Type`, `River Location`),
      annotation_row = c(group),
      annotation_colors = annColors,
      cluster_rows = clustRow,
      cluster_cols = clustCol,
      gaps_row = group,
      filename = paste(plotsDir, "Heat map all cluster row ", clustRow, " col ",
                       clustCol, ".pdf", sep = ""),
      show_colnames = FALSE
    )

    tidyheatmap(
      hmDataLongLog,
      rows = ID,
      columns =  Sample_ID,
      values =  Log_RA,
      #scale = "row",
      color_legend_n = 6,
      color_legend_min = -5,
      color_legend_max = 1,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      filename = "/Users/u1971521/Dropbox/Alberti/DataAnalysis/TestPlot6.pdf",
      show_colnames = FALSE
    )

    waterfiltLongLog %>%
    arrange(group) %>%
    tidyheatmap(
      waterfiltLongLog,
      rows = gene,
      columns =  Sample_ID,
      values =  Log_RA,
      #scale = "row",
      color_legend_n = 6,
      color_legend_min = -5,
      color_legend_max = 1,
      annotation_col = c(`River Location`, River),
      annotation_row = c(group),
      annotation_colors = annColors,
      cluster_rows = clustRow,
      cluster_cols = clustCol,
      gaps_row = group,
      filename = paste(plotsDir, "Heat map water cluster row ", clustRow, " col ",
                       clustCol, ".pdf", sep = ""),
      show_colnames = FALSE
    )

  sedfiltLongLog %>%
    arrange(group) %>%
    tidyheatmap(
      sedfiltLongLog,
      rows = gene,
      columns =  Sample_ID,
      values =  Log_RA,
      #scale = "row",
      color_legend_n = 6,
      color_legend_min = -5,
      color_legend_max = 1,
      annotation_col = c(`River Location`, River),
      annotation_row = c(group),
      annotation_colors = annColors,
      cluster_rows = clustRow,
      cluster_cols = clustCol,
      gaps_row = group,
      filename = paste(plotsDir, "Heat map sediment cluster row ", clustRow, " col ",
                       clustCol, ".pdf", sep = ""),
      show_colnames = FALSE
    )

  #  NMDS and PCoA / alpha diversity ----
  # CHECK THIS LINK https://www.rpubs.com/RGrieger/545184 for further details of these steps

# Water only start

#remove taxonomic data
waterfiltLongLogNoTx <- filter(waterfiltLongLog, group != "Taxonomic")

# transform dataframe with sample in row, arg RA in column and NA into 0
waterfiltLongLogNoTxWide <- pivot_wider(waterfiltLongLogNoTx, id_cols = Sample_ID,
                                        names_from = gene, values_from = relative_abundance) %>% replace(is.na(.), 0)
waterfiltLongLogNoTxWide <- as.data.frame(waterfiltLongLogNoTxWide)
rownames(waterfiltLongLogNoTxWide) <- waterfiltLongLogNoTxWide$Sample_ID
waterfiltLongLogNoTxWide <- waterfiltLongLogNoTxWide[, -1]


# data ordination using bray curtis
argOrd <- metaMDS(waterfiltLongLogNoTxWide, distance = "bray", autotransform = F, k = 3)
argOrd # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
argSppFit <- envfit(argOrd, waterfiltLongLogNoTxWide, permutations = 999)

# to use ggplot the results form the ordination need to be saved into a separate dataframe.
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. this can be done this by calling the scores of the mds.
argOrdMd <- merge(argOrd$points, metaData, by.x = "row.names", by.y = "Sample_ID", all.x = TRUE)

# A new dataset containing species data also needs to be made to look at species vectors.
spp.scrs <- as.data.frame(scores(argSppFit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = argSppFit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
# spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

head(spp.scrs)

# To show environmental extrinsic variables another datasheet needs to be created. CHECK LINK https://www.rpubs.com/RGrieger/545184

# Now you have the relevant information for plotting the ordination in ggplot! Lets get plotting!
riverColors <- annColors$River[match(argOrdMd$River, names(annColors$River))]

#sets up the plot
nmds.plot.arg <- ggplot(argOrdMd, aes(x = MDS1, y = MDS2))
#adds site points to plot, shape determined by Landuse, colour determined by Management
nmds.plot.arg <- nmds.plot.arg + geom_point(aes(MDS1, MDS2, colour = factor(River),
                                                shape = factor(`River Location`)),
                                            size = 4, fill = riverColors)

#add or exclude this command if you want aspect ratio fixed or variable
# nmds.plot.arg <- nmds.plot.arg + coord_fixed()

# Use defined colors for River
# Add ellipses
nmds.plot.arg <- nmds.plot.arg + scale_color_manual(values = riverColors)
nmds.plot.arg <- nmds.plot.arg + stat_ellipse(aes(colour = factor(River)),
                                              level = 0.95, type = "norm")
nmds.plot.arg <- nmds.plot.arg + theme_light()
  #theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
# add legend labels for River and River location
nmds.plot.arg <- nmds.plot.arg + labs(colour = "River", shape = "River location")
# add legend at right of plot
nmds.plot.arg <- nmds.plot.arg + theme(legend.position = "right",
                                       legend.text = element_text(size = 18),
                                       legend.title = element_text(size = 21),
                                       axis.text = element_text(size = 18))
nmds.plot.arg <- nmds.plot.arg + theme(text = element_text(size = 21))

# labs(title = "Ordination with arg vectors")

nmds.plot.arg
# ADD significant species (i.e. arg driving diversity) with the following code
# add vector arrows of significant species
nmds.plot.arg <- nmds.plot.arg + geom_segment(data = sig.spp.scrs,
                                              aes(x = 0, xend = NMDS1, y = 0,
                                                  yend = NMDS2),
                                              arrow = arrow(length = unit(0.25, "cm")),
                                              colour = "grey10", lwd = 0.3)
nmds.plot.arg <- nmds.plot.arg + ggrepel::geom_text_repel(data = sig.spp.scrs,
                                                          aes(x = NMDS1, y = NMDS2,
                                                              label = Species), cex = 3,
                                                          direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
nmds.plot.arg

nmds.plot.arg2 <- ggplot(argOrdMd, aes(x = MDS2, y = MDS3)) + #sets up the plot
  geom_point(aes(MDS2, MDS3, colour = factor(River), shape = factor(`River Location`)),
             size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", linewidth = 1,
                                        linetype = "solid")) +
  labs(colour = "Site_ID", shape = "River Location") + # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size = 21))
nmds.plot.arg2


nmds.plot.arg3 <- ggplot(argOrdMd, aes(x = MDS1, y = MDS3)) + #sets up the plot
  geom_point(aes(MDS1, MDS3, colour = factor(River), shape = factor(`River Location`)),
             size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1,
                                        linetype = "solid"))+
  labs(colour = "River", shape = "River Location")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text = element_text(size=21))
nmds.plot.arg3

combined_plots <- ggarrange(nmds.plot.arg, nmds.plot.arg2, nmds.plot.arg3, ncol = 3)
ggsave(paste(plotsDir, "combined_plots.pdf", sep = ""), combined_plots,
       width = 30, height = 9)

# Open a PDF file to save the plots
pdf(paste(plotsDir, "combined_plots.pdf", sep = ""))

# Arrange plots in a 3x1 grid
par(mfrow = c(3, 1))

# Print the last three plots
print(nmds.plot.arg)
print(nmds.plot.arg2)
print(nmds.plot.arg3)

# Close the PDF file
dev.off()

# water only end

# All sample types start

#remove taxonomic data
filtLongLogNoTx <- filter(filterTabLongLog, group != "Taxonomic")

# transform dataframe with sample in row, arg RA in column and NA into 0
filtLongLogNoTxWide <- pivot_wider(filterTabLongLog, id_cols = Sample_ID,
                                   names_from = gene, values_from = relative_abundance) %>% replace(is.na(.), 0)
filtLongLogNoTxWide <- as.data.frame(filtLongLogNoTxWide)
rownames(filtLongLogNoTxWide) <- filtLongLogNoTxWide$Sample_ID
filtLongLogNoTxWide <- filtLongLogNoTxWide[, -1]


# data ordination using bray curtis
argOrdAll <- metaMDS(filtLongLogNoTxWide, distance = "bray", autotransform = F, k = 3)
argOrdAll # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
argSppFit.all <- envfit(argOrdAll, filtLongLogNoTxWide, permutations = 999)

# to use ggplot the results form the ordination need to be saved into a separate dataframe.
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. this can be done this by calling the scores of the mds.
argOrdMd.all <- merge(argOrdAll$points, metaData, by.x = "row.names",
                      by.y = "Sample_ID", all.x = TRUE)

# A new dataset containing species data also needs to be made to look at species vectors.
spp.scrs.all <- as.data.frame(scores(argSppFit.all, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs.all <- cbind(spp.scrs.all, Species = rownames(spp.scrs.all)) #add species names to dataframe
spp.scrs.all <- cbind(spp.scrs.all, pval = argSppFit.all$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
# spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs.all <- subset(spp.scrs.all, pval <= 0.05) #subset data to show species significant at 0.05

head(spp.scrs.all)

# To show environmental extrinsic variables another datasheet needs to be created. CHECK LINK https://www.rpubs.com/RGrieger/545184

# Now you have the relevant information for plotting the ordination in ggplot! Lets get plotting!
riverColors <- annColors$River[match(argOrdMd.all$River, names(annColors$River))]

#sets up the plot
nmds.plot.arg.all <- ggplot(argOrdMd.all, aes(x = MDS1, y = MDS2))
#adds site points to plot, shape determined by Landuse, colour determined by Management
nmds.plot.arg.all <- nmds.plot.arg.all + geom_point(aes(MDS1, MDS2, colour = factor(River),
                                                        shape = factor(`River Location`)),
                                                    size = 4, fill = riverColors)

#add or exclude this command if you want aspect ratio fixed or variable
# nmds.plot.arg <- nmds.plot.arg + coord_fixed()

# Use defined colors for River
# Add ellipses
nmds.plot.arg.all <- nmds.plot.arg.all + scale_color_manual(values = riverColors)
nmds.plot.arg.all <- nmds.plot.arg.all + stat_ellipse(aes(colour = factor(River)),
                                              level = 0.95, type = "norm")
nmds.plot.arg.all <- nmds.plot.arg.all + theme_light()
#theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
# add legend labels for River and River location
nmds.plot.arg.all <- nmds.plot.arg.all + labs(colour = "River", shape = "River location")
# add legend at right of plot
nmds.plot.arg.all <- nmds.plot.arg.all + theme(legend.position = "right",
                                               legend.text = element_text(size = 18),
                                               legend.title = element_text(size = 21),
                                       axis.text = element_text(size = 18))
nmds.plot.arg.all <- nmds.plot.arg.all + theme(text = element_text(size = 21))

# labs(title = "Ordination with arg vectors")

nmds.plot.arg.all
# ADD significant species (i.e. arg driving diversity) with the following code
# add vector arrows of significant species
nmds.plot.arg.all <- nmds.plot.arg.all + geom_segment(data = sig.spp.scrs.all,
                                              aes(x = 0, xend = NMDS1, y = 0,
                                                  yend = NMDS2),
                                              arrow = arrow(length = unit(0.25, "cm")),
                                              colour = "grey10", lwd = 0.3)
nmds.plot.arg.all <- nmds.plot.arg.all + ggrepel::geom_text_repel(data = sig.spp.scrs.all,
                                                          aes(x = NMDS1, y = NMDS2,
                                                              label = Species), cex = 3,
                                                          direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
nmds.plot.arg.all

nmds.plot.arg2.all <- ggplot(argOrdMd.all, aes(x = MDS2, y = MDS3)) + #sets up the plot
  geom_point(aes(MDS2, MDS3, colour = factor(River), shape = factor(`River Location`)),
             size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", linewidth = 1,
                                        linetype = "solid")) +
  labs(colour = "Site_ID", shape = "River Location") + # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size = 21))
nmds.plot.arg2.all


nmds.plot.arg3.all <- ggplot(argOrdMd.all, aes(x = MDS1, y = MDS3)) + #sets up the plot
  geom_point(aes(MDS1, MDS3, colour = factor(River), shape = factor(`River Location`)),
             size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1,
                                        linetype = "solid"))+
  labs(colour = "River", shape = "River Location")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text = element_text(size = 21))
nmds.plot.arg3.all

combined_plots.all <- ggarrange(nmds.plot.arg.all, nmds.plot.arg2.all, nmds.plot.arg3.all, ncol = 3)
ggsave(paste(plotsDir, "combined_plots all data.pdf", sep = ""), combined_plots.all,
       width = 30, height = 9)

# Open a PDF file to save the plots
pdf(paste(plotsDir, "combined_plots all.pdf", sep = ""))

# Arrange plots in a 3x1 grid
par(mfrow = c(3, 1))

# Print the last three plots
print(nmds.plot.arg.all)
print(nmds.plot.arg2.all)
print(nmds.plot.arg3.all)

# Close the PDF file
dev.off()

# all sample types end

# Sediment only start

#remove taxonomic data
sedfiltLongLogNoTx <- filter(sedfiltLongLog, group != "Taxonomic")

# transform dataframe with sample in row, arg RA in column and NA into 0
sedfiltLongLogNoTxWide <- pivot_wider(sedfiltLongLogNoTx, id_cols = Sample_ID,
                                      names_from = gene, values_from = relative_abundance) %>% replace(is.na(.), 0)
sedfiltLongLogNoTxWide <- as.data.frame(sedfiltLongLogNoTxWide)
rownames(sedfiltLongLogNoTxWide) <- sedfiltLongLogNoTxWide$Sample_ID
sedfiltLongLogNoTxWide <- sedfiltLongLogNoTxWide[, -1]


# data ordination using bray Curtis
argOrd.sed <- metaMDS(sedfiltLongLogNoTxWide, distance = "bray", autotransform = F, k = 3)
argOrd.sed # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
argSppFit.sed <- envfit(argOrd.sed, sedfiltLongLogNoTxWide, permutations = 999)

# to use ggplot the results form the ordination need to be saved into a separate dataframe.
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. this can be done this by calling the scores of the mds.
argOrdMd.sed <- merge(argOrd.sed$points, metaData, by.x = "row.names", by.y = "Sample_ID", all.x = TRUE)

# A new dataset containing species data also needs to be made to look at species vectors.
spp.scrs.sed <- as.data.frame(scores(argSppFit.sed, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs.sed <- cbind(spp.scrs.sed, Species = rownames(spp.scrs.sed)) #add species names to dataframe
spp.scrs.sed <- cbind(spp.scrs.sed, pval = argSppFit.sed$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
# spp.scrs.sed<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs.sed <- subset(spp.scrs.sed, pval <= 0.05) #subset data to show species significant at 0.05

head(spp.scrs.sed)

# To show environmental extrinsic variables another datasheet needs to be created. CHECK LINK https://www.rpubs.com/RGrieger/545184

# Now you have the relevant information for plotting the ordination in ggplot! Lets get plotting!
riverColors <- annColors$River[match(argOrdMd.sed$River, names(annColors$River))]

#sets up the plot
nmds.plot.arg.sed <- ggplot(argOrdMd.sed, aes(x = MDS1, y = MDS2))
#adds site points to plot, shape determined by Landuse, colour determined by Management
nmds.plot.arg.sed <- nmds.plot.arg.sed + geom_point(aes(MDS1, MDS2, colour = factor(River),
                                                shape = factor(`River Location`)),
                                            size = 4, fill = riverColors)

#add or exclude this command if you want aspect ratio fixed or variable
# nmds.plot.arg.sed <- nmds.plot.arg.sed + coord_fixed()

# Use defined colors for River
# Add ellipses
nmds.plot.arg.sed <- nmds.plot.arg.sed + scale_color_manual(values = riverColors)
nmds.plot.arg.sed <- nmds.plot.arg.sed + stat_ellipse(aes(colour = factor(River)),
                                              level = 0.95, type = "norm")
nmds.plot.arg.sed <- nmds.plot.arg.sed + theme_light()
#theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
# add legend labels for River and River location
nmds.plot.arg.sed <- nmds.plot.arg.sed + labs(colour = "River", shape = "River location")
# add legend at right of plot
nmds.plot.arg.sed <- nmds.plot.arg.sed + theme(legend.position = "right",
                                       legend.text = element_text(size = 18),
                                       legend.title = element_text(size = 21),
                                       axis.text = element_text(size = 18))
nmds.plot.arg.sed <- nmds.plot.arg.sed + theme(text = element_text(size = 21))

# labs(title = "Ordination with arg vectors")

nmds.plot.arg.sed
# ADD significant species (i.e. arg driving diversity) with the following code
# add vector arrows of significant species
nmds.plot.arg.sed <- nmds.plot.arg.sed + geom_segment(data = sig.spp.scrs.sed,
                                              aes(x = 0, xend = NMDS1, y = 0,
                                                  yend = NMDS2),
                                              arrow = arrow(length = unit(0.25, "cm")),
                                              colour = "grey10", lwd = 0.3)
nmds.plot.arg.sed <- nmds.plot.arg.sed + ggrepel::geom_text_repel(data = sig.spp.scrs.sed,
                                                          aes(x = NMDS1, y = NMDS2,
                                                              label = Species), cex = 3,
                                                          direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
nmds.plot.arg.sed

nmds.plot.arg2.sed <- ggplot(argOrdMd.sed, aes(x = MDS2, y = MDS3)) + #sets up the plot
  geom_point(aes(MDS2, MDS3, colour = factor(River), shape = factor(`River Location`)),
             size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, colour = "black", linewidth = 1,
                                        linetype = "solid")) +
  labs(colour = "Site_ID", shape = "River Location") + # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size = 21))
nmds.plot.arg2.sed


nmds.plot.arg3.sed <- ggplot(argOrdMd.sed, aes(x = MDS1, y = MDS3)) + #sets up the plot
  geom_point(aes(MDS1, MDS3, colour = factor(River), shape = factor(`River Location`)),
             size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1,
                                        linetype = "solid"))+
  labs(colour = "River", shape = "River Location")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text = element_text(size = 21))
nmds.plot.arg3.sed

combined_plots.sed <- ggarrange(nmds.plot.arg.sed, nmds.plot.arg2.sed, nmds.plot.arg3.sed, ncol = 3)
ggsave(paste(plotsDir, "combined_plots sediment data.pdf", sep = ""), combined_plots.sed,
       width = 30, height = 9)

# Open a PDF file to save the plots
pdf(paste(plotsDir, "combined_plots sediment.pdf", sep = ""))

# Arrange plots in a 3x1 grid
par(mfrow = c(3, 1))

# Print the last three plots
print(nmds.plot.arg.sed)
print(nmds.plot.arg2.sed)
print(nmds.plot.arg3.sed)

# Close the PDF file
dev.off()

# sediment only end

# function for ellipsess

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the River (location) factor
df_ell.arg.location <- data.frame() #sets up a data frame before running the function.

for(g in unique(site.scrs$River)){
  df_ell.arg.location <-
    rbind(df_ell.arg.location,
          cbind(as.data.frame(with(site.scrs [site.scrs$Site_ID == g, ],
                                   veganCovEllipse(cov.wt(cbind(NMDS1, NMDS2),
                                                          wt = rep(1 / length(NMDS1),
                                                                   length(NMDS1)))$cov,
                                                   center=c(mean(NMDS1),
                                                            mean(NMDS2))))), Site_ID = g))
}

# data for labeling the ellipse
NMDS.mean.arg = aggregate(site.scrs[ ,c("NMDS1", "NMDS2")],
                          list(group = site.scrs$Site_ID), mean)

# data for labeling the ellipse
NMDS.mean=aggregate(site.scrs[, c("NMDS1", "NMDS2")],
                    list(group = site.scrs$Site_ID), mean)
nmds.plot.arg +
  geom_path(data = df_ell.arg.location, aes(x = NMDS1, y = NMDS2, group = River,
                                            colour = River)) + #this is the ellipse, seperate ones by Site.
  theme(text = element_text(size = 21))

ggsave(
  filename = "nmds_arg_site.jpeg",
  plot = last_plot(),
  device = "jpeg",
  scale = 1,
  dpi = 300,
)



#-----------------------anosim args ------------------
#allCW_RA_6_arg_0 <- allCW_RA_6_arg
#allCW_RA_6_arg_0[is.na(allCW_RA_6_arg_0)] <- 0


#anosimARG = anosim(allCW_RA_6_arg_0, allCW_RA_8_arg$group, distance = "bray", permutations = 9999)
#anosimARG

#anosimARG = anosim(allCW_RA_6_arg_0, allCW_RA_8_argT, distance = "bray", permutations = 9999)
#anosimARG

#-----------------------alpha diveristy of resistomap data / arg only ------------------
alphaDiv_res.water <- metaData %>% filter(`Sample Type` != "River sediment")
alphaDiv_res.sed <- metaData %>% filter(`Sample Type` == "River sediment")
alphaDiv_res <- metaData

alphaDiv_res$invSimpson <- diversity(filtLongLogNoTxWide, index = "invsimpson", groups = alphaDiv_res$Sample_ID)
alphaDiv_res$Shannon <- diversity(filtLongLogNoTxWide, index =  "shannon", groups = alphaDiv_res$Sample_ID)
alphaDiv_res$Simpson <- diversity(filtLongLogNoTxWide, index = "simpson", groups = alphaDiv_res$Sample_ID)

alphaDiv_res.water$invSimpson <- diversity(waterfiltLongLogNoTxWide, index = "invsimpson", groups = alphaDiv_res.water$Sample_ID)
alphaDiv_res.water$Shannon <- diversity(waterfiltLongLogNoTxWide, index =  "shannon", groups = alphaDiv_res.water$Sample_ID)
alphaDiv_res.water$Simpson <- diversity(waterfiltLongLogNoTxWide, index = "simpson", groups = alphaDiv_res.water$Sample_ID)

alphaDiv_res.sed$invSimpson <- diversity(sedfiltLongLogNoTxWide, index = "invsimpson", groups = alphaDiv_res.sed$Sample_ID)
alphaDiv_res.sed$Shannon <- diversity(sedfiltLongLogNoTxWide, index =  "shannon", groups = alphaDiv_res.sed$Sample_ID)
alphaDiv_res.sed$Simpson <- diversity(sedfiltLongLogNoTxWide, index = "simpson", groups = alphaDiv_res.sed$Sample_ID)

hist(alphaDiv_res$invSimps)
hist(alphaDiv_res$Shannon)
hist(alphaDiv_res$Simps)

hist(alphaDiv_res.water$invSimps)
hist(alphaDiv_res.water$Shannon)
hist(alphaDiv_res.water$Simps)

hist(alphaDiv_res.sed$invSimps)
hist(alphaDiv_res.sed$Shannon)
hist(alphaDiv_res.sed$Simps)

color_mapping <- c("Nurpur Upstream" = "cyan2",  "Nurpur Midstream" = "skyblue1",
                   "Nurpur Downstream" ="skyblue3", "Shahdara Upstream" = "burlywood1",
                   "Shahdara Midstream" = "burlywood3", "Shahdara Upstream" = "blue",
                   "Korang Upstream" = "chocolate4", "Korang Midstream" = "slategray1",
                   "Korang Downstream" = "slategray3", "Rawal Lake Upstream" = "slateblue4",
                   "Rawal Lake Midstream" = "pink",  "Rawal Lake Downstream" = "plum",
                   "Jinnah Upstream" = "violetred3",  "Jinnah Midstream" = "magenta4",
                   "Jinnah Downstream" = "violet")

alphaDiv_res_long <- pivot_longer(alphaDiv_res,
                                  cols = c("invSimpson", "Shannon", "Simpson"),
                                  names_to = "alphaIndex",
                                  values_to = "ARGalphaDiv")
alphaDiv_res_long$Site_ID <- factor(alphaDiv_res_long$Site_ID , levels = unique(alphaDiv_res$Site_ID))

#ALL alpha index graphs
p_adiv2 <- ggplot(alphaDiv_res_long, aes(x = Site_ID, y = ARGalphaDiv, fill = Site_ID)) +
  geom_boxplot() +
  facet_wrap(~alphaIndex, scales = "free", nrow = 3) +
  scale_fill_manual(values = color_mapping) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(fill = "Site_ID")
p_adiv2

ggsave(
  filename = paste(plotsDir, "resistomap_alpha_div_arg.jpeg", sep = ""),
  plot = p_adiv2,
  device = "jpeg",
  scale = 1,
  dpi = 300,
  height = 12
)

test <- wilcox.test(alphaDiv_res$invSimpson ~ alphaDiv_res$`Sample Type`)
test

test <- wilcox.test(alphaDiv_res$Simpson ~ alphaDiv_res$`Sample Type`)
test

test <- wilcox.test(alphaDiv_res$Shannon ~ alphaDiv_res$`Sample Type`)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$River)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$River)
test

test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$River)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$`River Location`)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$`River Location`)
test

test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$`River Location`)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$Site_ID)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$Site_ID)
test

test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$Site_ID)
test


#----- sum relative abundance by siteID to get resistance score

summarized_df <- waterfiltLongLogNoTx %>%
  group_by(Sample_ID, group) %>%
  summarise(total_relative_abundance = sum(relative_abundance,  na.rm = TRUE))
summarized_dfm <- merge(summarized_df, metaData, by = "Sample_ID", all.x = TRUE)

color_mapping2 = c(Aminoglycoside = "#1a1a1a", "Beta Lactam" = "#d53e4f", Integrons = "#3288bd",MDR = "grey49", MGE = "#66c2a5", MLSB = "#abdda4",  Other = "#e6f598" , Phenicol ="#ffffbf",
                   Quinolone = "#fee08b", Sulfonamide ="#fdae61",Tetracycline  ="#f46d43", Trimethoprim ="#8c510a", Vancomycin = "#9e0142")

summarized_dfm$Site_ID <- factor(summarized_dfm$Site_ID , levels = unique(summarized_dfm$Site_ID))
summarized_dfm$site_color <- color_mapping[summarized_dfm$Site_ID]

summarized_dfm$group <- factor(summarized_dfm$group , levels=unique(summarized_dfm$group))
summarized_dfm$group_color <- color_mapping[summarized_dfm$group]
summarized_dfm$`Sample Type` <- factor(summarized_dfm$`Sample Type`, levels = unique(summarized_dfm$`Sample Type`))

plot <- ggplot(summarized_dfm, aes(x = Site_ID, y = total_relative_abundance, fill = Site_ID)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free") +
  scale_fill_manual(values = color_mapping) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(fill = "Site_ID")  # Rename the legend title

plot
ggsave(paste(plotsDir, "RM ARGs abundance.pdf", sep = ""), plot, height = 8, width = 12)

plot <- ggplot(summarized_dfm, aes(x = Site_ID, y = total_relative_abundance, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  #facet_wrap(~Samp_time, nrow = 3) +
  scale_fill_manual(values = color_mapping2) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 21),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(fill = "group")  # Rename the legend title
plot



#---- summarise dataset water only
summarized_dfmb <- filter(summarized_dfm, `Sample Type` == "River water") %>%
  group_by(Sample_ID ) %>%
  summarise(total_ARG_relAb = sum(total_relative_abundance,  na.rm = TRUE))  %>%
  mutate(Log_RA = log10(total_ARG_relAb))
summarized_dfmb2 <- merge(summarized_dfmb, metaData, by = "Sample_ID", all.x = TRUE)


summarized_dfmb2$Site_ID <- factor(summarized_dfmb2$Site_ID , levels = unique(summarized_dfmb2$Site_ID))
summarized_dfmb2$site_color <- color_mapping[summarized_dfmb2$Site_ID]

#create a static map of the area of interest

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid

# get_stadiamap now requires an API key. To register the API run the following command only once. This will save the API key into the R environment for future request.
#register_stadiamaps("API_Key_avaialble on your account", write = TRUE)

map <- qmplot(Longitude, Latitude, data = summarized_dfmb2, maptype = "stamen_terrain")
#map <- qmplot(Longitude, Latitude, data = summarized_dfmb2, maptype = "stamen_toner")

map

map +
  geom_point(data = summarized_dfmb2, mapping = aes(x = Longitude, y = Latitude, size = total_ARG_relAb, color = Site_ID)) +
  scale_color_manual(values = color_mapping) +
  theme_light() +
  #scale_color_distiller(name  ="HuBac", palette = "OrRd", direction = 1) +
  scale_size_continuous(name  ="ARG\n(total relative\nabundace)") +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16)) + # add legend at right of plot
  geom_text_repel(aes(x = Longitude, y = Latitude, label = Site_ID),
                  box.padding = 0.5, point.padding = 0.1, segment.color = 'black',
                  size = 5, color = "black") +
  facet_wrap(~ Samp_time, nrow = 2) +  theme(strip.text = element_text(size = 16))

map +
  geom_point(data = summarized_dfmb2, mapping = aes(x = Longitude, y = Latitude, size = Log_RA, color = Site_ID)) +
  scale_color_manual(values = color_mapping) +
  theme_light() +
  #scale_color_distiller(name  ="HuBac", palette = "OrRd", direction = 1) +
  scale_size_continuous(name  ="Log10 ARG\n(total relative\nabundace)") +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16)) + # add legend at right of plot
  geom_text_repel(aes(x = Longitude, y = Latitude, label = Site_ID),
                  box.padding = 0.5, point.padding = 0.1, segment.color = 'black',
                  size = 5, color = "black") +
  facet_wrap(~ Samp_time, nrow = 2) +  theme(strip.text = element_text(size = 16))

#ggsave(plot = last_plot(), filename = "map_arg_resistomap_Log10totalRA.jpeg", device = "jpeg")




#---- summarise dataset sediment only

summarized_dfmk <- filter(summarized_dfm, Location == "Kangra") %>%
  group_by(sample_ID ) %>%
  summarise(total_ARG_relAb = sum(total_relative_abundance,  na.rm = TRUE))  %>%
  mutate(Log_RA = log10(summarized_dfmk$total_ARG_relAb))
summarized_dfmk2 <- merge(summarized_dfmk, metadata, by = "sample_ID", all.x = TRUE)


summarized_dfmk2$Site_ID <- factor(summarized_dfmk2$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
summarized_dfmk2$site_color <- color_mapping[summarized_dfmk2$Site_ID]
summarized_dfmk2$Samp_time <- factor(summarized_dfmk2$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))

#create a static map of the area of interest

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid

# get_stadiamap now requires an API key. To register the API run the following command only once. This will save the API key into the R environment for future request.
#register_stadiamaps("API_Key_avaialble on your account", write = TRUE)

map <- qmplot(Longitude, Latitude, data = summarized_dfmk2, maptype = "stamen_terrain")
#map <- qmplot(Longitude, Latitude, data = summarized_dfmb2, maptype = "stamen_toner")

map

map +
  geom_point(data = summarized_dfmk2, mapping = aes(x = Longitude, y = Latitude, size = total_ARG_relAb, color = Site_ID)) +
  scale_color_manual(values = color_mapping) +
  theme_light() +
  #scale_color_distiller(name  ="HuBac", palette = "OrRd", direction = 1) +
  scale_size_continuous(name  ="ARG\n(total relative\nabundace)") +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16)) + # add legend at right of plot
  geom_text_repel(aes(x = Longitude, y = Latitude, label = Site_ID),
                  box.padding = 0.5, point.padding = 0.1, segment.color = 'black',
                  size = 5, color = "black") +
  facet_wrap(~ Samp_time, nrow = 2) +  theme(strip.text = element_text(size = 16))

map +
  geom_point(data = summarized_dfmk2, mapping = aes(x = Longitude, y = Latitude, size = Log_RA, color = Site_ID)) +
  scale_color_manual(values = color_mapping) +
  theme_light() +
  #scale_color_distiller(name  ="HuBac", palette = "OrRd", direction = 1) +
  scale_size_continuous(name  ="Log10 ARG\n(total relative\nabundace)") +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16)) + # add legend at right of plot
  geom_text_repel(aes(x = Longitude, y = Latitude, label = Site_ID),
                  box.padding = 0.5, point.padding = 0.1, segment.color = 'black',
                  size = 5, color = "black") +
  facet_wrap(~ Samp_time, nrow = 2) +  theme(strip.text = element_text(size = 16))


#------ exploring Bacteroidets differential distribution
tax_wmtab_2223_llog <- filter(waterfiltLongLog, group == "Taxonomic")
bacteroides_w2223 <- filter(tax_wmtab_2223_llog, gene == "Bacteroidetes")

#Test normality of data
ggdensity(bacteroides_w2223$relative_abundance, fill = "lightgray")
ggqqplot(bacteroides_w2223$relative_abundance)
shapiro.test(bacteroides_w2223$relative_abundance)
#Significant difference tests
t_test_result <- t.test(relative_abundance ~ Location,  data = bacteroides_w2223, var.equal = FALSE )
print(t_test_result)

aov2 <- aov(relative_abundance ~ Samp_time, bacteroides_w2223 )
summary(aov2)
TukeyHSD(aov2)

aov1 <- aov(relative_abundance ~ `Sample Type`, bacteroides_w2223)
summary(aov1)
TukeyHSD(aov1)

wilcox.test(bacteroides_w2223$relative_abundance ~ bacteroides_w2223$Location)
kruskal.test(bacteroides_w2223$relative_abundance ~ bacteroides_w2223$`Sample Type`)
kruskal.test(bacteroides_w2223$relative_abundance ~ bacteroides_w2223$Samp_time)
# Perform Dunn's test with Bonferroni correction
dunn.test(bacteroides_w2223$relative_abundance, g = bacteroides_w2223$`Sample Type`, method = "bonferroni")

#Use Log_RA
#Test normality of data
ggdensity(bacteroides_w2223$Log_RA, fill = "lightgray")
ggqqplot(bacteroides_w2223$Log_RA)
shapiro.test(bacteroides_w2223$Log_RA)
#Significant difference tests
wilcox.test(bacteroides_w2223$Log_RA ~ bacteroides_w2223$Location)
kruskal.test(bacteroides_w2223$Log_RA ~ bacteroides_w2223$`Sample Type`)
kruskal.test(bacteroides_w2223$Log_RA ~ bacteroides_w2223$Samp_time)

# Perform Dunn's test with Bonferroni correction
dunn.test(bacteroides_w2223$Log_RA, g = bacteroides_w2223$`Sample Type`, method = "bonferroni")


#------ mapping Bacteroidets and ARGs/MGEs on map BADDI WATER ONLY

summarized_bacteroides_b <- filter(bacteroides_w2223, Location == "Baddi")  %>%
  group_by(sample_ID ) %>%
  summarise(Bacteroides_relAb = sum(relative_abundance, na.rm = TRUE))

summarized_dfmb3 <- merge(summarized_dfmb2, summarized_bacteroides_b, by = "sample_ID", all.x = TRUE)

kruskal.test(summarized_dfmb3$Log_RA ~ summarized_dfmb3$`Sample Type`)
dunn.test(summarized_dfmb3$Log_RA, g = summarized_dfmb3$`Sample Type`, method = "bonferroni")

kruskal.test(summarized_dfmb3$Log_RA ~ summarized_dfmb3$Samp_time)


summarized_dfmb3$Site_ID <- factor(summarized_dfmb3$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
summarized_dfmb3$site_color <- color_mapping[summarized_dfmb3$Site_ID]
summarized_dfmb3$Samp_time <- factor(summarized_dfmb3$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))

map <- qmplot(Longitude, Latitude, data = summarized_dfmb3, maptype = "stamen_terrain")
#map <- qmplot(Longitude, Latitude, data = summarized_dfmb2, maptype = "stamen_toner")

map

map +
  geom_point(data = summarized_dfmb3, mapping = aes(x = Longitude, y = Latitude, size = total_ARG_relAb, color = Bacteroides_relAb)) +
  #scale_color_manual(values = color_mapping) +
  theme_light() +
  scale_color_distiller(name  ="Bacteroidetes\n(total relative\nabundace)", palette = "OrRd", direction = 1) +
  scale_size_continuous(name  ="ARG\n(total relative\nabundace)") +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16)) + # add legend at right of plot
  geom_text_repel(aes(x = Longitude, y = Latitude, label = Site_ID),
                  box.padding = 0.5, point.padding = 0.1, segment.color = 'black',
                  size = 5, color = "black") +
  facet_wrap(~ Samp_time, nrow = 1) +  theme(strip.text = element_text(size = 16))


#------ mapping Bacteroidets and ARGs/MGEs on map KANGRA WATER ONLY
summarized_bacteroides_k <- filter(bacteroides_w2223, Location == "Kangra")  %>%
  group_by(sample_ID ) %>%
  summarise(Bacteroides_relAb = sum(relative_abundance, na.rm = TRUE))

summarized_dfmk3 <- merge(summarized_dfmk2, summarized_bacteroides_k, by = "sample_ID", all.x = TRUE)


summarized_dfmk3$Site_ID <- factor(summarized_dfmk3$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
summarized_dfmk3$site_color <- color_mapping[summarized_dfmk3$Site_ID]
summarized_dfmk3$Samp_time <- factor(summarized_dfmk3$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))

map <- qmplot(Longitude, Latitude, data = summarized_dfmk3, maptype = "stamen_terrain")
#map <- qmplot(Longitude, Latitude, data = summarized_dfmb2, maptype = "stamen_toner")

map

map +
  geom_point(data = summarized_dfmk3, mapping = aes(x = Longitude, y = Latitude, size = total_ARG_relAb, color = Bacteroides_relAb)) +
  #scale_color_manual(values = color_mapping) +
  theme_light() +
  scale_color_distiller(name  ="Bacteroidetes\n(total relative\nabundace)", palette = "OrRd", direction = 1) +
  scale_size_continuous(name  ="ARG\n(total relative\nabundace)") +
  theme(legend.position = "right", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16)) + # add legend at right of plot
  geom_text_repel(aes(x = Longitude, y = Latitude, label = Site_ID),
                  box.padding = 0.5, point.padding = 0.1, segment.color = 'black',
                  size = 5, color = "black") +
  facet_wrap(~ Samp_time, nrow = 1) +  theme(strip.text = element_text(size = 16))


#------- IMPORT PHYSICO-CHECMICAL AND HEAVY METAL DATASETS FOR 2022/2023 -----# use tidyverse and dplyr
df5 <- read_delim("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/neeri_datasets2024/heavymetal_analysis2024.txt")
df6 <- read_delim("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/neeri_datasets2024/physicochemical_analysis2024.txt")

#df5_long <- pivot_longer(df5,
#                            cols = c("Cr", "Pb", "Cu","Zn",	"Mn",	"Co", "Cd",	"Ni"),
#                            names_to = "Heavy_Metals",
#                            values_to = "HM_ppm")
df5_long <- pivot_longer(df5,
                         cols = c("Cr", "Pb", "Cu","Zn", "Fe",	"Mn", "Al",	"Co", "Cd",	"Ni"),
                         names_to = "Heavy_Metals",
                         values_to = "HM_ppm")

df5_long$Samp_time <- factor(df5_long$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))
df5$Samp_time <- factor(df5$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))
df6$Samp_time <- factor(df6$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))

df5_long$Site_ID <- factor(df5_long$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
df5$Site_ID <- factor(df5$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
df6$Site_ID <- factor(df6$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))


# Define a named vector mapping heavy metals to colors
heavynetals_colors <- c(
  Cr = "green",
  Pb = "darkgrey",
  Cu = "blue",
  Zn = "blanchedalmond",
  Fe = "orange",
  Mn = "purple",
  Al = "azure3",
  Co = "pink",
  Cd = "yellow",
  Ni = "lightseagreen"
)

# Create stacked bar chart
ggplot(df5_long, aes(x = Site_ID, y = HM_ppm, fill = Heavy_Metals)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Samp_time) +
  scale_fill_manual(values = heavynetals_colors) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 21),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(x = "Site ID", y = "ppm", fill = "Heavy Metals") # +
coord_cartesian(ylim = c(0, 50)) # Adjust the limits according to your preference

# Create boxplots
boxplots2 <- lapply(names(df5)[10:19], function(var) {
  ggplot(df5) +
    geom_boxplot(aes(x = Site_ID, y = !!sym(var)), fill = "lightgrey") +
    labs(y = var) +
    theme_light()
})

grid.arrange(grobs = c(boxplots2), ncol = 4)


# Create boxplot for each variable in df6
boxplots <- lapply(names(df6)[10:16], function(var) {
  ggplot(df6) +
    geom_boxplot(aes(x = Site_ID, y = !!sym(var)), fill = "lightgrey") +
    labs(y = var) +
    theme_light()
})

grid.arrange(grobs = c(boxplots), ncol = 4)


# ---------- DATASET prep for modelling--------#


#----- for all the water sample 2022 and 2023, sum the relative abundance by sample_ID to get Relative abundance for mobilome and resistome.

#Add column to subgroup the groups into resistome and mobilome
w_mod_2223 <- waterfiltLongLogNoTx %>% mutate(hgroup = ifelse(group %in% c("Integrons", "MGE"), "Mobilome", "Resistome"))
#Summarise the relative abundance according to resistome and mobilome
sum_w_mod_2223 <- w_mod_2223 %>%
  group_by(sample_ID, hgroup) %>%
  summarise(total_RA = sum(relative_abundance,  na.rm = TRUE))
# Transform in wide format so it is easier to combine different datasets
sum_w_mod_2223a <- pivot_wider(sum_w_mod_2223,
                               names_from = hgroup,
                               values_from = total_RA,
                               names_prefix = c("total_RA_"))
sum_w_mod_2223b <- sum_w_mod_2223a %>% mutate(tot_Log_RA_Resistome = log10(total_RA_Resistome)) %>% mutate(tot_Log_RA_Mobilome = log10(total_RA_Mobilome))  #calculate Log10 of total_RA values

sum_w_mod_2223m <- merge(sum_w_mod_2223b, metadata, by = "sample_ID", all.x = TRUE) #merge metadata based on sample_ID
#Add heavymeatlas metadata
sum_w_mod_2223m1 <- merge(sum_w_mod_2223m, df5[ , c(1,10:19)], by.x = "sample_ID", by.y = "sample_ID", all.x = TRUE) #merge metadata based on sample_ID
#Add physico-chemicals as metadata
sum_w_mod_2223m2 <- merge(sum_w_mod_2223m1, df6[ , c(1,10:16)], by.x = "sample_ID", by.y = "sample_ID", all.x = TRUE) #merge metadata based on sample_ID

#TO add NMDS1 from mobiolme and resistome analysed together
sum_w_mod_2223m3 <- merge(sum_w_mod_2223m2, argOrdMd[ , 1:3], by.x = "sample_ID", by.y = "Row.names", all.x = TRUE) #merge metadata based on sample_ID
sum_w_mod_2223m3 <- sum_w_mod_2223m3 %>% rename(NMDS1_RM = MDS1) %>% rename(NMDS2_RM = MDS2)# RENAME MDS1 column in NMDS1_RM (resistome and mobilome based NMDS)

#TO add NMDS1 from mobiolme and resistome analysed separately
#resistome data only
resistome <- filter(waterfiltLongLogNoTx, !group  %in% c("Integrons", "MGE" ))
rownames(resistome) <- resistome

# transform dataframe with sample in row, arg RA in column and NA into 0
resistome_wide <- pivot_wider(resistome, id_cols = sample_ID ,names_from = gene, values_from = relative_abundance) %>% replace(is.na(.), 0)
resistome_wide <-as.data.frame(resistome_wide)
rownames(resistome_wide) <- resistome_wide$sample_ID
resistome_wide <- resistome_wide[,-1]


#data ordination using bray curtis
arg_ord1 <- metaMDS(resistome_wide, distance = "bray", autotransform = F, k = 3)
arg_ord1 # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good

# to use ggplot the results form the ordination need to be saved into a separate dataframe.
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. this can be done this by calling the scores of the mds.
arg_ord_md1 <- merge(arg_ord1$points, metadata, by.x = "row.names", by.y = "sample_ID", all.x = TRUE)

nmds_resistomeplot <- ggplot(arg_ord_md1, aes(x=MDS1, y=MDS2))+ #sets up the plot
  geom_point(aes(MDS1, MDS2, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  stat_ellipse(aes(colour = factor(Location)), level = 0.95, type = "norm") +  # Add ellipses
  theme_light()+
  #theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Site_ID", shape = "Year")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size=21))
#labs(title = "Ordination with arg vectors")
nmds_resistomeplot

#TO add NMDS1 from mobiolme and resistome analysed together
sum_w_mod_2223m4 <- merge(sum_w_mod_2223m3, arg_ord_md1[ , 1:3], by.x = "sample_ID", by.y = "Row.names", all.x = TRUE) #merge metadata based on sample_ID
sum_w_mod_2223m4 <- sum_w_mod_2223m4 %>% rename(NMDS1_R = MDS1) %>% rename(NMDS2_R = MDS2)# RENAME MDS1 column in NMDS1_R (resistome based NMDS)


#mobiolome data only
mobilome <- filter(waterfiltLongLogNoTx, group  %in% c("Integrons", "MGE" ))
rownames(mobilome) <- mobilome

# transform dataframe with sample in row, arg RA in column and NA into 0
mobilome_wide <- pivot_wider(mobilome, id_cols = sample_ID ,names_from = gene, values_from = relative_abundance) %>% replace(is.na(.), 0)
mobilome_wide <-as.data.frame(mobilome_wide)
rownames(mobilome_wide) <- mobilome_wide$sample_ID
mobilome_wide <- mobilome_wide[,-1]


#data ordination using bray curtis
mob_ord1 <- metaMDS(mobilome_wide, distance = "bray", autotransform = F, k = 3)
mob_ord1 # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good

# to use ggplot the results form the ordination need to be saved into a separate dataframe.
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. this can be done this by calling the scores of the mds.
mob_ord_md1 <- merge(mob_ord1$points, metadata, by.x = "row.names", by.y = "sample_ID", all.x = TRUE)

nmds_mobilomeplot <- ggplot(mob_ord_md1, aes(x=MDS1, y=MDS2))+ #sets up the plot
  geom_point(aes(MDS1, MDS2, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = riverColors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = riverColors) +  # Use defined colors for Site_ID
  stat_ellipse(aes(colour = factor(Location)), level = 0.95, type = "norm") +  # Add ellipses
  theme_light()+
  #theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Site_ID", shape = "Year")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size=21))
#labs(title = "Ordination with arg vectors")
nmds_mobilomeplot

#TO add NMDS1 from mobiolme and resistome analysed together
sum_w_mod_2223m5 <- merge(sum_w_mod_2223m4, mob_ord_md1[ , 1:3], by.x = "sample_ID", by.y = "Row.names", all.x = TRUE) #merge metadata based on sample_ID
sum_w_mod_2223m5 <- sum_w_mod_2223m5 %>% rename(NMDS1_M = MDS1) %>% rename(NMDS2_M = MDS2)# RENAME MDS1 column in NMDS1_R (resistome based NMDS)

#SAVE FINAL FILE for MODELLING
write_csv(sum_w_mod_2223m5, file = "respharm_model_table.csv")

}


