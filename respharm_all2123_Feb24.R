#---- RESISTOMAP Kangra and Baddi 2021/2023 all (Feb2024)----#

libraries <-c("vegan", "pheatmap", "ggplot2", "dplyr", "sf", "tidyverse",
              "ggmap", "modeest", "ggpubr", "tidyheatmaps", "maps",
              "ggrepel", "cowplot", "gridExtra", "readxl", "pals", "dunn.test")
lapply(libraries, require, character.only = TRUE)

#---- import metadata file and gene annotation file; create colour annotation per groups -----# use tidyverse and dplyr
metadata <- read.delim("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/metadata_geo_update_Feb2024.txt")
gene54 <- read.delim("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/gene54_metadata_Feb2024.txt") %>% column_to_rownames("gene")
gene72 <- read.delim("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/gene72_metadata_Feb2024.txt") %>% column_to_rownames("gene")


#---- import datasets -----# use tidyverse and dplyr
rmData <- read_csv(rmFile)
df2 <- read_csv("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/Relative_abundances_UoW156.csv") %>% select(-c("group", "gene"))
df3 <- read_csv("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/Relative_abundances_UoW229.csv") #%>% select(-c("group", "gene"))
df4 <- read_csv("~/OneDrive - University of Warwick/RESPHARM/resistomap_datasets/Relative_abundances_UoW350.csv") %>% select(-c("group", "gene"))

#---- merge all 4 tables using 54 targets and write a cvs file -----# use tidyverse and dplyr
merged_tables <- Reduce(function(x, y) merge(x, y, by = "assay", all = FALSE), list(df1, df2, df3, df4))
write.csv(merged_tables, file = "merged_tables.csv", row.names = FALSE)

#---- transformed merged_tables into long format -----# use tidyverse and dplyr
merged_tables_l <- pivot_longer(merged_tables,
                                   cols = starts_with("IND"),
                                   names_to = "sample_ID",
                                   values_to = "relative_abundance")

#---STATISTICAL ANALYSIS----
#--- test for normality----
library(ggpubr)
ggdensity(merged_tables_l$relative_abundance, fill = "lightgray")
ggqqplot(merged_tables_l$relative_abundance)
shapiro.test(merged_tables_l$relative_abundance)

# Shapiro-Wilk normality test
# data:  merged_tables_l$relative_abundance
# W = 0.35105, p-value < 2.2e-16
# DATA ARE NOT NORMALLY DISTRIBUTED

# use log10 and NA instead of 0 values and test again for normality
merged_tables_llog <- mutate(merged_tables_l, Log_RA = log10(merged_tables_l$relative_abundance))

ggdensity(merged_tables_llog$Log_RA, fill = "lightgray")
ggqqplot(merged_tables_llog$Log_RA)
shapiro.test(merged_tables_llog$Log_RA)

# Shapiro-Wilk normality test
# data:  merged_tables_llog$Log_RA
# W = 0.99787, p-value = 5.841e-05
# DATA ARE ALMOST PERFECTLY NROMALLY DISTRIBUTED (W=1)

# add sample metadata using sample ID to merge files
merged_tables_llog <- full_join(merged_tables_llog, metadata, by = "sample_ID")

#---BOXPLOTS all targets
plot <- ggboxplot(merged_tables_llog, x = "sample_ID", y = "Log_RA", color = "Sample_Type")
facet(plot, facet.by = "group", nrow = 4, ncol = 3) +
  rotate_x_text(angle = 90)

plot <- ggboxplot(merged_tables_llog, x = "Sample_Type", y = "Log_RA", color = "Site_ID")
facet(plot, facet.by = "group", nrow = 4, ncol = 3) +
  rotate_x_text(angle = 90)

plot <- ggboxplot(merged_tables_llog, x = "Sample_Type", y = "Log_RA", color = "Location")
facet(plot, facet.by = "group", nrow = 4, ncol = 3) +
  rotate_x_text(angle = 90)

plot <- ggboxplot(merged_tables_llog, x = "Site_ID", y = "Log_RA", color = "Sample_Type")
facet(plot, facet.by = "group", nrow = 4, ncol = 3) +
  rotate_x_text(angle = 90)

#---two-way anova---- MIGHT NOT BE THE BEST STATS THOUGH
res.aov2 <- aov(data = merged_tables_llog, Log_RA ~ group * Location * Site_ID * Sample_Type * Season * Year)
summary(res.aov2)
TukeyHSD(res.aov2)


#---- merge tables from 2022 and 2023 only, add the metadata, remove data from 2021, add LOG10_RA and put in long format -----# use tidyverse and dplyr
mtab_2223 <- Reduce(function(x, y) merge(x, y, by = "assay", all = FALSE), list(df3, df4))
write.csv(mtab_2223, file = "mtab_2223.csv", row.names = FALSE)

#remove genes that had no detection across all samples or were only detected in 1-2 samples
mtab_2223f <- filter(mtab_2223, !assay %in% c("AY475", "AY476","AY478","AY479", "AY480", "AY222","AY197", "AY198", "AY521" ))

mtab_2223_l <- pivot_longer(mtab_2223f,
                            cols = starts_with("IND"),
                            names_to = "sample_ID",
                            values_to = "relative_abundance")
mtab_2223_llog <- mutate(mtab_2223_l, Log_RA = log10(mtab_2223_l$relative_abundance))
mtab_2223_llog <- full_join(mtab_2223_llog, metadata, by = "sample_ID")
mtab_2223_llog <- filter(mtab_2223_llog, Year != "Y2021")

#---- divide water and sediment
wmtab_2223_llog <- filter(mtab_2223_llog, Sample_Type != "River_sediment")
smtab_2223_llog <- filter(mtab_2223_llog, Sample_Type == "River_sediment")

#---- divide water and sediment
bwmtab_2223_llog <- filter(wmtab_2223_llog, Location != "Kangra")


#---- additional objects
ann_colors <- list(
  Location = c(Baddi = "lightcyan3", Kangra = "lightpink2"),
  Sample_Type = c(Drinking_water = "#c7eae5", River_sediment = "#8c510a" , River_water = "#35978f", Sewage = "#bf812d", Waste_stream = "#dfc27d"),
  group = c(Aminoglycoside = "#1a1a1a", "Beta Lactam" = "#d53e4f", Integrons = "#3288bd",MDR = "grey49", MGE = "#66c2a5", MLSB = "#abdda4",  Other = "#e6f598" , Phenicol ="#ffffbf",
                       Quinolone = "#fee08b", Sulfonamide ="#fdae61",Tetracycline  ="#f46d43", Trimethoprim ="#8c510a", Vancomycin = "#9e0142", Taxonomic = "deeppink"),
  Site_ID = c(B1 = "cyan2",  B7 = "skyblue1",  B14 ="skyblue3", B21 = "burlywood1", B23 = "burlywood3", B27 = "blue", B28 = "chocolate4", B30 = "slategray1", B34 = "slategray3", B36 = "slateblue4", K1 = "pink",  K2 = "plum", K3 = "violetred3",  K4 = "magenta4",K5 = "violet"),
  Year = c(Y2023 = 'snow2', Y2022 = 'navy'),
  Season = c(Wet = 'lightblue2', Dry = 'darkorange3'))

ann_colors_w <- list(
  Location = c(Baddi = "lightcyan3", Kangra = "lightpink2"),
  Sample_Type = c(Drinking_water = "#c7eae5", River_water = "#35978f", Sewage = "#bf812d", Waste_stream = "#dfc27d"),
  group = c(Aminoglycoside = "#1a1a1a", "Beta Lactam" = "#d53e4f", Integrons = "#3288bd",MDR = "grey49", MGE = "#66c2a5", MLSB = "#abdda4",  Other = "#e6f598" , Phenicol ="#ffffbf",
            Quinolone = "#fee08b", Sulfonamide ="#fdae61",Tetracycline  ="#f46d43", Trimethoprim ="#8c510a", Vancomycin = "#9e0142", Taxonomic = "deeppink"),
  Site_ID = c(B1 = "cyan2",  B7 = "skyblue1",  B14 ="skyblue3", B21 = "burlywood1", B23 = "burlywood3", B27 = "blue", B28 = "chocolate4", B30 = "slategray1", B34 = "slategray3", B36 = "slateblue4", K1 = "pink",  K2 = "plum", K3 = "violetred3",  K4 = "magenta4",K5 = "violet"),
  Year = c(Y2023 = 'snow2', Y2022 = 'navy'),
  Season = c(Wet = 'lightblue2', Dry = 'darkorange3'))

ann_colors_s <- list(
  group = c(Aminoglycoside = "#1a1a1a", "Beta Lactam" = "#d53e4f", Integrons = "#3288bd",MDR = "grey49", MGE = "#66c2a5", MLSB = "#abdda4",  Other = "#e6f598" , Phenicol ="#ffffbf",
            Quinolone = "#fee08b", Sulfonamide ="#fdae61",Tetracycline  ="#f46d43", Trimethoprim ="#8c510a", Vancomycin = "#9e0142", Taxonomic = "deeppink"),
  Site_ID = c(B1 = "cyan2",  B7 = "skyblue1",  B14 ="skyblue3", B21 = "burlywood1", B23 = "burlywood3", B27 = "blue", B28 = "chocolate4", B30 = "slategray1", B34 = "slategray3", B36 = "slateblue4"),
  Season = c(Wet = 'lightblue2', Dry = 'darkorange3'))

ann_colors_wb <- list(
  Sample_Type = c(Drinking_water = "#c7eae5", River_water = "#35978f", Sewage = "#bf812d", Waste_stream = "#dfc27d"),
  group = c(Aminoglycoside = "#1a1a1a", "Beta Lactam" = "#d53e4f", Integrons = "#3288bd",MDR = "grey49", MGE = "#66c2a5", MLSB = "#abdda4",  Other = "#e6f598" , Phenicol ="#ffffbf",
            Quinolone = "#fee08b", Sulfonamide ="#fdae61",Tetracycline  ="#f46d43", Trimethoprim ="#8c510a", Vancomycin = "#9e0142", Taxonomic = "deeppink"),
  Site_ID = c(B1 = "cyan2",  B7 = "skyblue1",  B14 ="skyblue3", B21 = "burlywood1", B23 = "burlywood3", B27 = "blue", B28 = "chocolate4", B30 = "slategray1", B34 = "slategray3", B36 = "slateblue4"),
  Year = c(Y2023 = 'snow2', Y2022 = 'navy'),
  Season = c(Wet = 'lightblue2', Dry = 'darkorange3'))


colorRampPalette(c('#2166ac', '#67a9cf','#d1e5f0','#fddbc7', '#ef8a62','#b2182b'))
colorRampPalette(c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac'))

mtab_2223_llog %>%
  arrange(group) %>%
tidyheatmap(
  mtab_2223_llog,
  rows = gene,
  columns =  sample_ID,
  values =  Log_RA,
  #scale = "row",
  color_legend_n = 6,
  color_legend_min = -5,
  color_legend_max = 1,
  annotation_col = c(Location, Sample_Type, Year, Season, Site_ID),
  annotation_row = c(group),
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  gaps_row = group
)

wmtab_2223_llog %>%
  arrange(group) %>%
  tidyheatmap(
    wmtab_2223_llog,
    rows = gene,
    columns =  sample_ID,
    values =  Log_RA,
    #scale = "row",
    color_legend_n = 6,
    color_legend_min = -5,
    color_legend_max = 1,
    annotation_col = c(Location, Sample_Type, Year, Season, Site_ID),
    annotation_row = c(group),
    annotation_colors = ann_colors_w,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    #gaps_row = group,
    cellwidth = 12,
    cellheight = 7
  )

smtab_2223_llog %>%
  arrange(group) %>%
  tidyheatmap(
    smtab_2223_llog,
    rows = gene,
    columns =  sample_ID,
    values =  Log_RA,
    #scale = "row",
    color_legend_n = 6,
    color_legend_min = -5,
    color_legend_max = 1,
    annotation_col = c( Season, Site_ID),
    annotation_row = c(group),
    annotation_colors = ann_colors_s,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    gaps_row = group,
    cellwidth = 12,
    cellheight = 7
  )


bwmtab_2223_llog %>%
  arrange(group) %>%
  tidyheatmap(
    bwmtab_2223_llog,
    rows = gene,
    columns =  sample_ID,
    values =  Log_RA,
    #scale = "row",
    color_legend_n = 6,
    color_legend_min = -5,
    color_legend_max = 1,
    annotation_col = c(Sample_Type, Year, Season, Site_ID),
    annotation_row = c(group),
    annotation_colors = ann_colors_wb,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    gaps_row = group,
    cellwidth = 12,
    cellheight = 7
  )

bwmtab_2223_llog %>%
  arrange(group) %>%
  tidyheatmap(
    bwmtab_2223_llog,
    rows = gene,
    columns =  sample_ID,
    values =  Log_RA,
    #scale = "row",
    color_legend_n = 6,
    color_legend_min = -5,
    color_legend_max = 1,
    annotation_col = c(Sample_Type, Year, Season, Site_ID),
    annotation_row = c(group),
    annotation_colors = ann_colors_wb,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    gaps_row = group,
    gaps_col = Season,
    cellwidth = 12,
    cellheight = 7,
    #show_selected_row_labels = c("intI1_1", "intI1_2", "intI3")
  )





#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  ---- NMDS and PCoA / alpha diversity ----
# CHECK THIS LINK https://www.rpubs.com/RGrieger/545184 for further details of these steps

library(vegan)
#remove taxomonomic data
wmtab_2223_llog_notax <- filter(wmtab_2223_llog, group != "Taxonomic")
rownames(wmtab_2223_llog_notax) <- wmtab_2223_llog_notax

# transform dataframe with sample in row, arg RA in column and NA into 0
wmtab_2223_notax_wide <- pivot_wider(wmtab_2223_llog_notax, id_cols = sample_ID ,names_from = gene, values_from = relative_abundance) %>% replace(is.na(.), 0)
wmtab_2223_notax_wide <-as.data.frame(wmtab_2223_notax_wide)
rownames(wmtab_2223_notax_wide) <- wmtab_2223_notax_wide$sample_ID
wmtab_2223_notax_wide <- wmtab_2223_notax_wide[,-1]


#data ordination using bray curtis
arg_ord <- metaMDS(wmtab_2223_notax_wide, distance = "bray", autotransform = F, k = 3)
arg_ord # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
arg_spp_fit <- envfit(arg_ord, wmtab_2223_notax_wide, permutations = 999)

# to use ggplot the results form the ordination need to be saved into a separate dataframe.
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. this can be done this by calling the scores of the mds.
arg_ord_md <- merge(arg_ord$points, metadata, by.x = "row.names", by.y = "sample_ID", all.x = TRUE)

# A new dataset containing species data also needs to be made to look at species vectors.
spp.scrs <- as.data.frame(scores(arg_spp_fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = arg_spp_fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

head(spp.scrs)

#To show environmental extrinsic variables another datasheet needs to be created. CHECK LINK https://www.rpubs.com/RGrieger/545184

#Now you have the relevant information for plotting the ordination in ggplot! Lets get plotting!
site_colors <- ann_colors_w$Site_ID[match(arg_ord_md$Site_ID, names(ann_colors_w$Site_ID))]

nmds.plot.arg <- ggplot(arg_ord_md, aes(x=MDS1, y=MDS2))+ #sets up the plot
  geom_point(aes(MDS1, MDS2, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = site_colors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = site_colors) +  # Use defined colors for Site_ID
  stat_ellipse(aes(colour = factor(Location)), level = 0.95, type = "norm") +  # Add ellipses
  theme_light()+
  #theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Site_ID", shape = "Year")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size=21))
  #labs(title = "Ordination with arg vectors")
nmds.plot.arg
#ADD significate species (i.e. arg driving diversity) with the following code
nmds.plot.arg +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap


nmds.plot.arg2 <- ggplot(arg_ord_md, aes(x=MDS2, y=MDS3))+ #sets up the plot
  geom_point(aes(MDS2, MDS3, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = site_colors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = site_colors) +  # Use defined colors for Site_ID
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Site_ID", shape = "Year")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size=21))
nmds.plot.arg2


nmds.plot.arg3 <- ggplot(arg_ord_md, aes(x=MDS1, y=MDS3))+ #sets up the plot
  geom_point(aes(MDS1, MDS3, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = site_colors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = site_colors) +  # Use defined colors for Site_ID
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Site_ID", shape = "Year")+ # add legend labels for Location and Month
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 21), axis.text = element_text(size = 18)) + # add legend at right of plot
  theme(text=element_text(size=21))
nmds.plot.arg3

combined_plots <- ggarrange(nmds.plot.arg, nmds.plot.arg2, nmds.plot.arg3, ncol = 3)
ggsave("combined_plots.jpeg", combined_plots, width = 30, height = 9)

# Open a PDF file to save the plots
pdf("combined_plots.pdf")

# Arrange plots in a 3x1 grid
par(mfrow = c(3, 1))

# Print the last three plots
print(nmds.plot.arg)
print(nmds.plot.arg2)
print(nmds.plot.arg3)

# Close the PDF file
dev.off()

# function for ellipsess

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the Site_ID (location) factor
df_ell.arg.location <- data.frame() #sets up a data frame before running the function.

for(g in unique(site.scrs$Site_ID)){
  df_ell.arg.location <- rbind(df_ell.arg.location, cbind(as.data.frame(with(site.scrs [site.scrs$Site_ID==g,],
                                                                             veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) , Site_ID=g))
}

# data for labelling the ellipse
NMDS.mean.arg=aggregate(site.scrs[ ,c("NMDS1", "NMDS2")],
                        list(group = site.scrs$Site_ID), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(site.scrs[,c("NMDS1", "NMDS2")],
                    list(group = site.scrs$Site_ID), mean)
nmds.plot.arg+
  geom_path(data = df_ell.arg.location, aes(x = NMDS1, y = NMDS2, group = Site_ID, colour = Site_ID)) + #this is the ellipse, seperate ones by Site.
  theme(text=element_text(size=21))

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
alphaDiv_res <- metadata %>% filter(Year != "Y2021") %>% filter(Sample_Type != "River_sediment")

alphaDiv_res$invSimpson <- diversity(wmtab_2223_notax_wide, index = "invsimpson", groups = alphaDiv_res$sample_ID)
alphaDiv_res$Shannon <- diversity(wmtab_2223_notax_wide, index =  "shannon", groups = alphaDiv_res$sample_ID)
alphaDiv_res$Simpson <- diversity(wmtab_2223_notax_wide, index = "simpson", groups = alphaDiv_res$sample_ID)

hist(alphaDiv_res$invSimps)
hist(alphaDiv_res$Shannon)
hist(alphaDiv_res$Simps)

alphaDiv_res_long <- pivot_longer(alphaDiv_res,
                                  cols = c("invSimpson", "Shannon", "Simpson"),
                                  names_to = "alphaIndex",
                                  values_to = "ARGalphaDiv")
alphaDiv_res_long$Site_ID <- factor(alphaDiv_res_long$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))

#ALL alpha index graphs
p_adiv2 <- ggplot(alphaDiv_res_long, aes(x = Site_ID, y = ARGalphaDiv, fill = Site_ID)) +
  geom_boxplot() +
  facet_wrap(~alphaIndex, scales = "free", nrow = 3) +
  scale_fill_manual(values = color_mapping) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 21),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(fill = "Site_ID")
p_adiv2

ggsave(
  filename = "resistomap_alpha_div_arg.jpeg",
  plot = p_adiv2,
  device = "jpeg",
  scale = 1,
  dpi = 300,
)

test <- wilcox.test(alphaDiv_res$invSimpson ~ alphaDiv_res$Location)
test

test <- wilcox.test(alphaDiv_res$Simpson ~ alphaDiv_res$Location)
test

test <- wilcox.test(alphaDiv_res$Shannon ~ alphaDiv_res$Location)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$Sample_Type)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$Sample_Type)
test

test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$Sample_Type)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$Samp_time)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$Samp_time)
test
test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$Samp_time)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$Sample_Type)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$Sample_Type)
test

test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$Sample_Type)
test

test <- kruskal.test(alphaDiv_res$Shannon ~ alphaDiv_res$Site_ID)
test

test <- kruskal.test(alphaDiv_res$Simpson ~ alphaDiv_res$Site_ID)
test
test <- kruskal.test(alphaDiv_res$invSimpson ~ alphaDiv_res$Site_ID)
test


#----- sum relative abundance by siteID to get resistance score

summarized_df <- wmtab_2223_llog_notax %>%
  group_by(sample_ID, group) %>%
  summarise(total_relative_abundance = sum(relative_abundance,  na.rm = TRUE))
summarized_dfm <- merge(summarized_df, metadata, by = "sample_ID", all.x = TRUE)

color_mapping <- c(B1 = "cyan2",  B7 = "skyblue1",  B14 ="skyblue3", B21 = "burlywood1", B23 = "burlywood3", B27 = "blue", B28 = "chocolate4", B30 = "slategray1", B34 = "slategray3", B36 = "slateblue4", K1 = "pink",  K2 = "plum", K3 = "violetred3",  K4 = "magenta4",K5 = "violet")
color_mapping2 = c(Aminoglycoside = "#1a1a1a", "Beta Lactam" = "#d53e4f", Integrons = "#3288bd",MDR = "grey49", MGE = "#66c2a5", MLSB = "#abdda4",  Other = "#e6f598" , Phenicol ="#ffffbf",
          Quinolone = "#fee08b", Sulfonamide ="#fdae61",Tetracycline  ="#f46d43", Trimethoprim ="#8c510a", Vancomycin = "#9e0142")

summarized_dfm$Site_ID <- factor(summarized_dfm$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
summarized_dfm$site_color <- color_mapping[summarized_dfm$Site_ID]

summarized_dfm$group <- factor(summarized_dfm$group , levels=c("Aminoglycoside","Beta Lactam","Integrons","MDR","MGE","MLSB","Other","Phenicol","Quinolone","Sulfonamide","Tetracycline","Trimethoprim","Vancomycin" ))
summarized_dfm$group_color <- color_mapping[summarized_dfm$group]
summarized_dfm$Samp_time <- factor(summarized_dfm$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))


plot <- ggplot(summarized_dfm, aes(x = Site_ID, y = total_relative_abundance, fill = Site_ID)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free") +
  scale_fill_manual(values = color_mapping) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 21),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(fill = "Site_ID")  # Rename the legend title

plot

plot <- ggplot(summarized_dfm, aes(x = Site_ID, y = total_relative_abundance, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  #facet_wrap(~Samp_time, nrow = 3) +
  scale_fill_manual(values = color_mapping2) +  # Specify the colors manually
  theme_light() +
  theme(text = element_text(size = 21),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Adjust angle and horizontal justification
  labs(fill = "group")  # Rename the legend title
plot



#---- summarise dataset BADDI water only

summarized_dfmb <- filter(summarized_dfm, Location == "Baddi") %>%
  group_by(sample_ID ) %>%
  summarise(total_ARG_relAb = sum(total_relative_abundance,  na.rm = TRUE))  %>%
  mutate(Log_RA = log10(summarized_dfmb$total_ARG_relAb))
summarized_dfmb2 <- merge(summarized_dfmb, metadata, by = "sample_ID", all.x = TRUE)


summarized_dfmb2$Site_ID <- factor(summarized_dfmb2$Site_ID , levels=c("B1", "B7", "B14", "B21", "B23","B27", "B28", "B30", "B34", "B36", "K1", "K2", "K3", "K4", "K5"))
summarized_dfmb2$site_color <- color_mapping[summarized_dfmb2$Site_ID]
summarized_dfmb2$Samp_time <- factor(summarized_dfmb2$Samp_time , levels=c("Dry_2022", "Wet_2022", "Dry_2023"))

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




#---- summarise dataset KANGRA water only

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
tax_wmtab_2223_llog <- filter(wmtab_2223_llog, group == "Taxonomic")
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

aov1 <- aov(relative_abundance ~ Sample_Type, bacteroides_w2223)
summary(aov1)
TukeyHSD(aov1)

wilcox.test(bacteroides_w2223$relative_abundance ~ bacteroides_w2223$Location)
kruskal.test(bacteroides_w2223$relative_abundance ~ bacteroides_w2223$Sample_Type)
kruskal.test(bacteroides_w2223$relative_abundance ~ bacteroides_w2223$Samp_time)
# Perform Dunn's test with Bonferroni correction
dunn.test(bacteroides_w2223$relative_abundance, g = bacteroides_w2223$Sample_Type, method = "bonferroni")

#Use Log_RA
#Test normality of data
ggdensity(bacteroides_w2223$Log_RA, fill = "lightgray")
ggqqplot(bacteroides_w2223$Log_RA)
shapiro.test(bacteroides_w2223$Log_RA)
#Significant difference tests
wilcox.test(bacteroides_w2223$Log_RA ~ bacteroides_w2223$Location)
kruskal.test(bacteroides_w2223$Log_RA ~ bacteroides_w2223$Sample_Type)
kruskal.test(bacteroides_w2223$Log_RA ~ bacteroides_w2223$Samp_time)

# Perform Dunn's test with Bonferroni correction
dunn.test(bacteroides_w2223$Log_RA, g = bacteroides_w2223$Sample_Type, method = "bonferroni")


#------ mapping Bacteroidets and ARGs/MGEs on map BADDI WATER ONLY

summarized_bacteroides_b <- filter(bacteroides_w2223, Location == "Baddi")  %>%
  group_by(sample_ID ) %>%
  summarise(Bacteroides_relAb = sum(relative_abundance, na.rm = TRUE))

summarized_dfmb3 <- merge(summarized_dfmb2, summarized_bacteroides_b, by = "sample_ID", all.x = TRUE)

kruskal.test(summarized_dfmb3$Log_RA ~ summarized_dfmb3$Sample_Type)
dunn.test(summarized_dfmb3$Log_RA, g = summarized_dfmb3$Sample_Type, method = "bonferroni")

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
w_mod_2223 <- wmtab_2223_llog_notax %>% mutate(hgroup = ifelse(group %in% c("Integrons", "MGE"), "Mobilome", "Resistome"))
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
sum_w_mod_2223m3 <- merge(sum_w_mod_2223m2, arg_ord_md[ , 1:3], by.x = "sample_ID", by.y = "Row.names", all.x = TRUE) #merge metadata based on sample_ID
sum_w_mod_2223m3 <- sum_w_mod_2223m3 %>% rename(NMDS1_RM = MDS1) %>% rename(NMDS2_RM = MDS2)# RENAME MDS1 column in NMDS1_RM (resistome and mobilome based NMDS)

#TO add NMDS1 from mobiolme and resistome analysed separately
#resistome data only
resistome <- filter(wmtab_2223_llog_notax, !group  %in% c("Integrons", "MGE" ))
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
  geom_point(aes(MDS1, MDS2, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = site_colors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = site_colors) +  # Use defined colors for Site_ID
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
mobilome <- filter(wmtab_2223_llog_notax, group  %in% c("Integrons", "MGE" ))
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
  geom_point(aes(MDS1, MDS2, colour = factor(Site_ID), shape = factor(Year)), size = 4, fill = site_colors) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  # coord_fixed()+ #add or exclude this command if you want aspect ratio fixed or variable
  scale_color_manual(values = site_colors) +  # Use defined colors for Site_ID
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




