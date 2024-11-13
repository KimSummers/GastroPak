# NMDSGenus
#
# Add sampling site and season to the data with sample id
#
# file NMDSGenus
#
# inputs
# 	metagenomicDataFile - File with data to analyse
#   metaDataEnvFile     - File with metadata

# Version    Author       Date      Affiliation
# 1.00       J K Summers  08/05/24  Wellington Lab - School of Life Sciences - University of Warwick
NMDSGenus <- function(metagenomicDataFile, metaDataEnvFile) {

# Load necessary packages
library(vegan)
library(ggplot2)
library(ggpubr)

# Read in data
mgGenusWide <- read.csv(metagenomicDataFile, sep = "\t")

metaEnv <- read.csv(metaDataEnvFile)

# 1st column is id of sample (by no 1-16) so need to remove this
mgWideGenus <- mgWideGenus[, -1]

# Calculate Shannon diversity index
shannon <- diversity(mgWideGenus, index = "shannon")

diversity_genus<-cbind(wimp_wide_env, shannon=shannon)

# Convert the Shannon diversity index to a dataframe
shannon_df <- as.data.frame(diversity_genus)

# Calculate Shannon diversity index
simpson <- diversity(mgWideGenus, index = "simpson")

diversity_genus<-cbind(diversity_genus, simpson=simpson)

# Convert the Shannon diversity index to a dataframe
shannon_df <- as.data.frame(diversity_genus)


# Create the bar plot
ggplot(shannon_df, aes(x = treatment, y = shannon, fill=River_loc)) +
  geom_col(position = "dodge") +
  theme_minimal()+
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4", "coral", "coral4")) + # Set colors
  labs(x = "Treatment", y = "Shannon Diversity", fill = "Location")

# Create the bar plot
ggplot(shannon_df, aes(x = treatment, y = simpson, fill=River_loc)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("deepskyblue", "deepskyblue4", "coral", "coral4")) + # Set colors
  labs(x = "Treatment", y = "Simpson Diversity", fill = "Location") +
  theme_minimal()

#Shannon Index:
#Higher values indicate more diverse communities, suggesting a balanced distribution of different genera.
#Changes in Shannon index might reflect changes in community structure or ecosystem stability.
#A decrease in Shannon index might indicate disturbances, such as pollution, habitat destruction, or shifts in environmental conditions.
#Simpson Index:
#  Higher values indicate less diverse communities, suggesting dominance by one or a few genera.
#It can help identify highly abundant or dominant genera within the community.
#An increase in Simpson index might indicate the loss of diversity or dominance by a few genera.


#rarefaction curves







# Convert data to a numeric matrix
data_matrix <- as.matrix(mgWideGenus)

# Produce rarefaction curve
rare_curve <- rarecurve(data_matrix, step = 1, col = "blue")

# Compute rarefaction curve
rare_curve <- rarefy(data_matrix, step = 1)

# Convert rarefaction curve to a dataframe
rare_df <- as.data.frame(rare_curve$accumresult)

# Plot with ggplot
ggplot(rare_df, aes(x = Sample, y = Sampled, group = Sample)) +
  geom_line() +
  labs(x = "Number of Samples", y = "Number of Species", title = "Rarefaction Curve") +
  theme_minimal()

rarecurve(mgWideGenus, col = "blue")

# Plot the rarefaction curve

plot(rare_curve, xlab = "Number of Samples", ylab = "Number of Species", main = "Rarefaction Curve")


#mgWideGenus is the species abundance data
wimp_genus.rel <-         #transforming to relative abundance
  decostand(mgWideGenus, method = "total")
# hellinger might be better- square root of total and accounts for more high and low numbers

#wimp_wide_env is the meta data associated with this

# Calculate distance matrix - make distance matrix data frame - write .csv
wimp_genus_distmat <-
  vegdist(wimp_genus.rel, method = "bray")
wimp_genus_DistanceMatrix <-
  as.matrix(wimp_genus_distmat, labels = T)
write.csv(wimp_genus_DistanceMatrix, "wimp_genus_distmat.csv")

# Running NMDS in vegan (metaMDS)
wimp_genus_NMS <-
  metaMDS(wimp_genus_distmat,
          distance = "bray",
          k = 2,
          maxit = 999,
          trymax = 250, autotransform=FALSE)
wimp_genus_scrs <-
  `sppscores<-`(wimp_genus_NMS, mgWideGenus)
wimp_genus_cor <-
  cor(mgWideGenus,
      wimp_genus_NMS$points,
      use = "complete.obs",
      method = "pearson")
write.csv(wimp_genus_cor, file = "wimp_genus_PearsonCor.csv")

# Shepards test/goodness of fit
goodness(wimp_genus_NMS)
stressplot(wimp_genus_NMS)


# Extract NMDS coordinates
nmds_coordinates <- scores(wimp_genus_NMS, choices = c(1, 2))

NMDS_genus<-cbind(wimp_wide_env, nmds_coordinates)

genus_spp_fit <- envfit(wimp_genus_NMS, wimp_genus.rel, permutations = 999)

#vectors
genus_scores<- as.data.frame(scores(wimp_genus_scrs, display="species"))
genus_p_val<-as.data.frame(scores(genus_spp_fit, display="vectors"))
genus_scores <- cbind(genus_scores, Species = rownames(genus_scores)) #add species names to dataframe
genus_pval<-as.data.frame(genus_spp_fit$vectors$pvals)
genus_scores_2 <- cbind(genus_scores, pval=genus_pval$`genus_spp_fit$vectors$pvals`) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.genus.scrs <- subset(genus_scores_2, pval<=0.01) #subset data to show species significant at 0.01


#calculate vector (distance from 0,0)

sig.genus.scrs$distance <- sqrt(sig.genus.scrs$NMDS1^2 + sig.genus.scrs$NMDS2^2)

top20_sig.genus.scrs<- head(sig.genus.scrs[order(-sig.genus.scrs$distance), ], 20)


# Define the desired order of levels and corresponding labels for treatment
desired_order <- c("control", "Ab", "PFOA", "PFOA_Ab")
new_labels <- c("Control", "Antibiotics", "PFOA", "Antibiotics + PFOA")

# Define the desired order of levels and corresponding labels for River_loc
desired_order1 <- c("Aire_up", "Aire_down", "Sowe_up", "Sowe_down")
new_labels1 <- c("River Aire upstream", "River Aire downstream", "River Sowe upstream", "River Sowe downstream")

# Convert categorical variable to factor with desired order of levels and new labels

wimp_wide_env<-wimp_wide_env3
wimp_wide_env$treatment<-factor(wimp_wide_env$treatment, levels=desired_order, labels=new_labels)
wimp_wide_env$River_loc<-factor(wimp_wide_env$River_loc, levels=desired_order1, labels=new_labels1)

# Extract NMDS coordinates
#nmds_coordinates1 <- scores(genus_NMS, choices = c(1, 2))


NMDS_genus<-cbind(wimp_wide_env, nmds_coordinates)
NMDS_arg<-cbind(wimp_wide_env, nmds_coordinates1)


# nmds colored by river location and shape by treatment
nmds.plot.genus <- ggplot(NMDS_genus, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(River_loc),
                 shape = factor(treatment)), size = 4)+ #adds site points to plot, shape determined by month, colour determined by site
  stat_ellipse(geom = "polygon", type="t", alpha=0.1,
               mapping = aes(fill=River_loc, color = River_loc))+
  labs(color = "location", shape = "treatment")+ # add legend labels for Location and Month
  ggtitle("Bacterial genus NMDS - Grouping by location")+
  theme_bw()+
  guides(fill=FALSE)+
  scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral", "coral4")) + # Set colors
  scale_shape_manual(values = c(19, 17, 15, 18)) + # Set shapes (19: circle, 17: triangle, 15: square) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))


nmds.plot.genus

#ADD significant species (i.e. arg driving diversity)
nmds.plot.genus3 <- nmds.plot.genus+
  geom_segment(data = sig.genus.scrs,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.genus.scrs,
                           aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25) # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

nmds.plot.genus3


#too many genus driving diversity so how to reduce this? Reduced sig to <0.01 not <0.05 - maybe also use size of vector??! / combine CTX-M etc?

## -------Perform PERMANOVA to check significance of clustering from ordination
# Produce distance matrix
# Load data
# Assuming row names are sample IDs in both data.matrix and sample_metadata

# Calculating distance matrix using vegan function vegdist
bray_dist_mat_arg <- vegdist(genus_wide_b, method = "bray")
genus_dist<-as.data.frame(genus_DistanceMatrix)
genus_dist2<-cbind(barcode=rownames(genus_dist), genus_dist)
meta_dist<-inner_join(wimp_wide_env, genus_dist2, by="barcode")

#permutations for tests
perm<-999
test1<-adonis2(genus_dist~River_loc, data=meta_dist, permutations=perm, method=euclidean)
# sig p=0.001
test2<-adonis2(genus_dist~River*wwtp, data=meta_dist, permutations=perm, method=euclidean)
#River p=0.004, wwtp p=0.001, River*wwtp p=0.003
test3<-adonis2(genus_dist~treatment, data=meta_dist, permutations=perm, method=euclidean)
#not sig p=0.994


# Run pairwise adonis tests (change location/month to check for each)
#pairwise_adonis_results <- pairwise.perm.manova(genus_dist, wimp_wide_env$River_loc, nperm = 999, p.method = "bonferroni")
# pairwise perm manova not working! can't find function
# View the results
#print(pairwise_adonis_results)

}
