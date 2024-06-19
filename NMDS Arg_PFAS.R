# Load necessary packages
library(vegan)
library(ggplot2)
library(ggvegan)
library(ggpubr)

# Set working directory
setwd("/Users/lfslbe/Library/CloudStorage/OneDrive-UniversityofWarwick/PFAS ECR project/R")



# Read in your data

amr_wide_a <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofWarwick/PFAS ECR project/R/amr_wide_a.csv")
View(amr_wide_a)

wimp_wide_env <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofWarwick/PFAS ECR project/R/wimp_wide_env.csv")
View(wimp_wide_env)


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

# preparing ARG data

amr_wide_a[is.na(amr_wide_a)] <- 0
amr_wide_b<-amr_wide_a %>% column_to_rownames("barcode")
# Read in your data
amr.rel <-         #transforming to relative abundance
  decostand(amr_wide_b, method = "total")
# hellinger might be better- square root of total and accounts for more high and low numbers

#wimp_wide_env is the meta data associated with this

# Calculate distance matrix - make distance matrix data frame - write .csv
amr_distmat <- 
  vegdist(amr.rel, method = "bray")
amr_DistanceMatrix <- 
  as.matrix(amr_distmat, labels = T)
write.csv(amr_DistanceMatrix, "amr_distmat.csv")

# Running NMDS in vegan (metaMDS)
amr_NMS <-
  metaMDS(amr_distmat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 250, autotransform=FALSE)
amr_scrs <-
  `sppscores<-`(amr_NMS, amr_wide_b)
amr_cor <-
  cor(amr_wide_b,
      amr_NMS$points,
      use = "complete.obs",
      method = "pearson")
write.csv(amr_cor, file = "amr_PearsonCor.csv")

# Shepards test/goodness of fit
goodness(amr_NMS)
stressplot(amr_NMS)


# Extract NMDS coordinates
nmds_coordinates1 <- scores(amr_NMS, choices = c(1, 2))


NMDS_arg<-cbind(wimp_wide_env, nmds_coordinates1)


arg_spp_fit <- envfit(amr_NMS, amr.rel, permutations = 999)

# obtain all info from NMDS as dataframe


amr_scores<- as.data.frame(scores(amr_scrs, display="species"))
amr_p_val<-as.data.frame(scores(arg_spp_fit, display="vectors"))
amr_scores <- cbind(amr_scores, Species = rownames(amr_scores)) #add species names to dataframe
amr_pval<-as.data.frame(arg_spp_fit$vectors$pvals)
amr_scores <- cbind(amr_scores, pval=amr_pval$`arg_spp_fit$vectors$pvals`) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.arg.scrs <- subset(amr_scores, pval<=0.01) #subset data to show species significant at 0.05

#calculate vector (distance from 0,0)

sig.arg.scrs$distance <- sqrt(sig.arg.scrs$NMDS1^2 + sig.arg.scrs$NMDS2^2)

top20_sig.arg.scrs<- head(sig.arg.scrs[order(-sig.arg.scrs$distance), ], 20)


#sig_vars<-arg_spp_fit$signif

# nmds colored by river location and shape by treatment
nmds.plot.arg <- ggplot(NMDS_arg, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(River_loc), 
                 shape = factor(treatment)), size = 4)+ #adds site points to plot, shape determined by month, colour determined by site
  stat_ellipse(geom = "polygon", type="t", alpha=0.1, 
               mapping = aes(fill=River_loc, color = River_loc))+
  labs(color = "location", shape = "treatment")+ # add legend labels for Location and Month
  ggtitle("Antibiotic resistance gene NMDS - Grouping by location")+
  theme_bw()+ 
  guides(fill=FALSE)+
  scale_color_manual(values = c("deepskyblue", "deepskyblue4", "coral", "coral4", "red","blue", "chartreuse3","lightblue","cornflowerblue","darkslategray","cadetblue4","orange","darkolivegreen","grey","darkblue")) + # Set colors
  scale_shape_manual(values = c(19, 17, 15, 18)) + # Set shapes (19: circle, 17: triangle, 15: square) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))


nmds.plot.arg


#ADD significant species (i.e. arg driving diversity)
nmds.plot.arg3 <- nmds.plot.arg+
  geom_segment(data = sig.arg.scrs,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.arg.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25) # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

nmds.plot.arg3
#too many arg driving diversity so how to reduce this? Reduced sig to <0.01 not <0.05 - maybe also use size of vector??! / combine CTX-M etc?

#ADD significant species (i.e. arg driving diversity)
nmds.plot.arg4 <- nmds.plot.arg+
  geom_segment(data = top20_sig.arg.scrs,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = top20_sig.arg.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25) # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

nmds.plot.arg4

nmds.plot.arg5<-nmds.plot.arg4 + ylim(-1,0.5)
nmds.plot.arg5


nmds.plot.arg6<-nmds.plot.arg3 + ylim(-0.25,0.5)+xlim(0,0.5)
nmds.plot.arg6


# link the name of arg to drug_class etc 

# Create new column 'name' with the same data as 'Species'
sig.arg.scrs$name <- sig.arg.scrs$Species

# add metadata information to the amr_count file
sig.arg.scrs.meta <- inner_join(sig.arg.scrs, amr_groups_metadata2, by = "name")

#ADD significant species (i.e. arg driving diversity)
nmds.plot.arg7 <- nmds.plot.arg+
  geom_segment(data = sig.arg.scrs.meta,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, colour=drug_class), 
               arrow = arrow(length = unit(0.25, "cm")), 
                lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.arg.scrs.meta, 
                           aes(x=NMDS1, y=NMDS2, label = Species, colour=drug_class), cex = 3, direction = "both", segment.size = 0.25) # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

nmds.plot.arg7


nmds.plot.arg8<-nmds.plot.arg7 + ylim(-0.25,0.5)+xlim(0,0.5)
nmds.plot.arg8


# Filter out a particular class, for example, class 'A'
filt_sig.arg.scrs.meta<- subset(sig.arg.scrs.meta, drug_class != "MDR efflux pump")

nmds.plot.arg9 <- nmds.plot.arg+
  geom_segment(data = filt_sig.arg.scrs.meta,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, colour=drug_class), 
               arrow = arrow(length = unit(0.25, "cm")), 
               lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = filt_sig.arg.scrs.meta, 
                           aes(x=NMDS1, y=NMDS2, label = Species, colour=drug_class), cex = 3, direction = "both", segment.size = 0.25) # #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

nmds.plot.arg9

## -------Perform PERMANOVA to check significance of clustering from ordination
# Produce distance matrix
# Load data
# Assuming row names are sample IDs in both data.matrix and sample_metadata

# Calculating distance matrix using vegan function vegdist
bray_dist_mat_arg <- vegdist(amr_wide_b, method = "bray")
amr_dist<-as.data.frame(amr_DistanceMatrix)
amr_dist2<-cbind(barcode=rownames(amr_dist), amr_dist)
meta_dist<-inner_join(wimp_wide_env, amr_dist2, by="barcode")

#permutations for tests
perm<-999
test1<-adonis2(amr_dist~River_loc, data=meta_dist, permutations=perm, method=euclidean) 
# sig p=0.001
test2<-adonis2(amr_dist~River*wwtp, data=meta_dist, permutations=perm, method=euclidean)
#River p=0.004, wwtp p=0.001, River*wwtp p=0.003
test3<-adonis2(amr_dist~treatment, data=meta_dist, permutations=perm, method=euclidean)
#not sig p=0.994


# Run pairwise adonis tests (change location/month to check for each)
#pairwise_adonis_results <- pairwise.perm.manova(amr_dist, wimp_wide_env$River_loc, nperm = 999, p.method = "bonferroni")
# pairwise perm manova not working! can't find function
# View the results
#print(pairwise_adonis_results)

