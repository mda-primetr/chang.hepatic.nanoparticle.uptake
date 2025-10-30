library(phyloseq)
library(dplyr)

# Load data as a phyloseq object
physeq <- readRDS("data/physeq_rarefied.rds")

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame()

# Create Brays Curtis dissimilarity matrix
library(vegan)
dist <- vegdist(t(as.matrix(otu_table(physeq))))

# PCOA analysis and plot ----
## PCOA ellipse plot data preparation and analysis ----
#Conduct PCoA on the dissimilarities
library(ape)
pcoa <- ape::pcoa(dist)

#Merge first two pcoa axes with metadata
pcoa_df <- data.frame(sample = row.names(pcoa$vectors),
                         axis_1 = pcoa$vectors[,1],
                         axis_2 = pcoa$vectors[,2])%>%
  right_join(metadata, "sample")


#Add points for centroids for plotting
pcoa_df <- pcoa_df%>%
  group_by(treat)%>%
  mutate(centroid.x = mean(axis_1),
         centroid.y = mean(axis_2))%>%
  ungroup()

#Run permanova
adonis_out <- adonis2(dist ~ treat, pcoa_df)


#Label for PERMANOVA results
adonis_annot <- paste("Adonis2 test: R2 =", formatC(adonis_out$R2[1], digits = 2), 
                      "p =", adonis_out$`Pr(>F)`[1] )



## Create PCOA ellipse plot ----
##set colorblind friendly palette for plot 
cbpalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#F781BF", "#A65628" )

# Create plot 
pcoa_plot <- pcoa_df%>%
  mutate(treat = factor(treat, 
                        levels = c("Untreated",
                                   "Colistin",
                                   "Gentamycin",
                                   "Kanamycin",
                                   "Metronidazole",
                                   "Vancomycin",
                                   "All"))) %>%
  ggplot(aes(axis_1, axis_2, color = treat, fill = treat, group = treat)) +
  geom_point(size = 2) +
  stat_ellipse() + 
  theme_classic() +
  ggtitle("Bray-Curtis dissimilarity") +
  labs(x = paste("PCoA Axis 1 (", formatC(pcoa$values$Relative_eig[1] * 100, digits = 2, format = "f"),"%)"),
       y = paste("PCoA Axis 2 (", formatC(pcoa$values$Relative_eig[2] * 100, digits = 2, format = "f"),"%)")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = cbpalette) +
  scale_color_manual(values = cbpalette) +
  guides(color = guide_legend(title = "Treatment"),
         fill = guide_legend(title = "Treatment")) +
  annotate(label = adonis_annot, geom = "text", y = 0.6, x = 0, size = 4)






# Pairwise dissimilarities ----

## Data preparation ----

# Extract distance matrix as a matrix
dist_mat <- dist%>% as.matrix()

# Remove the upper triangle to get rid of duplicates
dist_mat[upper.tri(dist_mat)] <- NA

# Melt to data frame and add metadata for each variable
library(reshape2)
dist_df <- dist_mat%>%
  melt(na.rm = T) %>%
  rename(sample1 = Var1, sample2 = Var2, distance = value) %>%
  left_join(metadata, by = c("sample1" = "sample"))%>%
  dplyr::select(-c(collection_date, source, sample_type, pool, pi, exp_number, sample_number, sample_description, primer)) %>%
  rename(treat1 = treat, mouse1 = source_id) %>%
  left_join(metadata, by = c("sample2" = "sample")) %>%
  dplyr::select(-c(collection_date, source, sample_type, pool, pi, exp_number, sample_number, sample_description, primer)) %>%
  rename(treat2 = treat, mouse2 = source_id) %>%
  filter(treat1 == treat2, mouse1 != mouse2) %>%
  dplyr::select(-treat2) %>%
  rename(treat = treat1)

## Analysis ----

# Run model
library(lme4)
dist_lm <- lmer(distance ~ treat + (1|mouse1) + (1|mouse2), data = dist_df) 

# Get p value of main effect
library(car)
dist_anova <- Anova(dist_lm)

# Create label with test statistic and p value for annotating plot
dist_annot <- paste0("Chisq = ", formatC(dist_anova$Chisq, digits = 3),
                     " p = ", formatC(dist_anova$'Pr(>Chisq)', digits = 3))

# Calculate means for plotting
library(emmeans)
dist_means <- emmeans(dist_lm, "treat") %>%
  data.frame() %>%
  rename(distance = emmean)


## Dissimilarities plot ----
dist_plot <- dist_df %>%
  mutate(treat = factor(treat, 
                        levels = c("Untreated",
                                   "Colistin",
                                   "Gentamycin",
                                   "Kanamycin",
                                   "Metronidazole",
                                   "Vancomycin",
                                   "All"))) %>%
  ggplot(aes(x = treat, y = distance, color = treat)) +
  geom_jitter(alpha = 0.2, size = 2, width = 0.1) +
  geom_point(data = dist_means, size = 4) +
  geom_errorbar(data = dist_means, aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) + 
  theme_classic() +
  scale_color_manual(values = cbpalette) +
  labs(x = "Treatment", y = "Distance ± 95% CI") +
  guides(color = guide_legend(title = "Treatment")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13, angle = 90),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  annotate(geom = "text", label = dist_annot, x = 4, y = 0.9, size = 4.5)



