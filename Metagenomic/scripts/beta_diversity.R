library(phyloseq)
library(dplyr)

# Load data as phyloseq object
physeq <- readRDS("data/processed_data/physeq_rarefied.rds")

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame()



#Create Brays Curtis dissimilarity matrix
library(vegan)
dist <- vegdist(t(as.matrix(otu_table(physeq))))


# Calculate PCOA and permanova ----
library(ape)
#Conduct PCoA on the dissimilarities
pcoa <- ape::pcoa(dist)

#Merge first two pcoa axes with metadata
pcoa_df <- data.frame(sample = row.names(pcoa$vectors),
                         axis_1 = pcoa$vectors[,1],
                         axis_2 = pcoa$vectors[,2])%>%
  right_join(metadata, "sample")


#Add points for centroids for plotting
pcoa_df <- pcoa_df%>%
  group_by(treatment, timepoint)%>%
  mutate(centroid.x = mean(axis_1),
         centroid.y = mean(axis_2))%>%
  ungroup()

# Run permanova
adonis_out <- adonis2(dist ~ treatment * timepoint +mouse,
                         pcoa_df,
                         by = "terms")

adonis_annot <- paste("Adonis2 test: treatment * timepoint R2 =", formatC(adonis_out$R2[3], digits = 2), 
                      "p =", adonis_out$`Pr(>F)`[3] )



#PCOA plots ----

# Save pcoa_df data
write.csv(pcoa_df,
          "data/figure_data/figure_4i_4j_data.csv",
          row.names = FALSE)


# Set colorblind friendly palette
treat_pal <-   c( '#737373', '#EE6677',  '#4477AA', '#AA3377')


##  Create plot with ellipses ----
pcoa_ellipse <- pcoa_df%>%
  ggplot(aes(axis_1, axis_2, color = treatment, shape = timepoint, fill = treatment, group = timepoint)) +
  geom_point(size = 2) +
  stat_ellipse(aes(linetype = timepoint)) + 
  theme_classic() +
  ggtitle("Bray Curtis dissimilarity") +
  labs(x = paste("PCoA Axis 1 (", formatC(pcoa$values$Relative_eig[1] * 100, digits = 2, format = "f"),"%)"),
       y = paste("PCoA Axis 2 (", formatC(pcoa$values$Relative_eig[2] * 100, digits = 2, format = "f"),"%)")) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) +
  scale_fill_manual(values = treat_pal) +
  scale_color_manual(values = treat_pal) +
  facet_wrap(~ treatment) +
  guides(color = guide_legend(title = "Treatment"),
         fill = guide_legend(title = "Treatment"),
         linetype = guide_legend(title = "Timepoint"),
         shape = guide_legend(title = "Timepoint")) +
  #annotate(label = adonis_annot, geom = "text", y = 0.5, x = 0, size = 5) +
  scale_shape_manual(values = c(1, 16)) +
  scale_linetype_manual(values = c(2, 1))


# Print plot
pcoa_ellipse


## Create plot showing lines across timepoints ----
pcoa_line <- pcoa_df%>%
  ggplot(aes(axis_1, axis_2, color = treatment, shape = timepoint, fill = treatment, group = timepoint)) +
  geom_point(size = 2) +
  #stat_ellipse(aes(linetype = timepoint)) + 
  geom_line(aes(group = mouse), alpha = 0.5) + 
  theme_classic() +
  labs(x = paste("PCoA Axis 1 (", formatC(pcoa$values$Relative_eig[1] * 100, digits = 2, format = "f"),"%)"),
       y = paste("PCoA Axis 2 (", formatC(pcoa$values$Relative_eig[2] * 100, digits = 2, format = "f"),"%)")) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) +
  scale_fill_manual(values = treat_pal) +
  scale_color_manual(values = treat_pal) +
  facet_wrap(~ treatment) +
  guides(color = guide_legend(title = "Treatment"),
         fill = guide_legend(title = "Treatment"),
         linetype = guide_legend(title = "Timepoint"),
         shape = guide_legend(title = "Timepoint")) +
  #annotate(label = adonis_annot, geom = "text", y = 0.5, x = 0, size = 5) +
  scale_shape_manual(values = c(1, 16)) +
  scale_linetype_manual(values = c(2, 1))

# Print plot
pcoa_line





# Pairwise dissimilarities ----

## Distances across timepoints ----
### Prepare data ----
dist_mat <- dist%>% as.matrix()

## Remove the upper triangle to get rid of duplicates
dist_mat[upper.tri(dist_mat)] <- NA

## Melt to data frame and add metadata for each variable
library(reshape2)
dist_time_df <- dist_mat%>%
  melt(na.rm = T) %>%
  rename(sample1 = Var1, sample2 = Var2, distance = value) %>%
  left_join(metadata, by = c("sample1" = "sample"))%>%
  dplyr::select(-c(principal_investigator, protocol_number, 
                   project_name, sample_weight_mg, sample_source,
                   date_sent_to_prime_tr, sequencing_requred, core_sample_id)) %>%
  rename(treatment1 = treatment, mouse1 = mouse, timepoint1 = timepoint) %>%
  left_join(metadata, by = c("sample2" = "sample")) %>%
  dplyr::select(-c(principal_investigator, protocol_number, 
                   project_name, sample_weight_mg, sample_source,
                   date_sent_to_prime_tr, sequencing_requred)) %>%
  rename(treatment2 = treatment, mouse2 = mouse, timepoint2 = timepoint) %>%
  filter(mouse1 == mouse2, timepoint1 != timepoint2) %>%
  dplyr::select(-c(treatment2, sample2, mouse2, timepoint2, core_sample_id)) %>%
  rename(treatment = treatment1, timepoint = timepoint1, mouse = mouse1)

### Run model using nonparametric approach ----

dist_time_stat <- kruskal.test(distance ~ treatment, data = dist_time_df)


### Create label with test statistic and p value for annotating plot

dist_time_annot <- paste0("Chisq = ", formatC(dist_time_stat$statistic, digits = 3),
                          ",p = ", formatC(dist_time_stat$p.value, digits = 3))


### Plot----

# Save distance time plot data
write.csv(dist_time_df,
          "data/figure_data/figure_4k_data.csv",
          row.names = FALSE)


dist_time_plot <- dist_time_df %>%
  ggplot(aes(x = treatment, y = distance, color = treatment)) +
  geom_boxplot() +
  # geom_point(alpha = 0.3, size = 3) +
  # geom_point(data = dist_time_means, size = 5) +
  # geom_errorbar(data = dist_time_means, aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) + 
  theme_classic() +
  scale_color_manual(values = treat_pal) +
  labs(x = "Treatment", y = "Distance ± 95% CI") +
  guides(color = guide_legend(title = "Treatment")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13, angle = 90),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  annotate(geom = "text", label = dist_time_annot, x = 2, y = 0.9, size = 4) +
  ggtitle("Within-mouse pairwise distances between timepoints")

# Print plot
dist_time_plot


## Distances within group x timepoint ----

### Prepare data ----
## Melt to data frame and add metadata for each variable
dist_treat_df <- dist_mat%>%
  melt(na.rm = T) %>%
  rename(sample1 = Var1, sample2 = Var2, distance = value) %>%
  left_join(metadata, by = c("sample1" = "sample"))%>%
  dplyr::select(-c(principal_investigator, protocol_number, 
                   project_name, sample_weight_mg, sample_source,
                   date_sent_to_prime_tr, sequencing_requred, core_sample_id)) %>%
  rename(treatment1 = treatment, mouse1 = mouse, timepoint1 = timepoint) %>%
  left_join(metadata, by = c("sample2" = "sample")) %>%
  dplyr::select(-c(principal_investigator, protocol_number, 
                   project_name, sample_weight_mg, sample_source,
                   date_sent_to_prime_tr, sequencing_requred, core_sample_id)) %>%
  rename(treatment2 = treatment, mouse2 = mouse, timepoint2 = timepoint) %>%
  filter(timepoint1 == timepoint2 & treatment1 == treatment2 & mouse1 != mouse2) %>%
  dplyr::select(-c(timepoint2, treatment2)) %>%
  rename(treatment = treatment1, timepoint = timepoint1, sample = sample1)


### Run model. Data allow parametric approach ----

dist_treat_lm <- lmer(log(distance) ~ treatment*timepoint + (1|mouse1) + (1|mouse2), data = dist_treat_df) 

### Get p values
dist_treat_anova <- Anova(dist_treat_lm)

### Calculate means for plotting

dist_treat_means <- emmeans(dist_treat_lm, pairwise ~ treatment + timepoint, type = "response")$emmeans %>%
  data.frame() %>%
  rename(distance = response)

### Generate plot

# Save treatment distance plot
write.csv(dist_treat_df,
          "data/figure_data/figure_4l_data.csv",
          row.names = FALSE)


dist_treat_plot <- dist_treat_df %>%
  ggplot(aes(x = timepoint, y = distance, color = treatment, shape = timepoint)) +
  geom_jitter(alpha = 0.2, 
              size = 2,
              width = 0.2) +
  geom_point(data = dist_treat_means, size = 5) +
  geom_errorbar(data = dist_treat_means, 
                aes(ymin =lower.CL, 
                    ymax = upper.CL,
                    linetype = timepoint), 
                width = 0.2) + 
  theme_classic() +
  facet_wrap(vars(treatment), nrow = 1) +
  scale_color_manual(values = treat_pal) +
  labs(x = "Treatment", y = "Distance ± 95% CI") +
  guides(color = guide_legend(title = "Treatment")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13, angle = 90),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  #annotate(geom = "text", label = dist_treat_annot, x = 3, y = 0.9, size = 4.5) +
  ggtitle("Pairwise distances within treatment x timepoint combination") +
  scale_linetype_manual(values = c(2, 1)) +
  scale_shape_manual(values = c(1, 16)) 


# Print plot
dist_treat_plot
