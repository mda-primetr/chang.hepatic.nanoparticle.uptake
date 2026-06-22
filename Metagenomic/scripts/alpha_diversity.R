library(phyloseq)
library(dplyr)


# Load data as phyloseq object
physeq <- readRDS("data/processed_data/physeq_rarefied.rds")

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame()

library(tibble)
#Calculate alpha diversities
alpha <- estimate_richness(physeq, measures = c("Observed", "InvSimpson", "Shannon"))%>%
  rownames_to_column("sample")%>%
  rename(Richness = "Observed") %>%
  left_join(metadata, by = "sample")





# Analyze diversity metrics ----

## Richness ----

### Run glm
library(lme4)
rich_lm <- lmer(Richness ~ treatment * timepoint + (1|mouse), data = alpha)


# Extract anova table
library(car)
rich_anova <- Anova(rich_lm)




## InvSimpson ----

### Glm
simp_lm <- glmer(InvSimpson ~ treatment * timepoint + (1|mouse), family = "Gamma", data = alpha)


### create Anova table
simp_anova <- Anova(simp_lm)



## Shannon ----

### Shannon glm
shan_lm <- lmer(Shannon ~ treatment * timepoint + (1|mouse), data = alpha)

### Create anova table
shan_anova <- Anova(shan_lm)




## Create data frames with test statistics for plot tables

rich_anova_df <- rich_anova %>%
  data.frame() %>%
  mutate(across(where(is.numeric),  ~ formatC(.x, digits = 3)))

simp_anova_df <- simp_anova %>%
  data.frame() %>%
  mutate(across(where(is.numeric),  ~ formatC(.x, digits = 3)))

shan_anova_df <- shan_anova %>%
  data.frame() %>%
  mutate(across(where(is.numeric),  ~ formatC(.x, digits = 3)))



## Create data frames of estimated marginal means for each metric
library(emmeans)
rich_means <- data.frame(emmeans(rich_lm, pairwise ~ treatment * timepoint, type = "response")$emmeans) %>%
  rename(Richness = emmean)
simp_means <- emmeans(simp_lm, pairwise ~ treatment * timepoint, type = "response")$emmeans%>%
  data.frame()%>%
  rename(InvSimpson = response) 
shan_means <- emmeans(shan_lm, pairwise ~ treatment * timepoint, type = "response")$emmeans %>%
  data.frame()%>%
  rename(Shannon = emmean)




# Plots ----
## Create object of font sizes for plots
library(ggplot2)
alpha_textsize =   theme(axis.text.y = element_text(size = 14),
                         axis.title = element_text(size = 16),
                         axis.text.x = element_text(size = 12, angle = 45, hjust = 0.5, vjust = 0.5),
                         legend.title = element_text(size = 14),
                         strip.text = element_text(size = 12),
                         legend.position = "none") 

# Set colorblind friendly palette
treat_pal <-   c( '#737373', '#EE6677',  '#4477AA', '#AA3377')

# Save plot data
write.csv(alpha,
          "data/figure_data/figure_4f_4g_4h_data.csv",
          row.names = FALSE)
##Create plots

rich_plot <- alpha %>%
  ggplot(aes(x = timepoint, y = Richness, color = treatment, shape = timepoint)) +
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.2) +
  geom_line(aes(group = mouse), alpha = 0.2) +
  facet_wrap(vars(treatment), nrow = 1) +
  theme_classic() +
  scale_color_manual(values = treat_pal) +
  alpha_textsize +
  labs(x = "Timepoint") +
  geom_point(data = rich_means, position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(data = rich_means, 
                position = position_dodge(width = 0.5), 
                aes(ymin = lower.CL, ymax = upper.CL, linetype = timepoint),
                width = 0.4)+
  scale_linetype_manual(values = c(2, 1)) +
  scale_shape_manual(values = c(1, 16)) 

simp_plot <- alpha %>%
  ggplot(aes(x = timepoint, y = InvSimpson, color = treatment, shape = timepoint)) +
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.2) +
  geom_line(aes(group = mouse), alpha = 0.2) +
  facet_wrap(vars(treatment), nrow = 1) +
  theme_classic() +
  scale_color_manual(values = treat_pal) +
  alpha_textsize +
  labs(x = "Timepoint") +
  geom_point(data = simp_means, position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(data = simp_means, 
                position = position_dodge(width = 0.5), 
                aes(ymin = asymp.LCL, ymax = asymp.UCL, linetype = timepoint),
                width = 0.4) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_shape_manual(values = c(1, 16)) 

shan_plot <- alpha %>%
  ggplot(aes(x = timepoint, y = Shannon, color = treatment, shape = timepoint)) +
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.2) +
  geom_line(aes(group = mouse), alpha = 0.2) +
  facet_wrap(vars(treatment), nrow = 1) +
  theme_classic() +
  alpha_textsize +
  scale_color_manual(values = treat_pal) +
  labs(x = "Timepoint") +
  geom_point(data = shan_means, position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(data = shan_means, 
                position = position_dodge(width = 0.5), 
                aes(ymin = lower.CL, ymax = upper.CL, linetype = timepoint),
                width = 0.4)+
  scale_linetype_manual(values = c(2, 1)) +
  scale_shape_manual(values = c(1, 16)) 


# Print plots and anova tables
ggarrange(rich_plot,
          ggtexttable(rich_anova_df),
          ncol = 1,
          heights = c(2, 1)
          )
ggarrange(simp_plot,
          ggtexttable(simp_anova_df),
          ncol = 1,
          heights = c(2, 1)
)
ggarrange(shan_plot,
          ggtexttable(shan_anova_df),
          ncol = 1,
          heights = c(2, 1)
)
