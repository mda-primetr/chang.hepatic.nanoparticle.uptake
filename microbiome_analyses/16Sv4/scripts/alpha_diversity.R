# Load data as phyloseq object
library(dplyr)
library(phyloseq)

physeq <- readRDS("data/physeq_rarefied.RDS") 

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame() 

#Calculate alpha diversities
library(tibble)
alpha <- estimate_richness(physeq, measures = c("Observed", "InvSimpson", "Shannon"))%>%
  rownames_to_column("sample")%>%
  rename(Richness = "Observed") %>%
  right_join(metadata, by = "sample")




# Analyses ----

##Richness
rich_lm <- glm(Richness ~ treat, data = alpha, family = "poisson")
rich_anova <- anova(rich_lm)

## InvSimpson index
simp_lm <- lm(log10(InvSimpson) ~ treat, data = alpha)
simp_anova <- anova(simp_lm)

## Shannon index
shan_lm <- lm(Shannon ~ treat, data = alpha)
shan_anova <- anova(shan_lm, test = "F")


## Create data frame with test statistics for plot labels

alpha_annot_df <- data.frame(metric = c("Richness", "InvSimpson", "Shannon"),
                             label = c(paste0("Chisq = ", formatC(rich_anova$Deviance[2], digits = 3),
                                              " p = ", formatC(rich_anova$`Pr(>Chi)`[2], digits = 3)),
                                       paste0("F = ", formatC(simp_anova$`F value`[1], digits = 3),
                                              " p = ", formatC(simp_anova$`Pr(>F)`[1], digits = 3)),
                                       paste0("F = ", formatC(shan_anova$`F value`[1], digits = 3),
                                              " p = ", formatC(shan_anova$`Pr(>F)`[1], digits = 3))))

## Create data frames of estimated marginal means for each metric
library(emmeans)
rich_means <- data.frame(emmeans(rich_lm, "treat", type = "response")) %>%
  mutate(metric = "Richness") %>%
  rename(emmean = rate, 
         lower.CL = asymp.LCL,
         upper.CL = asymp.UCL)
simp_means <- emmeans(simp_lm, "treat", type = "response")%>%
  data.frame()%>%
  mutate(metric = "InvSimpson") %>%
  rename(emmean = "response")
shan_means <- emmeans(shan_lm, "treat", type = "response") %>%
  data.frame()%>%
  mutate(metric = "Shannon")

##Combine into one emmeans object
emmeans_df <- rbind(rich_means, simp_means, shan_means)




# Plots ----

## Create object of font sizes for plots
alpha_textsize = theme(axis.title = element_text(size = 15),
                       axis.text.y = element_text(size = 13),
                       axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5),
                       legend.position = "none")

##set colorblind friendly palette
cbpalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#F781BF", "#A65628" )

## Load custom function for plotting alpha diversities
source("src/custom_functions.R")


## Plot Richness
rich_plot <- plot_alpha(df = alpha, response = "Richness", x_axis = "treat",
                        point = "no", color_var = "treat", est_mean = emmeans_df,
                        annot = alpha_annot_df, annot_size = 3) +
  ylab("Mean Richness ± 95% CI")

## Plot InvSimpson index
simp_plot <- plot_alpha(df = alpha, response = "InvSimpson", x_axis = "treat",
                        point = "no", color_var = "treat", est_mean = emmeans_df,
                        annot = alpha_annot_df, annot_size = 3) +
  ylab("Mean InvSimpson ± 95% CI")


## Plot Shannon index
shan_plot <- plot_alpha(df = alpha, response = "Shannon", x_axis = "treat",
                        point = "no", color_var = "treat", est_mean = emmeans_df,
                        annot = alpha_annot_df, annot_size = 3) +
  ylab("Mean Shannon ± 95% CI")

# Combine the three plots into final figure
library(ggpubr)
alpha_combined_plot <- ggarrange(rich_plot, simp_plot, shan_plot, nrow = 1)
