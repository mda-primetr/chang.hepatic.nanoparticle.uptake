library(phyloseq)
library(dplyr)

# Load data as phyloseq object
physeq <- readRDS("data/processed_data/physeq_rarefied.rds")

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame()


# Agglomerate phyloseq object to Species level
physeq_glom <- physeq%>%
  tax_glom(taxrank = "Species")


#Create empty lists to insert volcano plots and ancombc results tables
volc_list <- list()
ancombc_res_list <- list()


##For loop to run differential abundance analysis comparing each treatment to control
library(ANCOMBC)
library(forcats)
library(ggrepel)
for(trt in levels(metadata$treatment)) {
  
  #Filter physeq_glom to just include control and selected treatments
  physeq2 <- physeq_glom%>%
    subset_samples(treatment == trt)
  
  #Run analysis
  ancombc_out<- ancombc(physeq2, tax_level = "Species", formula = "timepoint")
  
  # Dynamically create the object name for ancombc_out_trt
  assign(paste0("ancombc_out_", trt), ancombc_out)
  
  #Create df of the desired ancombc results, including some wrangling
  ancombc_res_df <- ancombc_out$res$q_val%>%
    mutate(treatment = trt) %>%
    rename(qval = timepointafter)%>%
    left_join(rename(ancombc_out$res$lfc, lfc = timepointafter), by = "taxon") %>%
    left_join(physeq2 %>%
                tax_table() %>%
                data.frame() %>%
                rownames_to_column("taxon"),
              by = "taxon") %>%
    rename(species = Species)
  
  #Create labels for significant species for graphing
  ancombc_res_df <- ancombc_res_df%>%
    mutate(label = case_when(qval >= 0.05 ~ NA,
                             qval < 0.05 ~ species))
  
  #Format the labels
  ancombc_res_df$label <- gsub("s__", " ", ancombc_res_df$label)
  ancombc_res_df$label <- gsub("g__", " ", ancombc_res_df$label)
  ancombc_res_df$label <- gsub("f__", " ", ancombc_res_df$label)
  ancombc_res_df$label <- gsub("p__", " ", ancombc_res_df$label)
  ancombc_res_df$label <- gsub("_", " ", ancombc_res_df$label)
  
  #Add logchange direction values for plot legend
  ancombc_res_df <- ancombc_res_df%>%
    mutate(direction = case_when(lfc < 0 ~ "Before",
                                 lfc > 0 ~ "After"))%>%
    mutate(direction = case_when(qval >= 0.05 ~ "No change",
                                 qval < 0.05 ~ direction)) %>%
    mutate(direction = fct_relevel(direction, c("Before", "After", "No change")))
  
  # #Extract the position of the treatment in the vector of treatment levels for assigning colors
  # # Your vector of levels
  # levels_vector <-levels(metadata$treatment)
  # position <- which(levels_vector == trt)
  
  #Volcano plot
  volcano <- ancombc_res_df%>%
    ggplot(aes(x = lfc, y = -log10(qval), label = label, color = direction)) +
    geom_point(size = 1) +
    geom_hline(yintercept = -log10(0.05), col = "red", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    #geom_text_repel(point.size = 2, max.overlaps = Inf, size = 4, force = 30) + #Will add these individually for best formatting
    theme_classic() +
    scale_color_manual(values = c("Before" = "#E69F00" , "After" = "#009E73", "No change" = "grey")) +
    guides(color = guide_legend(title = "Timepoint")) +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    ggtitle(trt)
  
  
  # Dynamically create the plot object name and ancombc results table
  assign(paste0("volcano_", trt), volcano) 
  assign(paste0("ancombc_res_", trt), ancombc_res_df)
  
  volc_list[[trt]] <-  assign(paste0("volcano_", trt), volcano) 
  ancombc_res_list[[trt]] <- assign(paste0("ancombc_res_", trt), ancombc_res_df)
}

## Add labels

# volc_list$control <- volc_list$control +
#   geom_text_repel(point.size = 2, 
#                   max.overlaps = Inf, 
#                   show.legend = F,
#                   segment.alpha = 0.5,
#                   point.size = 3,
#                   force = 50,
#                   min.segment.length = 2)
# 
# volc_list$doxil <- volc_list$doxil +
#   geom_text_repel(point.size = 2, 
#                   force = 60, 
#                   show.legend = F,
#                   segment.alpha = 0.5,
#                   nudge_y = -0.1) +
#   ylim(-20, 105)
# 
# volc_list$metro <- volc_list$metro +
#   geom_text_repel(show.legend = F,
#                   segment.alpha = 0.5)
# 
# volc_list$'metro+doxil' <- volc_list$`metro+doxil` +
#   geom_text_repel(point.size = 2, 
#                   show.legend = F,
#                   segment.alpha = 0.5) +
#   ylim(-0.5, 3) +
#   xlim(-7, 15)


# Save data files
for(trt in levels(metadata$treatment)) {
  
  # Extract appropriate data frame from list
  df <- ancombc_res_list[[trt]]
  
  # Set filename
  filename <- paste0("data/figure_data/ext_data_fig10_", trt, "_data.csv")
  # save
  write.csv(df,
            filename,
            row.names = FALSE)
  
}
## Add species labels to plots
library(ggrepel)
volc_list$control <- volc_list$control +
  geom_text_repel(point.size = 2, 
                  max.overlaps = Inf, 
                  show.legend = F,
                  segment.alpha = 0.5,
                  point.size = 3,
                  force = 50,
                  min.segment.length = 2)

volc_list$doxil <- volc_list$doxil +
  geom_text_repel(point.size = 2, 
                  force = 60, 
                  show.legend = F,
                  segment.alpha = 0.5,
                  nudge_y = -0.1) +
  ylim(-20, 105)

volc_list$metro <- volc_list$metro +
  geom_text_repel(show.legend = F,
                  segment.alpha = 0.5)

volc_list$'metro+doxil' <- volc_list$`metro+doxil` +
  geom_text_repel(point.size = 2, 
                  show.legend = F,
                  segment.alpha = 0.5) +
  ylim(-0.5, 3) +
  xlim(-7, 15)


# Print plots
volc_list$control
volc_list$doxil
volc_list$metro
volc_list$'metro+doxil'
