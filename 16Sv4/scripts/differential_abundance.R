library(phyloseq)
library(dplyr)

# Load data as phyloseq object
physeq <- readRDS("data/processed_data/physeq_rarefied.rds")

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame()

#Agglomerate data to genus level
physeq_glom <- physeq%>%
  tax_glom("Genus") 

# Extract taxonomic table
tax <- tax_table(physeq_glom) %>%
  data.frame() %>%
  rownames_to_column("asv")

#Create empty lists to insert volcano plots and ancombc results tables
volc_list <- list()
ancombc_res_list <- list()

#set colorblind friendly palette for plot s
cbpalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#F781BF", "#A65628" )

#For loop to run differential abundance analyses comparing each treatment to control and create plots
library(ANCOMBC)
library(forcats)

for(trt in levels(metadata$treat)[2:length(levels(metadata$treat))]) {
  
  #Filter physeq to just include control and selected treatments
  physeq2 <- physeq_glom%>%
    subset_samples(treat == "Untreated" | treat == trt)
  
  #Run analysis
  ancombc_out<- ancombc(physeq2, tax_level = "Genus", formula = "treat")
  
  # Dynamically create the object name for ancombc_out_trt
  assign(paste0("ancombc_out_", trt), ancombc_out)
  
  
  
  #Create df of the desired ancombc results, including some wrangling
  ancombc_res_df <- ancombc_out$res$q_val%>%
    mutate(treat = trt) %>%
    rename(qval = paste0("treat", trt))%>%
    left_join(rename(ancombc_out$res$lfc, lfc = paste0("treat", trt)), by = "taxon") %>%
    left_join(tax, by = c("taxon" = "asv")) %>%
    rename(genus = Genus)
  
  
  #Create labels for significant species for graphing
  ancombc_res_df <- ancombc_res_df%>%
    mutate(label = case_when(qval >= 0.05 ~ NA,
                             qval < 0.05 ~ genus))
  
  #Format the labels
  ancombc_res_df$label <- gsub("__", " ", ancombc_res_df$label)
  ancombc_res_df$label <- gsub("_", " ", ancombc_res_df$label)
  
  #Add logchange direction values for plot legend
  ancombc_res_df <- ancombc_res_df%>%
    mutate(direction = case_when(lfc < 0 ~ "Untreated",
                                 lfc > 0 ~ trt))%>%
    mutate(direction = case_when(qval >= 0.05 ~ "No change",
                                 qval < 0.05 ~ direction)) %>%
    mutate(direction = fct_relevel(direction, c("Untreated", trt, "No change")))
  
  #Extract the position of the treatment in the vector of treatment levels for assigning colors
  # Your vector of levels
  levels_vector <-levels(metadata$treat)
  position <- which(levels_vector == trt)
  
  #Volcano plot
  volcano <- ancombc_res_df%>%
    ggplot(aes(x = lfc, y = -log10(qval), label = label, color = direction)) +
    geom_point(size = 2) +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    #geom_text_repel(point.size = 2, max.overlaps = Inf, size = 4, force = 30) + #Will add these individually for best formatting
    theme_classic() +
    scale_color_manual(values = c(cbpalette[c(1,position)], "grey")) +
    guides(color = guide_legend(title = "Log change direction")) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  # Dynamically create the plot object name and ancombc results table
  assign(paste0("volcano_", trt), volcano) 
  assign(paste0("ancombc_res_", trt), ancombc_res_df)
  
  volc_list[[trt]] <-  assign(paste0("volcano_", trt), volcano) 
  ancombc_res_list[[trt]] <- assign(paste0("ancombc_res_", trt), ancombc_res_df)
}


# Add genus labels to plots
library(ggrepel)
volc_list$Colistin <- volc_list$Colistin +
  geom_text_repel(point.size = 2, max.overlaps = Inf, size = 3.5, force = 30)

volc_list$Gentamycin <- volc_list$Gentamycin +
  geom_text_repel(point.size = 2, max.overlaps = Inf, size = 3.5, force = 30) +
  annotate(x = -2.5, y = 190, geom = "text", label = "Bifidobacterium", color = cbpalette[1])  

volc_list$Kanamycin <- volc_list$Kanamycin +
  geom_text_repel(point.size = 2, max.overlaps = Inf, size = 3.5, force = 30) +
  ylim(-20, 290)

volc_list$Metronidazole<- volc_list$Metronidazole +
  geom_text_repel(point.size = 2, max.overlaps = Inf, size = 3.5, force = 30) 

volc_list$Vancomycin <- volc_list$Vancomycin +
  geom_text_repel(point.size = 2, max.overlaps = Inf, size = 3.5, force = 30) +
  ylim(-20,230)

volc_list$All <- volc_list$All +
  geom_text_repel(point.size = 2, max.overlaps = Inf, size = 3.5, force = 40) +
  annotate(x = 5.5, y = 300, geom = "text", label = "Dubosiella", color = cbpalette[7]) +
  ylim(-25,300)

volc_list
