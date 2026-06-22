# Load data as phyloseq object
library(dplyr)
library(phyloseq)
library(microViz)
physeq <- readRDS("data/processed_data/physeq_rarefied.RDS") 

# Create composition plot
library(ggplot2)
library(stringr)

# Save phyloseq object as figure_1I data
saveRDS(physeq,
        "data/figure_data/figure_1l_data.rds")
top_taxa_plot <- physeq%>%
  comp_barplot(tax_level = "Genus",
               taxon_renamer = function(x) str_replace_all(x, c( "_" = " ")),
               sample_order = "asis",
               n_taxa = 20,
               label = "source_id") +
  coord_flip() +
  ggtitle("Top taxa") +
  guides(fill = guide_legend(ncol = 1)) +
  xlab("Mouse") +
  facet_wrap(~treat, ncol = 1, scales = "free_y")


# View plot
top_taxa_plot

