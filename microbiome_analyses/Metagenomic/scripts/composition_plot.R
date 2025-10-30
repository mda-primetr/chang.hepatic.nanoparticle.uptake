library(phyloseq)
library(dplyr)

# Load data as phyloseq object
physeq <- readRDS("data/physeq_rarefied.rds")


# Generate top taxa plot
library(microViz)
library(ggplot2)
library(stringr)
top_taxa_plot <- physeq%>%
  comp_barplot(tax_level = "Species",
               taxon_renamer = function(x) str_replace_all(x, c("s__" = "", 
                                                                "g__" = "",
                                                                "f__" = "",
                                                                "p__" = "",
                                                                "_" = " ")),
               sample_order = "asis",
               n_taxa = 20,
               label = "mouse") +
  coord_flip() +
  ggtitle("Top taxa") +
  guides(fill = guide_legend(ncol = 1)) +
  xlab("Mouse") +
  facet_wrap(~treatment*timepoint, ncol = 2, scales = "free_y")

# Print plot
top_taxa_plot
