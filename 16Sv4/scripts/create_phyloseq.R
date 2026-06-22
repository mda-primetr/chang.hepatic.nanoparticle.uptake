library(rbiom)

# Load .biom file with feature counts

biom <- as_rbiom("data/raw_data/wjiang.liver-macrophages.batch2.biom")

# Extract feature table and transform for phyloseq creation
library(phyloseq)
feature_table <- biom$counts %>%
  as.matrix() %>%
  otu_table(taxa_are_rows = TRUE)
# Extract taxonomy table and transform for phyloseq creation
library(tibble)
tax <- biom$taxonomy %>%
  data.frame() %>%
  column_to_rownames(".otu") %>%
  as.matrix() %>%
  tax_table()

# Extract 
# Load metadata and transform for phyloseq creation

metadata <- read.csv("data/metadata/16sv4_metadata.csv") %>%
  mutate(
    treat = factor(treat, 
                   levels = c("Untreated", "Colistin", "Gentamycin", "Kanamycin", "Metronidazole", "Vancomycin", "All")),
    rownames = sample) %>% # Duplicate sample column to set sample values as row names while preserving the column
  column_to_rownames("rownames") %>%
  sample_data()

# Extract phylogenetic tree
phylogeny <- biom$tree

# Create phyloseq object
physeq_unrarefied <- phyloseq(feature_table,
                              tax,
                              metadata,
                              phylogeny)

# Rarefy phyloseq object to sample with fewest counts
physeq_rarefied <- rarefy_even_depth(physeq_unrarefied, rngseed = 9281)


# Save phyloseq object
saveRDS(physeq_rarefied,
        "data/processed_data/physeq_rarefied.rds")
