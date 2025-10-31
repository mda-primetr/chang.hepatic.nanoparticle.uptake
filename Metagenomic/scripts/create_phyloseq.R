library(dplyr)
library(readr)
library(tidyr)

# Load and transform metaphlan readstats
readstats <- read_tsv("data/raw_data/wjiang.liver-macrophages.batch1.metaphlan4.readStats.txt") %>%
  separate(col = c(clade), into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ID"), sep = "\\|") %>%
  mutate(ID = gsub("t__SGB", "", ID)) %>%
  mutate(ID = gsub("_group", "", ID)) %>%
  filter(!is.na(ID))

# Convert readstats to feature table
library(phyloseq)
feature_table <- readstats %>%
  dplyr::select(-c(Kingdom:Species)) %>%
  # dplyr::mutate(TaxID = paste0("TaxID", 1:nrow(t4) - 1)) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix() %>%
  otu_table(taxa_are_rows = T)

# Create tax table from readstats

library(tibble)
tax <- readstats %>%
  dplyr::select(Kingdom:ID) %>%
  dplyr::mutate(Kingdom = ifelse(Kingdom == "", "Unknown", Kingdom)) %>%
  rowwise() %>%
  dplyr::mutate(Phylum = case_when(
    grepl("p__PGB", Phylum) ~ last(c_across(Kingdom)[!is.na(c_across(Kingdom))]),
    TRUE ~ as.character(Phylum)
  )) %>%
  dplyr::mutate(Class = case_when(
    grepl("c__CFGB", Class) ~ last(c_across(Phylum)[!is.na(c_across(Phylum))]),
    TRUE ~ as.character(Class)
  )) %>%
  dplyr::mutate(Order = case_when(
    grepl("o__OFGB", Order) ~ last(c_across(Class)[!is.na(c_across(Class))]),
    TRUE ~ as.character(Order)
  )) %>%
  dplyr::mutate(Family = case_when(
    grepl("f__FGB", Family) ~ last(c_across(Order)[!is.na(c_across(Order))]),
    TRUE ~ as.character(Family)
  )) %>%
  dplyr::mutate(Genus = case_when(
    grepl("g__GGB", Genus) ~ last(c_across(Family)[!is.na(c_across(Family))]),
    TRUE ~ as.character(Genus)
  )) %>%
  dplyr::mutate(Species = case_when(
    grepl("s__GGB", Species) ~ last(c_across(Genus)[!is.na(c_across(Genus))]),
    TRUE ~ as.character(Species)
  )) %>%
  dplyr::mutate(Phylum = case_when(
    !grepl("p__", Phylum) ~ paste0("LKT_", Phylum),
    TRUE ~ as.character(Phylum)
  )) %>% # This is for labelling unknown taxa levels
  dplyr::mutate(Class = case_when(
    !grepl("c__", Class) ~ paste0("LKT_", Class),
    TRUE ~ as.character(Class)
  )) %>%
  dplyr::mutate(Order = case_when(
    !grepl("o__", Order) ~ paste0("LKT_", Order),
    TRUE ~ as.character(Order)
  )) %>%
  dplyr::mutate(Family = case_when(
    !grepl("f__", Family) ~ paste0("LKT_", Family),
    TRUE ~ as.character(Family)
  )) %>%
  dplyr::mutate(Genus = case_when(
    !grepl("g__", Genus) ~ paste0("LKT_", Genus),
    TRUE ~ as.character(Genus)
  )) %>%
  dplyr::mutate(Species = case_when(
    !grepl("s__", Species) ~ paste0("LKT_", Species),
    TRUE ~ as.character(Species)
  )) %>%
  ungroup() %>%
  column_to_rownames(var = "ID") %>%
  as.matrix() %>% # Convert to matrix for tax_table function
  tax_table() # Convert to tax table for phyloseq function



# Import metadata
metadata <- read.csv("data/metadata/metagenomic_metadata.csv") %>%
  mutate(rownames = sample) %>% # Duplicate sample column to assign to rownames while preserving column
  column_to_rownames("sample")%>% # Set sample values as rownames for phyloseq function
  sample_data() # Convert to sample_data object for phyloseq function

# Load phylogenetic tree
library(ape)
tree <- ape::read.tree("data/raw_data/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk")


# Create phyloseq object
physeq_unrare <- phyloseq(feature_table, 
                          tax, 
                          tree, 
                          metadata)

# Rarefy
physeq_rare <- rarefy_even_depth(physeq_unrare, rngseed = 9281)

# Save phyloseq object
saveRDS(physeq_rare, "data/processed_data/physeq_rarefied.rds")