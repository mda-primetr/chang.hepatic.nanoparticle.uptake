library(dplyr)
library(purrr)
library(stringr)
library(janitor)
library(readr)
library(phyloseq)

# Load phyloseq object
physeq <- readRDS("data/processed_data/physeq_rarefied.rds")

# Extract metadata
metadata <- physeq %>%
  sample_data() %>%
  data.frame()

# Load phyloseq object
# Get list of RGI sample file names in the folder
rgi_files <- list.files("data/raw_data/rgi_output", pattern = "gene_mapping_data.txt", full.names = TRUE)

# Load and compile rgi output files
rgi_raw <- map_dfr(rgi_files, function(f) {
  # Extract the sample name
  sample_id <- basename(f) %>%
    str_remove("\\.gene_mapping_data\\.txt") %>%
    str_remove("Sample_")
  # Import the file
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = sample_id) %>%
    clean_names()
}) 

# Plot distribution of percent coverage to determine quality cutoff
ggplot(rgi_raw, aes(x = average_percent_coverage)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 80, color = "red", linetype = "dashed") +
  labs(title = "Distribution of gene coverage",
       x = "Average percent coverage", y = "Count")

# Filter out reads with < 80% coverage
rgi_clean <- rgi_raw %>%
  filter(average_percent_coverage >= 80) %>%
  mutate(rpk = all_mapped_reads / (reference_length / 1000)) %>%
  group_by(sample) %>%
  mutate(tpm = rpk / sum(rpk) * 1e6) %>%
  ungroup() %>%
  left_join(metadata,
            by = "sample") %>%
  dplyr::select(sample, 
                mouse,
                treatment,
                timepoint,
                aro_term,
                amr_gene_family,
                drug_class,
                resistance_mechanism,
                tpm)



# Save this dataframe
write.csv(rgi_clean,
          "data/processed_data/rgi_data_clean.csv",
          row.names = FALSE)






