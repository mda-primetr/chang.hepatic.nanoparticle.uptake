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
            by = "sample") 



library(rstatix)
library(ggpubr)

# Get mapped reads per sample pre and post filter
mapped_reads_prefilter <- rgi_raw %>%
  group_by(sample) %>%
  summarise(mapped_reads_pre = sum(all_mapped_reads), .groups = "drop")

mapped_reads_postfilter <- rgi_clean %>%
  group_by(sample) %>%
  summarise(mapped_reads_post = sum(all_mapped_reads), .groups = "drop")

# Load total reads
stats_files <- list.files("data/raw_data/rgi_output",
                          pattern = "overall_mapping_stats.txt",
                          full.names = TRUE)

total_reads <- map_dfr(stats_files, function(f) {
  sample <- basename(f) %>%
    str_remove("^Sample_") %>%
    str_remove("\\.overall_mapping_stats\\.txt")
  total <- read_lines(f) %>%
    str_subset("Total reads:") %>%
    str_extract("[0-9]+") %>%
    as.numeric()
  tibble(sample = sample, total_reads = total)
})

# Build reads summary
reads_summary <- total_reads %>%
  left_join(mapped_reads_prefilter, by = "sample") %>%
  left_join(mapped_reads_postfilter, by = "sample") %>%
  left_join(rgi_clean %>% distinct(sample, treatment, timepoint), by = "sample") %>%
  mutate(
    rate_pre  = (mapped_reads_pre  / total_reads) * 100,
    rate_post = (mapped_reads_post / total_reads) * 100
  ) %>%
  pivot_longer(cols = c(rate_pre, rate_post),
               names_to = "filter_stage",
               values_to = "mapping_rate") %>%
  mutate(filter_stage = factor(filter_stage,
                               levels = c("rate_pre", "rate_post"),
                               labels = c("Pre-filter", "Post-filter")))

# Stats
stats_mapping <- reads_summary %>%
  group_by(treatment, filter_stage) %>%
  wilcox_test(mapping_rate ~ timepoint, paired = TRUE) %>%
  add_significance() %>%
  rstatix::add_xy_position(x = "timepoint", dodge = 0.8) %>%
  mutate(y.position = y.position * 1.2)


# Save data for plot
write.csv(reads_summary, 
          "data/figure_data/si_figure_23a_data",
          row.names = FALSE)


# Plot
reads_summary %>%
  ggplot(aes(x = timepoint, y = mapping_rate, fill = filter_stage)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 1.5, alpha = 0.6) +
  stat_pvalue_manual(stats_mapping,
                     label = "{filter_stage}: P = {round(p, 3)}",
                     tip.length = 0.01,
                     step.increase = 0.08) +
  facet_wrap(~treatment) +
  labs(y = "% reads mapping to CARD", x = NULL, fill = NULL) +
  theme_classic()





# Create final rgi object and save
rgi_final <- rgi_clean %>%
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
write.csv(rgi_final,
          "data/processed_data/rgi_data_clean.csv",
          row.names = FALSE)






