source("src/libraries.r")

# Load phyloseq object
physeq <- readRDS("data/processed_data/physeq_rarefied.RDS")

#Extract metadata

metadata <- sample_data(physeq) %>%
  data.frame() %>%
  rownames_to_column("sample")
# Load RGI data


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


# Load total reads file
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

rgi_tpm <- rgi_raw %>%
  filter(average_percent_coverage >= 80) %>%
  mutate(rpk = all_mapped_reads / (reference_length / 1000)) %>%
  group_by(sample) %>%
  mutate(tpm = rpk / sum(rpk) * 1e6) %>%
  ungroup()

# Merge rgi file with metdata
rgi_final <- rgi_tpm %>%
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


# Save rgi file
write.csv(rgi_final,
          "data/processed_data/rgi_data_clean.csv",
          row.names = FALSE)

# Visualize read depth before and after filtering and 

# Pre-filtering mapped reads (all detections regardless of coverage)
mapped_reads_prefilter <- rgi_raw %>%
  group_by(sample) %>%
  summarise(mapped_reads_pre = sum(all_mapped_reads), .groups = "drop")

# Post-filtering mapped reads (≥80% coverage)
mapped_reads_postfilter <- rgi_tpm %>%
  group_by(sample) %>%
  summarise(mapped_reads_post = sum(all_mapped_reads), .groups = "drop")

# Join everything together
reads_summary <- total_reads %>%
  left_join(mapped_reads_prefilter, by = "sample") %>%
  left_join(mapped_reads_postfilter, by = "sample") %>%
  left_join(rgi_df %>% distinct(sample, treatment, timepoint), by = "sample") %>%
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

stats_mapping <- reads_summary %>%
  group_by(treatment, filter_stage) %>%
  wilcox_test(mapping_rate ~ timepoint, paired = TRUE) %>%
  add_significance() %>%
  add_xy_position(x = "timepoint", dodge = 0.8, scales = "free_y") %>%
  mutate(y.position = case_when(treatment == "control" | treatment == "doxil" ~ y.position + 0.02,
                                treatment == "metro+doxil" ~ y.position - 0.01,
                                .default = y.position))

card_read_summary_plot <- reads_summary %>%
  ggplot(aes(x = timepoint, y = mapping_rate, fill = filter_stage)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 1.5, alpha = 0.6) +
  stat_pvalue_manual(stats_mapping, 
                     label = "{filter_stage}: p = {round(p, 3)}", 
                     tip.length = 0.01,
                     step.increase = 0.08) +
  facet_wrap(~treatment) +
  labs(y = "% reads mapping to CARD", x = NULL, fill = NULL) +
  theme_classic()

card_read_summary_plot

ggsave("output/figures/card/read_summary/read_summary_plot.pdf",
       card_read_summary_plot,
       height = 6.09,
       width = 5.3)

# Make phyloseq object

# otu table
rgi_wide <- rgi_tpm %>%
  pivot_wider(id_cols = aro_term,
              names_from = sample,
              values_from = tpm,
              values_fill = 0) %>%
  column_to_rownames("aro_term")

# Tax table
rgi_class <- rgi_final %>%
  dplyr::select(aro_term,
                amr_gene_family,
                drug_class,
                resistance_mechanism) %>%
  distinct(aro_term,
           .keep_all = TRUE) %>%
  mutate(rownames = aro_term) %>%
  column_to_rownames("rownames") %>%
  as.matrix()

# Create phyloseq object
physeq_card <- phyloseq(
  sample_data(metadata %>% column_to_rownames("sample")),
  tax_table(rgi_class),
  otu_table(rgi_wide, taxa_are_rows =TRUE)
)

# Save CARD phyloseq object
saveRDS(physeq_card,
        "data/processed_data/physeq_card.rds")


