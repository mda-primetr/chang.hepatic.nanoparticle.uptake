source("src/libraries.r")

# Load CARD data frame
rgi_df <- read.csv("data/processed_data/rgi_data_clean.csv") %>%
  mutate(mouse = as.factor(mouse),
         timepoint = factor(timepoint,
                            levels = c("before", "after"))) %>%
  arrange(mouse, timepoint)

# Load card phyloseq object
physeq_card <- readRDS("data/processed_data/physeq_card.rds")

# Extract card classification
card_classes <- physeq_card %>%
  tax_table() %>%
  data.frame()
# Prepare input — MaAsLin2 needs samples as rows, features as columns
maaslin_input <- rgi_df %>%
  dplyr::select(sample, aro_term, tpm) %>%
  pivot_wider(names_from = aro_term, values_from = tpm, values_fill = 0) %>%
  column_to_rownames("sample")

# Metadata
maaslin_meta <- rgi_df %>%
  dplyr::select(sample, mouse, treatment, timepoint) %>%
  distinct() %>%
  column_to_rownames("sample")


# Gene level ----

# Run maaslin comparing before vs after within each treatment
treatments <- unique(maaslin_meta$treatment)

maaslin_results_gene <- map_dfr(treatments, function(tx) {
  
  
  # Subset metadata and input data to this treatment
  meta_sub <- maaslin_meta %>% filter(treatment == tx)
  data_sub <- maaslin_input[rownames(maaslin_input) %in% rownames(meta_sub), ]
  
  # Run MaAsLin2
  fit <- Maaslin2(
    input_data     = data_sub,
    input_metadata = meta_sub,
    output         = paste0("output/maaslin/card/gene/", tx),
    fixed_effects  = "timepoint",
    random_effects = "mouse",
    normalization  = "NONE",
    transform      = "LOG",
    min_prevalence = 0.1,
    reference      = "timepoint,before"
  )
  
  # Read results and add treatment label
  read_tsv(paste0("output/maaslin/card/gene/", tx, "/all_results.tsv"),
           show_col_types = FALSE) %>%
    mutate(treatment = tx)
})


# Create lookup table
name_lookup <- card_classes %>%
  dplyr::select(aro_term, amr_gene_family, drug_class, resistance_mechanism) %>%
  mutate(aro_term_clean = make.names(aro_term) %>% str_remove("^X")) %>%
  as_tibble()

# Clean maaslin results and join
maaslin_results_gene_clean <- maaslin_results_gene %>%
  mutate(aro_term_clean = gsub("^X", "", feature)) %>%
  left_join(name_lookup, by = "aro_term_clean") %>%
  mutate(log2fc = log2(exp(coef))) %>%
  dplyr::select(-feature, -aro_term_clean)

# Save this
write.csv(maaslin_results_gene_clean %>%
            arrange(qval),
          "output/maaslin/card/gene/maaslin_results_gene.csv"
          )


# Drug class level ----

# Aggregate to drug_class level with separate_rows
maaslin_input_drug_class <- rgi_df %>%
  separate_rows(drug_class, sep = "; ") %>%
  group_by(sample, drug_class) %>%
  summarise(tpm = sum(tpm), .groups = "drop") %>%
  pivot_wider(names_from = drug_class, values_from = tpm, values_fill = 0) %>%
  column_to_rownames("sample")


maaslin_results_drug_class <- map_dfr(treatments, function(tx) {
  meta_sub <- maaslin_meta %>% filter(treatment == tx)
  data_sub <- maaslin_input_drug_class[rownames(maaslin_input_drug_class) %in% rownames(meta_sub), ]
  fit <- Maaslin2(
    input_data     = data_sub,
    input_metadata = meta_sub,
    output         = paste0("output/maaslin/card/drug_class/", tx),
    fixed_effects  = "timepoint",
    random_effects = "mouse",
    normalization  = "NONE",
    transform      = "LOG",
    min_prevalence = 0.1,
    reference      = "timepoint,before"
  ) 
  read_tsv(paste0("output/maaslin/card/drug_class/", tx, "/all_results.tsv"),
           show_col_types = FALSE) %>%
    mutate(treatment = tx)
}) %>%
  mutate(log2fc = log2(exp(coef)),
          feature = gsub("^X", "", feature),
         feature = str_replace_all(feature, "\\.", " ")) %>%
  rename(drug_class = feature)


# Save this
write.csv(maaslin_results_drug_class %>%
            arrange(qval),
          "output/maaslin/card/drug_class/maaslin_results_drug_class.csv"
)



# Resistance mechanism level ----
maaslin_input_mechanism <- rgi_df %>%
  separate_rows(resistance_mechanism, sep = "; ") %>%
  group_by(sample, resistance_mechanism) %>%
  summarise(tpm = sum(tpm), .groups = "drop") %>%
  pivot_wider(names_from = resistance_mechanism, values_from = tpm, values_fill = 0) %>%
  column_to_rownames("sample")


maaslin_results_mechanism <- map_dfr(treatments, function(tx) {
  meta_sub <- maaslin_meta %>% filter(treatment == tx)
  data_sub <- maaslin_input_mechanism[rownames(maaslin_input_mechanism) %in% rownames(meta_sub), ]
  fit <- Maaslin2(
    input_data     = data_sub,
    input_metadata = meta_sub,
    output         = paste0("output/maaslin/card/resistance_mechanism/", tx),
    fixed_effects  = "timepoint",
    random_effects = "mouse",
    normalization  = "NONE",
    transform      = "LOG",
    min_prevalence = 0.1,
    reference      = "timepoint,before"
  )
  read_tsv(paste0("output/maaslin/card/resistance_mechanism/", tx, "/all_results.tsv"),
           show_col_types = FALSE) %>%
    mutate(treatment = tx)
}) %>%
  mutate(log2fc = log2(exp(coef)),
         feature = gsub("^X", "", feature),
         feature = str_replace_all(feature, "\\.", " ")) %>%
  rename(resistance_mechanism = feature)


# Save this
write.csv(maaslin_results_mechanism%>%
            arrange(qval),
          "output/maaslin/card/resistance_mechanism/maaslin_results_resistance_mechanism.csv"
)




# Plots ----

# Prep function — adds label and direction columns
prep_volcano <- function(data, feature_col, label_col) {
  data %>%
    mutate(
      direction = case_when(log2fc < 0 ~ "before",
                            log2fc > 0 ~ "after"),
      point_color = case_when(qval < 0.05 & log2fc < 0 ~ "before",
                              qval < 0.05 & log2fc > 0 ~ "after",
                              TRUE ~ "ns"),
      point_color = factor(point_color,
                           levels = c('before',
                                      'after',
                                      'ns')),
      label = ifelse(qval < 0.05,
                     paste0(as.character(.data[[label_col]]), " (qval = ", round(qval, digits = 4), ")"),
                     NA),
    )
}

# Updated label for gene level — include drug class and resistance mechanism
maaslin_results_gene_clean <- maaslin_results_gene_clean %>%
  mutate(gene_label = ifelse(qval < 0.05,
                             paste0(aro_term, "\n", drug_class, "\n", resistance_mechanism),
                             NA))

# Updated volcano function using point_color
make_volcano <- function(data, feature_col, title) {
  data %>%
    ggplot(aes(x = log2fc, y = -log10(pval),
               label = label, color = point_color)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), col = "red", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_text_repel(point.size = 2,
                    max.overlaps = Inf,
                    size = 4,
                    force = 50,
                    segment.alpha = 0.3,
                    show.legend = FALSE) +
    theme_classic() +
    guides(color = guide_legend("timepoint")) +
    scale_color_manual(values = time_pal) +
    ggtitle(title)
}

# Gene level plots — use gene_label instead of aro_term
volc_gene <- maaslin_results_gene_clean %>%
  prep_volcano(feature_col = "aro_term", label_col = "gene_label") %>%
  group_by(treatment) %>%
  group_map(~ make_volcano(.x, "aro_term",
                           paste0("AMR genes - ", .y$treatment)),
            .keep = TRUE) %>%
  setNames(treatments)

volc_gene
# Drug class and mechanism plots unchanged
volc_drug_class <- maaslin_results_drug_class %>%
  prep_volcano(feature_col = "drug_class", label_col = "drug_class") %>%
  group_by(treatment) %>%
  group_map(~ make_volcano(.x, "drug_class",
                           paste0("Drug class - ", .y$treatment)),
            .keep = TRUE) %>%
  setNames(treatments)

volc_drug_class
volc_mechanism <- maaslin_results_mechanism %>%
  prep_volcano(feature_col = "resistance_mechanism", label_col = "resistance_mechanism") %>%
  group_by(treatment) %>%
  group_map(~ make_volcano(.x, "resistance_mechanism",
                           paste0("Resistance mechanism - ", .y$treatment)),
            .keep = TRUE) %>%
  setNames(treatments)
volc_mechanism



# Create directories
dir.create("output/figures/card/differential_abundance/gene", recursive = TRUE)
dir.create("output/figures/card/differential_abundance/drug_class", recursive = TRUE)
dir.create("output/figures/card/differential_abundance/resistance_mechanism", recursive = TRUE)

# Save gene level plots
walk(treatments, function(tx) {
  ggsave(
    filename = paste0("output/figures/card/differential_abundance/gene/volc_gene_", tx, ".pdf"),
    plot = volc_gene[[tx]],
    width = 8, height = 6, dpi = 400
  )
})

# Save drug class level plots
walk(treatments, function(tx) {
  ggsave(
    filename = paste0("output/figures/card/differential_abundance/drug_class/volc_drug_class_", tx, ".pdf"),
    plot = volc_drug_class[[tx]],
    width = 8, height = 6, dpi = 400
  )
})

# Save resistance mechanism level plots
walk(treatments, function(tx) {
  ggsave(
    filename = paste0("output/figures/card/differential_abundance/resistance_mechanism/volc_mechanism_", tx, ".pdf"),
    plot = volc_mechanism[[tx]],
    width = 8, height = 6, dpi = 400
  )
})
