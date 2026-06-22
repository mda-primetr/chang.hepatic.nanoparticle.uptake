source("src/libraries.r")

# Load CARD data frame
rgi_df <- read.csv("data/processed_data/rgi_data_clean.csv") %>%
  mutate(mouse = as.factor(mouse),
         timepoint = factor(timepoint,
                            levels = c("before", "after"))) %>%
  arrange(mouse, timepoint) %>%
  unite(mouse_time, mouse, timepoint) %>%
  mutate(mouse_time = fct_inorder(mouse_time))
# Load CARD phyloseq object
physeq_card <- readRDS("data/processed_data/physeq_card.rds")



# By aro_term ----



comp_plot_card_aro <- rgi_df %>%
  ggplot(aes(x = mouse_time, y = tpm, fill = aro_term)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~treatment, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = ggsci::pal_d3("category20")(20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

comp_plot_card_aro
ggsave("output/figures/composition_plots/card/aro_comp.pdf",
       comp_plot_card_aro,
       height = 8,
       width = 9)


# By gene family ----

# Sum aro's in each gene family
gene_family_df <- rgi_df %>%
  group_by(mouse_time, treatment, amr_gene_family) %>%
  summarise(tpm = sum(tpm), .groups = "drop")

comp_plot_card_gene_family <- gene_family_df %>%
  ggplot(aes(x = mouse_time, y = tpm, fill = amr_gene_family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~treatment , ncol = 2, scales = "free_x") +
  scale_fill_manual(values = ggsci::pal_d3("category20")(20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

comp_plot_card_gene_family
ggsave("output/figures/composition_plots/card/gene_family_comp.pdf",
       comp_plot_card_gene_family,
       height = 8,
       width = 9)



# By drug class ----


# Sum aro's in each gene family
drug_class_df <- rgi_df %>%
  group_by(mouse_time, treatment,drug_class) %>%
  summarise(tpm = sum(tpm), .groups = "drop")

comp_plot_card_drug_class <- drug_class_df %>%
  ggplot(aes(x = mouse_time, y = tpm, fill = drug_class)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~treatment, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = ggsci::pal_d3("category20")(20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

comp_plot_card_drug_class
ggsave("output/figures/composition_plots/card/drug_class_comp.pdf",
       comp_plot_card_drug_class,
       height = 8,
       width = 10)




# By resistance mechanism ----


# Sum aro's in each gene family
resistance_mechanism_df <- rgi_df %>%
  group_by(mouse_time, treatment,  resistance_mechanism) %>%
  summarise(tpm = sum(tpm), .groups = "drop")

comp_plot_card_resistance_mechanism <- resistance_mechanism_df %>%
  ggplot(aes(x = mouse_time, y = tpm, fill = resistance_mechanism)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~treatment , ncol = 2, scales = "free_x") +
  scale_fill_manual(values = ggsci::pal_d3("category20")(20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

comp_plot_card_resistance_mechanism
ggsave("output/figures/composition_plots/card/resistance_mechanism_comp.pdf",
       comp_plot_card_resistance_mechanism,
       height = 8,
       width = 10)
