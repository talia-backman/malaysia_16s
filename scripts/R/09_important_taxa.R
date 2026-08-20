# 09_important_taxa.R
#
# identify habitat-associated and coast-associated genera
# run corncob differential abundance models
# save figure 3

set.seed(666)

library(tidyverse)
library(phyloseq)
library(janitor)
library(patchwork)
library(corncob)

# load helper functions
source("./scripts/R/functions.R")

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/important_taxa"
fig_dir <- "./figures/figure3_important_taxa"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# load phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# prepare metadata factors
sample_data(ps)$habitat <- factor(sample_data(ps)$habitat)
sample_data(ps)$east_or_west <- factor(sample_data(ps)$east_or_west)
sample_data(ps)$site_code <- factor(sample_data(ps)$site_code)
sample_data(ps)$protected_area <- factor(sample_data(ps)$protected_area)

message("samples: ", nsamples(ps))
message("asvs: ", ntaxa(ps))

# agglomerate to genus level
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 0, ps_genus)
ps_genus@refseq <- taxa_names(ps_genus) %>% Biostrings::DNAStringSet()
taxa_names(ps_genus) <- paste0("Taxon_", seq_len(ntaxa(ps_genus)))

saveRDS(ps_genus, file.path(out_dir, "phyloseq_genus_level.RDS"))
message("genus-level taxa: ", ntaxa(ps_genus))

# prevalence filter genus-level taxa
prev_threshold <- ceiling(0.05 * nsamples(ps_genus))
tax_prev <- apply(otu_table(ps_genus), 2, function(x) sum(x > 0))
ps_genus_filt <- prune_taxa(tax_prev >= prev_threshold, ps_genus)
ps_genus_filt <- prune_taxa(taxa_sums(ps_genus_filt) > 0, ps_genus_filt)

saveRDS(ps_genus_filt, file.path(out_dir, "phyloseq_genus_level_prev5pct.RDS"))

message("prevalence threshold: ", prev_threshold, " samples")
message("genus-level taxa before filtering: ", ntaxa(ps_genus))
message("genus-level taxa after filtering: ", ntaxa(ps_genus_filt))

# prepare taxonomy table
retained_taxa <- tax_table(ps_genus_filt) %>% as("matrix") %>% as.data.frame(check.names = FALSE)
colnames(retained_taxa) <- make.unique(colnames(retained_taxa))

# test habitat-associated taxa
habitat_da <- corncob::differentialTest(formula = ~ habitat, phi.formula = ~ 1, formula_null = ~ 1, phi.formula_null = ~ 1,
                                        test = "Wald", boot = FALSE, data = ps_genus_filt, fdr_cutoff = 0.05, full_output = TRUE)

saveRDS(habitat_da, file.path(out_dir, "corncob_habitat_diff_analysis.RDS"))
names(habitat_da$full_output) <- taxa_names(ps_genus_filt)

habitat_sig_taxa <- habitat_da$significant_taxa
message("significant habitat-associated taxa: ", length(habitat_sig_taxa))

habitat_sig_table <- retained_taxa %>% rownames_to_column("tax_row") %>% mutate(taxon_id = taxa_names(ps_genus_filt)) %>%
  filter(taxon_id %in% habitat_sig_taxa) %>% dplyr::select(taxon_id, Phylum, Class, Order, Family, Genus) %>%
  arrange(Phylum, Family, Genus) %>% mutate(taxonomy = case_when(!is.na(Genus) ~ paste(Family, Genus, sep = ": "),
                                                                 !is.na(Family) ~ Family, !is.na(Order) ~ Order, TRUE ~ Phylum))

write_csv(habitat_sig_table, file.path(out_dir, "habitat_significant_taxa.csv"))

# summarize habitat-associated taxa
habitat_sig_melt <- ps_genus_filt %>% transform_sample_counts(function(x) x / sum(x)) %>%
  subset_taxa(taxa_names(ps_genus_filt) %in% habitat_sig_taxa) %>% psmelt() %>%
  left_join(habitat_sig_table %>% dplyr::select(OTU = taxon_id, taxonomy), by = "OTU")

habitat_sig_summary <- habitat_sig_melt %>% group_by(OTU, taxonomy, habitat) %>%
  summarise(mean_rel_abund = mean(Abundance, na.rm = TRUE), median_rel_abund = median(Abundance, na.rm = TRUE),
            prevalence = mean(Abundance > 0), .groups = "drop") %>%
  group_by(OTU, taxonomy) %>% mutate(max_habitat = habitat[which.max(mean_rel_abund)],
                                     max_mean_rel_abund = max(mean_rel_abund)) %>% ungroup() %>% arrange(max_habitat, desc(max_mean_rel_abund))

write_csv(habitat_sig_summary, file.path(out_dir, "habitat_significant_taxa_relative_abundance_summary.csv"))

habitat_sig_winners <- habitat_sig_summary %>% group_by(OTU, taxonomy) %>% slice_max(mean_rel_abund, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(OTU, taxonomy, max_habitat = habitat, max_mean_rel_abund = mean_rel_abund, prevalence)

write_csv(habitat_sig_winners, file.path(out_dir, "habitat_significant_taxa_highest_habitat.csv"))

# test coast-associated taxa
coast_da <- corncob::differentialTest(formula = ~ east_or_west, phi.formula = ~ 1, formula_null = ~ 1, phi.formula_null = ~ 1,
                                      test = "Wald", boot = FALSE, data = ps_genus_filt, fdr_cutoff = 0.05, full_output = TRUE)

saveRDS(coast_da, file.path(out_dir, "corncob_eastwest_diff_analysis.RDS"))
names(coast_da$full_output) <- taxa_names(ps_genus_filt)

coast_sig_taxa <- coast_da$significant_taxa
message("significant east-west-associated taxa: ", length(coast_sig_taxa))

coast_sig_table <- retained_taxa %>% rownames_to_column("tax_row") %>% mutate(taxon_id = taxa_names(ps_genus_filt)) %>%
  filter(taxon_id %in% coast_sig_taxa) %>% dplyr::select(taxon_id, Phylum, Class, Order, Family, Genus) %>%
  arrange(Phylum, Family, Genus) %>% mutate(taxonomy = case_when(!is.na(Genus) ~ paste(Family, Genus, sep = ": "),
                                                                 !is.na(Family) ~ Family, !is.na(Order) ~ Order, TRUE ~ Phylum))

write_csv(coast_sig_table, file.path(out_dir, "eastwest_significant_taxa.csv"))

# summarize coast-associated taxa
coast_sig_melt <- ps_genus_filt %>% transform_sample_counts(function(x) x / sum(x)) %>%
  subset_taxa(taxa_names(ps_genus_filt) %in% coast_sig_taxa) %>% psmelt() %>%
  left_join(coast_sig_table %>% dplyr::select(OTU = taxon_id, taxonomy), by = "OTU")

coast_sig_summary <- coast_sig_melt %>% group_by(OTU, taxonomy, east_or_west) %>%
  summarise(mean_rel_abund = mean(Abundance, na.rm = TRUE), median_rel_abund = median(Abundance, na.rm = TRUE),
            prevalence = mean(Abundance > 0), .groups = "drop") %>%
  group_by(OTU, taxonomy) %>% mutate(max_coast = east_or_west[which.max(mean_rel_abund)],
                                     max_mean_rel_abund = max(mean_rel_abund)) %>% ungroup() %>% arrange(max_coast, desc(max_mean_rel_abund))

write_csv(coast_sig_summary, file.path(out_dir, "eastwest_significant_taxa_relative_abundance_summary.csv"))

# make heatmaps
heatmap_df <- habitat_sig_summary %>% dplyr::select(taxonomy, habitat, mean_rel_abund, max_habitat) %>% distinct() %>%
  mutate(max_habitat = factor(max_habitat, levels = c("Coral", "Mangrove", "Seagrass"))) %>%
  arrange(max_habitat, desc(mean_rel_abund)) %>% mutate(taxonomy = factor(taxonomy, levels = unique(taxonomy)))

coast_heatmap_df <- coast_sig_summary %>% dplyr::select(taxonomy, east_or_west, mean_rel_abund, max_coast) %>% distinct() %>%
  mutate(max_coast = factor(max_coast, levels = c("East", "West"))) %>%
  arrange(max_coast, desc(mean_rel_abund)) %>% mutate(taxonomy = factor(taxonomy, levels = unique(taxonomy)))

shared_fill_limits <- c(0, max(heatmap_df$mean_rel_abund, coast_heatmap_df$mean_rel_abund, na.rm = TRUE))

p_habitat_heatmap <- ggplot(heatmap_df, aes(x = habitat, y = taxonomy, fill = mean_rel_abund)) +
  geom_tile(color = "white") + scale_fill_viridis_c(option = "magma", trans = "sqrt", limits = shared_fill_limits) +
  labs(x = NULL, y = NULL, fill = "mean\nrelative\nabundance", tag = "A") + theme_bw() +
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), panel.grid = element_blank())

p_coast_heatmap <- ggplot(coast_heatmap_df, aes(x = east_or_west, y = taxonomy, fill = mean_rel_abund)) +
  geom_tile(color = "white") + scale_fill_viridis_c(option = "magma", trans = "sqrt", limits = shared_fill_limits) +
  labs(x = NULL, y = NULL, fill = "mean\nrelative\nabundance", tag = "B") + theme_bw() +
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), panel.grid = element_blank())

figure3_taxa <- (p_habitat_heatmap | p_coast_heatmap) + plot_layout(guides = "collect", widths = c(1.5, 1)) &
  theme(legend.position = "right", panel.grid = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

figure3_taxa

ggsave(file.path(fig_dir, "figure3_important_taxa.pdf"), figure3_taxa, width = 15, height = 11)
ggsave(file.path(fig_dir, "figure3_important_taxa.png"), figure3_taxa, width = 15, height = 11, dpi = 600)

# combined habitat and coast model
combined_da <- corncob::differentialTest(formula = ~ habitat + east_or_west, phi.formula = ~ 1, formula_null = ~ 1,
                                         phi.formula_null = ~ 1, test = "Wald", boot = FALSE, data = ps_genus_filt,
                                         fdr_cutoff = 0.05, full_output = TRUE)

saveRDS(combined_da, file.path(out_dir, "corncob_combined_habitat_coast_analysis.RDS"))
names(combined_da$full_output) <- taxa_names(ps_genus_filt)

combined_sig_taxa <- combined_da$significant_taxa
message("significant taxa in combined habitat + coast model: ", length(combined_sig_taxa))

# compare significant taxa across models
sig_overlap_tbl <- tibble(taxon_id = unique(c(habitat_sig_taxa, coast_sig_taxa, combined_sig_taxa))) %>%
  mutate(sig_habitat = taxon_id %in% habitat_sig_taxa, sig_coast = taxon_id %in% coast_sig_taxa,
         sig_combined = taxon_id %in% combined_sig_taxa,
         category = case_when(sig_habitat & sig_coast ~ "habitat_and_coast", sig_habitat & !sig_coast ~ "habitat_only",
                              !sig_habitat & sig_coast ~ "coast_only", TRUE ~ "combined_only")) %>%
  left_join(retained_taxa %>% rownames_to_column("tax_row") %>% mutate(taxon_id = taxa_names(ps_genus_filt)) %>%
              dplyr::select(taxon_id, Phylum, Class, Order, Family, Genus), by = "taxon_id") %>%
  mutate(taxonomy = case_when(!is.na(Genus) ~ paste(Family, Genus, sep = ": "), !is.na(Family) ~ Family,
                              !is.na(Order) ~ Order, TRUE ~ Phylum)) %>% arrange(category, taxonomy)

write_csv(sig_overlap_tbl, file.path(out_dir, "significant_taxa_model_overlap.csv"))

message("script 9 done. outputs in: ", out_dir)