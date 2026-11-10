# 09_important_taxa.R
#
# identify habitat-associated and coast-associated genera using corncob
# summarize significant taxa, relative abundance, prevalence, and model overlap
# provide habitat- and coast-associated taxa data for table s8
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
otu_prev <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) otu_prev <- t(otu_prev)
tax_prev <- colSums(otu_prev > 0)
ps_genus_filt <- prune_taxa(names(tax_prev)[tax_prev >= prev_threshold], ps_genus)
ps_genus_filt <- prune_taxa(taxa_sums(ps_genus_filt) > 0, ps_genus_filt)

stopifnot(ntaxa(ps_genus_filt) > 0)

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

# ============================================================
# site-level corncob sensitivity analysis
# ============================================================

# aggregate genus-level counts within each site
ps_site_da <- merge_samples(ps_genus_filt, "site_code")
ps_site_da <- prune_taxa(taxa_sums(ps_site_da) > 0, ps_site_da)

# restore site metadata from the original sample-level object
site_metadata <- as(sample_data(ps_genus_filt), "data.frame") %>%
  rownames_to_column("sample_name") %>%
  group_by(site_code) %>%
  summarise(habitat = first(as.character(habitat)), east_or_west = first(as.character(east_or_west)),
            n_samples = n(), .groups = "drop")

stopifnot(nrow(site_metadata) == 19, all(site_metadata$site_code %in% sample_names(ps_site_da)))

site_metadata <- site_metadata[match(sample_names(ps_site_da), site_metadata$site_code), ]
stopifnot(identical(as.character(site_metadata$site_code), sample_names(ps_site_da)))

sample_data(ps_site_da)$habitat <- factor(site_metadata$habitat)
sample_data(ps_site_da)$east_or_west <- factor(site_metadata$east_or_west)

message("site-level samples: ", nsamples(ps_site_da), " | taxa: ", ntaxa(ps_site_da))
print(table(sample_data(ps_site_da)$east_or_west, sample_data(ps_site_da)$habitat))

# save site-level phyloseq object and metadata
saveRDS(ps_site_da, file.path(out_dir, "phyloseq_genus_site_level.RDS"))

write_csv(site_metadata, file.path(out_dir, "site_level_da_metadata.csv"))

# check sequencing depth and replication across sites
sample_depth_df <- tibble(sample_name = sample_names(ps_genus_filt), site_code = as.character(sample_data(ps_genus_filt)$site_code),
                          habitat = as.character(sample_data(ps_genus_filt)$habitat),
                          east_or_west = as.character(sample_data(ps_genus_filt)$east_or_west),
                          sequencing_depth = as.numeric(sample_sums(ps_genus_filt)))

site_depth_summary <- sample_depth_df %>% group_by(site_code, habitat, east_or_west) %>%
  summarise(n_samples = n(), total_reads = sum(sequencing_depth), mean_reads = mean(sequencing_depth),
            min_reads = min(sequencing_depth), max_reads = max(sequencing_depth), .groups = "drop") %>%
  arrange(desc(total_reads))

write_csv(site_depth_summary, file.path(out_dir, "site_level_da_depth_summary.csv"))

message("site-level sequencing depth and replication:")
print(site_depth_summary, n = Inf)

message("range of retained samples per site: ", min(site_depth_summary$n_samples), " - ", max(site_depth_summary$n_samples))
message("range of pooled reads per site: ", min(site_depth_summary$total_reads), " - ", max(site_depth_summary$total_reads))

# run corncob models with sites as independent observations
set.seed(666)
site_habitat_da <- corncob::differentialTest(formula = ~ habitat, phi.formula = ~ 1, formula_null = ~ 1, phi.formula_null = ~ 1,
                                             test = "Wald", boot = FALSE, data = ps_site_da, fdr_cutoff = 0.05, full_output = TRUE)

set.seed(666)
site_coast_da <- corncob::differentialTest(formula = ~ east_or_west, phi.formula = ~ 1, formula_null = ~ 1, phi.formula_null = ~ 1,
                                           test = "Wald", boot = FALSE, data = ps_site_da, fdr_cutoff = 0.05, full_output = TRUE)

set.seed(666)
site_combined_da <- corncob::differentialTest(formula = ~ habitat + east_or_west, phi.formula = ~ 1, formula_null = ~ 1,
                                              phi.formula_null = ~ 1, test = "Wald", boot = FALSE, data = ps_site_da,
                                              fdr_cutoff = 0.05, full_output = TRUE)

# save full model objects
saveRDS(site_habitat_da, file.path(out_dir, "corncob_site_habitat_analysis.RDS"))
saveRDS(site_coast_da, file.path(out_dir, "corncob_site_coast_analysis.RDS"))
saveRDS(site_combined_da, file.path(out_dir, "corncob_site_combined_analysis.RDS"))

# compare significant taxa between sample-level and site-level analyses
site_sig_overlap <- tibble(taxon_id = taxa_names(ps_genus_filt)) %>%
  mutate(sample_habitat_sig = taxon_id %in% habitat_sig_taxa, site_habitat_sig = taxon_id %in% site_habitat_da$significant_taxa,
         sample_coast_sig = taxon_id %in% coast_sig_taxa, site_coast_sig = taxon_id %in% site_coast_da$significant_taxa,
         sample_combined_sig = taxon_id %in% combined_sig_taxa, site_combined_sig = taxon_id %in% site_combined_da$significant_taxa) %>%
  left_join(retained_taxa %>% rownames_to_column("tax_row") %>% mutate(taxon_id = taxa_names(ps_genus_filt)) %>%
              dplyr::select(taxon_id, Phylum, Class, Order, Family, Genus), by = "taxon_id") %>%
  mutate(taxonomy = case_when(!is.na(Genus) ~ paste(Family, Genus, sep = ": "), !is.na(Family) ~ Family,
                              !is.na(Order) ~ Order, TRUE ~ Phylum))

write_csv(site_sig_overlap, file.path(out_dir, "corncob_sample_vs_site_significance.csv"))

# summarize counts of significant taxa
site_da_summary <- tibble(model = c("habitat", "coast", "combined"),
                          sample_level_significant = c(length(habitat_sig_taxa), length(coast_sig_taxa), length(combined_sig_taxa)),
                          site_level_significant = c(length(site_habitat_da$significant_taxa), length(site_coast_da$significant_taxa),
                                                     length(site_combined_da$significant_taxa)))

write_csv(site_da_summary, file.path(out_dir, "corncob_site_sensitivity_summary.csv"))
print(site_da_summary)

message("script 9 done. outputs in: ", out_dir)