# 16_ancombc2_sensitivity.R
#
# evaluate habitat and coast differential abundance using ancom-bc2
# account for site-level clustering while retaining individual samples
# compare site-aware results with original corncob findings

set.seed(666)

library(tidyverse)
library(phyloseq)
library(ANCOMBC)
library(microbiome)

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
corncob_path <- "./output/important_taxa/significant_taxa_model_overlap.csv"
out_dir <- "./output/ancombc2_sensitivity"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# load cleaned phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# aggregate to genus level and retain taxa present in at least 5% of samples
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 0, ps_genus)

otu_mat <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) otu_mat <- t(otu_mat)

tax_prev <- colSums(otu_mat > 0)
ps_genus <- prune_taxa(names(tax_prev)[tax_prev >= ceiling(0.05 * nsamples(ps_genus))], ps_genus)

# prepare metadata
sample_data(ps_genus)$habitat <- factor(sample_data(ps_genus)$habitat, levels = c("Mangrove", "Coral", "Seagrass"))
sample_data(ps_genus)$east_or_west <- factor(sample_data(ps_genus)$east_or_west, levels = c("West", "East"))
sample_data(ps_genus)$site_code <- factor(sample_data(ps_genus)$site_code)

stopifnot(nsamples(ps_genus) == 74, nlevels(sample_data(ps_genus)$site_code) == 19)

# assign stable taxon IDs and save taxonomy lookup
tax_df <- as.data.frame(as(tax_table(ps_genus), "matrix"), check.names = FALSE)
colnames(tax_df) <- make.unique(colnames(tax_df))

tax_lookup <- tax_df %>% rownames_to_column("original_asv") %>%
  mutate(taxon_id = paste0("Taxon_", row_number())) %>%
  dplyr::select(taxon_id, original_asv, Phylum, Class, Order, Family, Genus)

taxa_names(ps_genus) <- tax_lookup$taxon_id
write_csv(tax_lookup, file.path(out_dir, "ancombc2_taxonomy_lookup.csv"))

# fit site-aware ancom-bc2 model
ancom_res <- ANCOMBC::ancombc2(data = ps_genus, tax_level = NULL, fix_formula = "habitat + east_or_west",
                               rand_formula = "(1 | site_code)", p_adj_method = "BH", prv_cut = 0,
                               lib_cut = 0, group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                               global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                               pseudo_sens = TRUE, n_cl = 1, verbose = TRUE)

saveRDS(ancom_res, file.path(out_dir, "ancombc2_site_random_effect_results.RDS"))

# export complete coefficient-level results
ancom_tbl <- ancom_res$res

write_csv(ancom_tbl, file.path(out_dir, "ancombc2_all_coefficients.csv"))

# summarize significant coefficient-level associations
ancom_sig <- ancom_tbl %>%
  pivot_longer(cols = starts_with("diff_"), names_to = "coefficient", values_to = "significant") %>%
  filter(significant %in% TRUE) %>%
  dplyr::select(taxon, coefficient) %>%
  left_join(tax_lookup, by = c("taxon" = "taxon_id"))

write_csv(ancom_sig, file.path(out_dir, "ancombc2_significant_coefficients.csv"))

message("significant ancom-bc2 coefficient-level associations: ", nrow(ancom_sig))
message("script 16 done. outputs in: ", out_dir)