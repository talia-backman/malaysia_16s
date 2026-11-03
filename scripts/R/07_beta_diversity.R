# 07_beta_diversity.R
#
# test microbial community differences by habitat, coast, site, and protection
# assess geographic envfit and multivariate dispersion
# test whether seagrass is intermediate using habitat centroid distances
# summarize community structure statistics for table s3 and centroid distances for table s4c
# save figure 1

set.seed(666)

library(phyloseq)
library(vegan)
library(tidyverse)
library(janitor)
library(broom)
library(patchwork)

# load helper functions
source("./scripts/R/functions.R")

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/beta_diversity"
fig_dir <- "./figures/figure1_beta_diversity"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# load rarefied object and add phylogenetic tree
# tree was built from the filtered ASV set in script 05
ps_rare <- readRDS(file.path("./output", "phyloseq_rarefied_16s.rds"))
ps_tree <- readRDS(file.path("./output", "16s_physeq_cleaned_w_tree.RDS"))

# add tree to rarefied object -- subset tree tips to rarefied ASV set
shared_taxa <- intersect(taxa_names(ps_rare), taxa_names(ps_tree))
ps <- merge_phyloseq(ps_rare, phy_tree(ps_tree))
ps <- prune_taxa(taxa_names(ps) %in% shared_taxa, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

message("samples after merging tree: ", nsamples(ps))
message("taxa after merging tree: ", ntaxa(ps))

# prepare metadata
meta_raw <- as(sample_data(ps), "data.frame")
meta_raw[] <- lapply(meta_raw, function(x) if (is.list(x)) unlist(x) else x)
meta <- meta_raw %>% janitor::clean_names() %>% rownames_to_column("sample_name") %>%
  mutate(east_or_west = factor(east_or_west), habitat = factor(habitat), site_code = factor(site_code), protected_area = factor(protected_area))

ps_sample_names <- as.character(sample_names(ps))
meta <- meta[match(ps_sample_names, meta$sample_name), ]
stopifnot(identical(meta$sample_name, ps_sample_names))

# calculate unifrac distance and nmds
dist_unif <- phyloseq::distance(ps, method = "unifrac")
ord_nmds <- ordinate(ps, method = "NMDS", distance = dist_unif)

scr <- as.data.frame(ord_nmds$points) %>% rownames_to_column("sample_name") %>% left_join(meta, by = "sample_name")

habitat_cols <- c(Coral = "#440154FF", Mangrove = "#21908CFF", Seagrass = "#FDE725FF")
coast_cols <- c(East = "#39568CFF", West = "#95D840FF")

# habitat centroid analysis using betadisper principal coordinate space
# tests whether seagrass is compositionally intermediate between mangrove and coral communities
bd_habitat_centroids <- vegan::betadisper(dist_unif, group = meta$habitat, type = "centroid")

centroid_coords <- as.data.frame(bd_habitat_centroids$centroids) %>% rownames_to_column("habitat") %>%
  filter(habitat %in% c("Coral", "Mangrove", "Seagrass"))

centroid_distances <- centroid_coords %>% column_to_rownames("habitat") %>%
  dist() %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("habitat_1") %>%
  pivot_longer(-habitat_1, names_to = "habitat_2", values_to = "unifrac_centroid_distance") %>%
  filter(habitat_1 < habitat_2)

write_csv(centroid_coords, file.path(out_dir, "habitat_unifrac_centroids.csv"))
write_csv(centroid_distances, file.path(out_dir, "habitat_unifrac_centroid_distances.csv"))

seagrass_intermediate_summary <- centroid_distances %>%
  mutate(comparison = paste(habitat_1, habitat_2, sep = "_vs_")) %>%
  dplyr::select(comparison, unifrac_centroid_distance) %>% pivot_wider(names_from = comparison,
              values_from = unifrac_centroid_distance) %>%
  mutate(seagrass_closer_to_coral_than_mangrove_is =
           Coral_vs_Seagrass < Coral_vs_Mangrove,
         seagrass_closer_to_mangrove_than_coral_is =
           Mangrove_vs_Seagrass < Coral_vs_Mangrove,
         seagrass_intermediate =
           seagrass_closer_to_coral_than_mangrove_is &
           seagrass_closer_to_mangrove_than_coral_is)

write_csv(seagrass_intermediate_summary, file.path(out_dir, "seagrass_intermediate_unifrac_summary.csv"))

site_scores <- as.data.frame(bd_habitat_centroids$vectors) %>% rownames_to_column("sample_name") %>%
  left_join(meta, by = "sample_name") %>% filter(habitat %in% c("Coral", "Mangrove", "Seagrass"))

p_habitat_unifrac_centroids <- ggplot(site_scores, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(fill = habitat), shape = 21, color = "black", size = 3, alpha = 0.6, stroke = 0.3) +
  geom_point(data = centroid_coords, aes(x = PCoA1, y = PCoA2, fill = habitat),
             shape = 23, color = "black", size = 5, stroke = 1) +
  geom_text(data = centroid_coords, aes(x = PCoA1, y = PCoA2, label = habitat),
            nudge_y = 0.04, size = 4) +
  scale_fill_manual(values = habitat_cols) + theme_bw() +
  labs(x = "PCoA1", y = "PCoA2", fill = "habitat")

p_habitat_unifrac_centroids

ggsave(file.path(fig_dir, "habitat_unifrac_centroids.pdf"), p_habitat_unifrac_centroids, width = 6, height = 5)
ggsave(file.path(fig_dir, "habitat_unifrac_centroids.png"), p_habitat_unifrac_centroids, width = 6, height = 5, dpi = 300)

# plot nmds panels
p_habitat <- ggplot(scr, aes(x = MDS1, y = MDS2)) + geom_point(aes(fill = habitat), shape = 21, color = "black", size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(color = habitat), level = 0.95, linewidth = 0.7) + scale_fill_manual(values = habitat_cols) +
  scale_color_manual(values = habitat_cols) + theme_bw() + labs(title = "A", x = "NMDS1", y = "NMDS2", fill = "habitat", color = "habitat")

p_coast <- ggplot(scr, aes(x = MDS1, y = MDS2)) + geom_point(aes(fill = east_or_west), shape = 21, color = "black", size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(color = east_or_west), level = 0.95, linewidth = 0.7) + scale_fill_manual(values = coast_cols) +
  scale_color_manual(values = coast_cols) + theme_bw() + labs(title = "B", x = "NMDS1", y = "NMDS2", fill = "coast", color = "coast")

p_habitat_by_coast <- ggplot(scr, aes(x = MDS1, y = MDS2)) + geom_point(aes(fill = habitat), shape = 21, color = "black", size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(color = habitat, group = habitat), level = 0.95, linewidth = 0.7) + facet_wrap(~ east_or_west) +
  scale_fill_manual(values = habitat_cols) + scale_color_manual(values = habitat_cols) + theme_bw() +
  labs(title = "C", x = "NMDS1", y = "NMDS2", fill = "habitat", color = "habitat")

figure1_beta <- (p_habitat | p_coast) / p_habitat_by_coast + plot_layout(heights = c(1, 1.1))
figure1_beta

ggsave(file.path(fig_dir, "figure1_beta_diversity.pdf"), figure1_beta, width = 10, height = 8)
ggsave(file.path(fig_dir, "figure1_beta_diversity.png"), figure1_beta, width = 10, height = 8, dpi = 300)

# fit latitude and longitude vectors to nmds ordination
geo_vars <- meta %>%
  transmute(longitude = as.numeric(as.character(longitude)),
            latitude = as.numeric(as.character(latitude)))
stopifnot(!anyNA(geo_vars$longitude), !anyNA(geo_vars$latitude))
geo_envfit <- envfit(ord_nmds, geo_vars, permutations = 999)
geo_envfit_tbl <- tibble(variable = rownames(geo_envfit$vectors$arrows),
                         r2 = geo_envfit$vectors$r,
                         p_value = geo_envfit$vectors$pvals)
write_csv(geo_envfit_tbl, file.path(out_dir, "geographic_envfit_summary.csv"))
geo_envfit_tbl

# sample-level PERMANOVA (reported alongside site-level as sensitivity check)
adon_eastwest <- adonis2(dist_unif ~ east_or_west, data = meta, permutations = 999)
adon_habitat <- adonis2(dist_unif ~ habitat, data = meta, permutations = 999)
adon_site <- adonis2(dist_unif ~ site_code, data = meta, permutations = 999)
adon_protected <- adonis2(dist_unif ~ protected_area, data = meta, permutations = 999)
adon_additive <- adonis2(dist_unif ~ east_or_west + habitat, data = meta, permutations = 999)
adon_interaction <- adonis2(dist_unif ~ east_or_west * habitat, data = meta, permutations = 999)
adon_full_additive <- adonis2(dist_unif ~ east_or_west + habitat + protected_area, data = meta, permutations = 999)
adon_spatial_prot <- adonis2(dist_unif ~ east_or_west * habitat + protected_area, data = meta, permutations = 999)
adon_prot_siteblock <- adonis2(dist_unif ~ protected_area, data = meta, permutations = 999, strata = meta$site_code)
capture.output(adon_eastwest, file = file.path(out_dir, "adonis2_unifrac_east_or_west.txt"))
capture.output(adon_habitat, file = file.path(out_dir, "adonis2_unifrac_habitat.txt"))
capture.output(adon_site, file = file.path(out_dir, "adonis2_unifrac_site_code.txt"))
capture.output(adon_protected, file = file.path(out_dir, "adonis2_unifrac_protected_area.txt"))
capture.output(adon_additive, file = file.path(out_dir, "adonis2_unifrac_east_or_west_plus_habitat.txt"))
capture.output(adon_interaction, file = file.path(out_dir, "adonis2_unifrac_east_or_west_by_habitat.txt"))
capture.output(adon_full_additive, file = file.path(out_dir, "adonis2_unifrac_east_or_west_plus_habitat_plus_protected_area.txt"))
capture.output(adon_spatial_prot, file = file.path(out_dir, "adonis2_unifrac_east_or_west_by_habitat_plus_protected_area.txt"))
capture.output(adon_prot_siteblock, file = file.path(out_dir, "adonis2_unifrac_protected_area_blocked_by_site.txt"))

# site-level PERMANOVA
# aggregate samples to site-level mean relative abundance profiles
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

site_meta <- as(sample_data(ps_rel), "data.frame") %>%
  janitor::clean_names() %>%
  rownames_to_column("sample_name") %>%
  distinct(site_code, east_or_west, habitat, protected_area)

otu_site <- as.data.frame(as(otu_table(ps_rel), "matrix"))
if (taxa_are_rows(ps_rel)) otu_site <- t(otu_site)

otu_site <- as.data.frame(otu_site) %>% rownames_to_column("sample_name") %>%
  left_join(as(sample_data(ps_rel), "data.frame") %>% rownames_to_column("sample_name") %>% 
              select(sample_name, site_code), by = "sample_name") %>%
  group_by(site_code) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("site_code")

site_meta <- site_meta[match(rownames(otu_site), site_meta$site_code), ]
stopifnot(identical(rownames(otu_site), site_meta$site_code))

dist_site <- vegdist(otu_site, method = "bray")

adon_site_coast <- adonis2(dist_site ~ east_or_west, data = site_meta, permutations = 999)
adon_site_habitat <- adonis2(dist_site ~ habitat, data = site_meta, permutations = 999)
adon_site_additive <- adonis2(dist_site ~ east_or_west + habitat, data = site_meta, permutations = 999)

message("site-level PERMANOVA results:")
print(adon_site_coast)
print(adon_site_habitat)
print(adon_site_additive)

write_csv(bind_rows(tidy_adonis(adon_site_coast, "site_level_coast"),
    tidy_adonis(adon_site_habitat, "site_level_habitat"), tidy_adonis(adon_site_additive, "site_level_additive")),
  file.path(out_dir, "permanova_site_level_summary.csv"))

perm_tbl <- bind_rows(tidy_adonis(adon_eastwest, "east_or_west"), tidy_adonis(adon_habitat, "habitat"),
                      tidy_adonis(adon_site, "site_code"), tidy_adonis(adon_protected, "protected_area"),
                      tidy_adonis(adon_additive, "east_or_west_plus_habitat"), tidy_adonis(adon_interaction, "east_or_west_by_habitat"),
                      tidy_adonis(adon_full_additive, "east_or_west_plus_habitat_plus_protected_area"),
                      tidy_adonis(adon_spatial_prot, "east_or_west_by_habitat_plus_protected_area"),
                      tidy_adonis(adon_prot_siteblock, "protected_area_blocked_by_site"))

write_csv(perm_tbl, file.path(out_dir, "permanova_unifrac_summary.csv"))

# betadisper tests
bd_eastwest <- run_betadisp(dist_unif, meta$east_or_west, "east_or_west")
bd_habitat <- run_betadisp(dist_unif, meta$habitat, "habitat")
bd_site <- run_betadisp(dist_unif, meta$site_code, "site_code")
bd_protected <- run_betadisp(dist_unif, meta$protected_area, "protected_area")

capture.output(bd_eastwest$pt, file = file.path(out_dir, "betadisper_unifrac_east_or_west.txt"))
capture.output(bd_habitat$pt, file = file.path(out_dir, "betadisper_unifrac_habitat.txt"))
capture.output(bd_site$pt, file = file.path(out_dir, "betadisper_unifrac_site_code.txt"))
capture.output(bd_protected$pt, file = file.path(out_dir, "betadisper_unifrac_protected_area.txt"))

betadisp_summary <- bind_rows(bd_eastwest$summary, bd_habitat$summary, bd_site$summary, bd_protected$summary)
betadisp_distances <- bind_rows(bd_eastwest$distances, bd_habitat$distances, bd_site$distances, bd_protected$distances)

write_csv(betadisp_summary, file.path(out_dir, "betadisper_unifrac_summary.csv"))
write_csv(betadisp_distances, file.path(out_dir, "betadisper_unifrac_distances.csv"))

# build table S3 community structure statistics
table_s3_permanova <- perm_tbl %>% filter(model %in% c("site_code", "habitat", "east_or_west", "protected_area"), term != "Residual", term != "Total") %>%
  transmute(analysis = "PERMANOVA", predictor = recode(model, site_code = "site identity", habitat = "habitat identity",
                                                       east_or_west = "east-west coast identity", protected_area = "protected-area status"),
            df, r2, statistic, p_value = p.value)
table_s3_envfit <- geo_envfit_tbl %>% transmute(analysis = "geographic envfit", predictor = variable, df = NA_real_,
                                                r2, statistic = NA_real_, p_value)
table_s3_betadisper <- betadisp_summary %>% filter(model %in% c("east_or_west", "habitat", "protected_area")) %>%
  transmute(analysis = "betadisper", predictor = recode(model, east_or_west = "east-west coast identity",
                                                        habitat = "habitat identity", protected_area = "protected-area status"),
            df = NA_real_, r2 = NA_real_, statistic = F, p_value = p)
table_s3 <- bind_rows(table_s3_permanova, table_s3_envfit, table_s3_betadisper)
write_csv(table_s3, file.path(out_dir, "table_s3_community_structure_statistics.csv"))
table_s3

message("script 7 done. outputs in: ", out_dir)
