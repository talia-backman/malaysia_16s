# 13_nestedness_turnover.R
#
# sort presence-absence beta diversity into turnover and nestedness components
# test whether community differences reflect taxon replacement or subset structure

set.seed(666)

library(tidyverse)
library(phyloseq)
library(vegan)
library(betapart)
library(janitor)

# load helper functions
source("./scripts/R/functions.R")

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/nestedness_turnover"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# load phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# prepare metadata
meta <- as(sample_data(ps), "data.frame") %>% janitor::clean_names() %>% rownames_to_column("sample_name") %>%
  mutate(east_or_west = factor(east_or_west), habitat = factor(habitat), site_code = factor(site_code), protected_area = factor(protected_area))

meta <- meta %>% slice(match(sample_names(ps), sample_name))
stopifnot(identical(meta$sample_name, sample_names(ps)))

# convert otu table to presence-absence matrix
otu_pa <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_pa <- t(otu_pa) }
otu_pa[otu_pa > 0] <- 1
otu_pa <- as.data.frame(otu_pa)

# sort sorensen beta diversity
beta_parts <- betapart::beta.pair(otu_pa, index.family = "sorensen")
dist_turnover <- beta_parts$beta.sim
dist_nestedness <- beta_parts$beta.sne
dist_total_sor <- beta_parts$beta.sor

# permanova models for turnover, nestedness, and total sorensen beta diversity
adon_turnover_protected <- adonis2(dist_turnover ~ protected_area, data = meta, permutations = 999)
adon_nestedness_protected <- adonis2(dist_nestedness ~ protected_area, data = meta, permutations = 999)
adon_total_sor_protected <- adonis2(dist_total_sor ~ protected_area, data = meta, permutations = 999)

adon_turnover_protected_siteblock <- adonis2(dist_turnover ~ protected_area, data = meta, permutations = 999, strata = meta$site_code)
adon_nestedness_protected_siteblock <- adonis2(dist_nestedness ~ protected_area, data = meta, permutations = 999, strata = meta$site_code)

adon_turnover_eastwest <- adonis2(dist_turnover ~ east_or_west, data = meta, permutations = 999)
adon_nestedness_eastwest <- adonis2(dist_nestedness ~ east_or_west, data = meta, permutations = 999)

adon_turnover_habitat <- adonis2(dist_turnover ~ habitat, data = meta, permutations = 999)
adon_nestedness_habitat <- adonis2(dist_nestedness ~ habitat, data = meta, permutations = 999)

# save full outputs
capture.output(adon_turnover_protected, file = file.path(out_dir, "adonis2_turnover_protected_area.txt"))
capture.output(adon_nestedness_protected, file = file.path(out_dir, "adonis2_nestedness_protected_area.txt"))
capture.output(adon_total_sor_protected, file = file.path(out_dir, "adonis2_total_sorensen_protected_area.txt"))
capture.output(adon_turnover_protected_siteblock, file = file.path(out_dir, "adonis2_turnover_protected_area_blocked_by_site.txt"))
capture.output(adon_nestedness_protected_siteblock, file = file.path(out_dir, "adonis2_nestedness_protected_area_blocked_by_site.txt"))
capture.output(adon_turnover_eastwest, file = file.path(out_dir, "adonis2_turnover_east_or_west.txt"))
capture.output(adon_nestedness_eastwest, file = file.path(out_dir, "adonis2_nestedness_east_or_west.txt"))
capture.output(adon_turnover_habitat, file = file.path(out_dir, "adonis2_turnover_habitat.txt"))
capture.output(adon_nestedness_habitat, file = file.path(out_dir, "adonis2_nestedness_habitat.txt"))

# make summary table
nestedness_tbl <- bind_rows(
  tidy_adonis(adon_turnover_protected, "turnover_protected_area"),
  tidy_adonis(adon_nestedness_protected, "nestedness_protected_area"),
  tidy_adonis(adon_total_sor_protected, "total_sorensen_protected_area"),
  tidy_adonis(adon_turnover_protected_siteblock, "turnover_protected_area_blocked_by_site"),
  tidy_adonis(adon_nestedness_protected_siteblock, "nestedness_protected_area_blocked_by_site"),
  tidy_adonis(adon_turnover_eastwest, "turnover_east_or_west"),
  tidy_adonis(adon_nestedness_eastwest, "nestedness_east_or_west"),
  tidy_adonis(adon_turnover_habitat, "turnover_habitat"),
  tidy_adonis(adon_nestedness_habitat, "nestedness_habitat")
)

write_csv(nestedness_tbl, file.path(out_dir, "nestedness_turnover_permanova_summary.csv"))

# pairwise turnover and nestedness summaries among focal habitats
pairwise_beta_df <- tibble(
  sample_1 = rownames(as.matrix(dist_total_sor))[row(as.matrix(dist_total_sor))[lower.tri(as.matrix(dist_total_sor))]],
  sample_2 = colnames(as.matrix(dist_total_sor))[col(as.matrix(dist_total_sor))[lower.tri(as.matrix(dist_total_sor))]],
  beta_sor = as.matrix(dist_total_sor)[lower.tri(as.matrix(dist_total_sor))],
  beta_sim_turnover = as.matrix(dist_turnover)[lower.tri(as.matrix(dist_turnover))],
  beta_sne_nestedness = as.matrix(dist_nestedness)[lower.tri(as.matrix(dist_nestedness))]) %>%
  left_join(meta %>% select(sample_name, habitat, east_or_west), by = c("sample_1" = "sample_name")) %>%
  rename(habitat_1 = habitat, coast_1 = east_or_west) %>%
  left_join(meta %>% select(sample_name, habitat, east_or_west), by = c("sample_2" = "sample_name")) %>%
  rename(habitat_2 = habitat, coast_2 = east_or_west) %>%
  filter(habitat_1 %in% c("Coral", "Mangrove", "Seagrass"), habitat_2 %in% c("Coral", "Mangrove", "Seagrass"),
         habitat_1 != habitat_2) %>% mutate(habitat_pair = map2_chr(habitat_1, habitat_2,
                            ~ paste(sort(c(as.character(.x), as.character(.y))), collapse = "_vs_")),
    same_coast = coast_1 == coast_2)

pairwise_beta_summary <- pairwise_beta_df %>% group_by(habitat_pair) %>%
  summarize(n_pairwise_comparisons = n(), mean_total_sorensen = mean(beta_sor, na.rm = TRUE),
    median_total_sorensen = median(beta_sor, na.rm = TRUE), mean_turnover = mean(beta_sim_turnover, na.rm = TRUE),
    median_turnover = median(beta_sim_turnover, na.rm = TRUE), mean_nestedness = mean(beta_sne_nestedness, na.rm = TRUE),
    median_nestedness = median(beta_sne_nestedness, na.rm = TRUE), .groups = "drop")

write_csv(pairwise_beta_df, file.path(out_dir, "pairwise_sample_beta_components_focal_habitats.csv"))
write_csv(pairwise_beta_summary, file.path(out_dir, "pairwise_beta_components_focal_habitat_summary.csv"))

seagrass_bridge_beta_summary <- pairwise_beta_summary %>%
  select(habitat_pair, mean_total_sorensen, mean_turnover, mean_nestedness) %>%
  pivot_wider(names_from = habitat_pair, values_from = c(mean_total_sorensen, mean_turnover, mean_nestedness)) %>%
  mutate(seagrass_lower_total_dissimilarity_than_coral_mangrove =
      mean_total_sorensen_Coral_vs_Seagrass < mean_total_sorensen_Coral_vs_Mangrove &
      mean_total_sorensen_Mangrove_vs_Seagrass < mean_total_sorensen_Coral_vs_Mangrove,
    seagrass_lower_turnover_than_coral_mangrove =
      mean_turnover_Coral_vs_Seagrass < mean_turnover_Coral_vs_Mangrove &
      mean_turnover_Mangrove_vs_Seagrass < mean_turnover_Coral_vs_Mangrove)

write_csv(seagrass_bridge_beta_summary, file.path(out_dir, "seagrass_bridge_beta_component_summary.csv"))

message("script 13 done. outputs in: ", out_dir)