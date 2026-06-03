# 13_nestedness_turnover.R
#
# partition presence-absence beta diversity into turnover and nestedness components
# test whether community differences reflect taxon replacement or subset structure
# save supplemental results

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

# partition sorensen beta diversity
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

message("done. outputs in: ", out_dir)