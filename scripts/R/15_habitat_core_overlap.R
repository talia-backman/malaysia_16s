# 15_habitat_core_overlap.R
#
# define within-study habitat core asvs using a 50% prevalence threshold
# quantify core asv overlap among coral, mangrove, and seagrass habitats
# test whether the seagrass core contains coral and mangrove core asvs
# provide habitat core overlap results for table s4e

set.seed(666)

library(tidyverse)
library(phyloseq)
library(janitor)

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/habitat_core_overlap"
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
otu_pa <- as.data.frame(otu_pa) %>% rownames_to_column("sample_name")

# prepare taxonomy table
tax <- as(tax_table(ps), "matrix") %>% as.data.frame(check.names = FALSE)
names(tax) <- make.unique(names(tax))
names(tax) <- janitor::make_clean_names(names(tax))
tax <- tax %>% tibble::rownames_to_column("asv")

# define habitat cores by prevalence
core_prevalence_cutoff <- 0.5
habitat_core_asvs <- otu_pa %>%
  left_join(meta %>% select(sample_name, habitat), by = "sample_name") %>%
  filter(habitat %in% c("Coral", "Mangrove", "Seagrass")) %>%
  pivot_longer(cols = -c(sample_name, habitat), names_to = "asv", values_to = "present") %>%
  group_by(habitat, asv) %>%
  summarize(prevalence = mean(present > 0), n_samples_present = sum(present > 0), n_samples_habitat = n(), .groups = "drop") %>%
  mutate(is_core = prevalence >= core_prevalence_cutoff)

write_csv(habitat_core_asvs, file.path(out_dir, "within_study_habitat_asv_prevalence.csv"))

# summarize habitat core sizes
habitat_core_summary <- habitat_core_asvs %>% group_by(habitat) %>%
  summarize(n_samples = first(n_samples_habitat), n_core_asvs = sum(is_core), .groups = "drop")
write_csv(habitat_core_summary, file.path(out_dir, "within_study_habitat_core_summary.csv"))

# extract core ASV sets
coral_core <- habitat_core_asvs %>% filter(habitat == "Coral", is_core) %>% pull(asv)
mangrove_core <- habitat_core_asvs %>% filter(habitat == "Mangrove", is_core) %>% pull(asv)
seagrass_core <- habitat_core_asvs %>% filter(habitat == "Seagrass", is_core) %>% pull(asv)

# compare core overlap among coral, mangrove, and seagrass
core_overlap_summary <- tibble(
  comparison = c("coral_core_in_seagrass_core", "mangrove_core_in_seagrass_core",
                 "coral_core_in_mangrove_core", "mangrove_core_in_coral_core"),
  focal_core_n = c(length(coral_core), length(mangrove_core), length(coral_core), length(mangrove_core)),
  overlap_n = c(length(intersect(coral_core, seagrass_core)), length(intersect(mangrove_core, seagrass_core)),
                length(intersect(coral_core, mangrove_core)), length(intersect(mangrove_core, coral_core)))) %>%
  mutate(overlap_fraction = overlap_n / focal_core_n)
write_csv(core_overlap_summary, file.path(out_dir, "within_study_habitat_core_overlap_summary.csv"))

# identify ASVs by core membership pattern
core_overlap_asvs <- tibble(asv = unique(c(coral_core, mangrove_core, seagrass_core))) %>%
  mutate(coral_core = asv %in% coral_core, mangrove_core = asv %in% mangrove_core, seagrass_core = asv %in% seagrass_core,
         core_pattern = case_when(
           coral_core & mangrove_core & seagrass_core ~ "core in coral, mangrove, and seagrass",
           coral_core & seagrass_core & !mangrove_core ~ "core in coral and seagrass",
           mangrove_core & seagrass_core & !coral_core ~ "core in mangrove and seagrass",
           coral_core & mangrove_core & !seagrass_core ~ "core in coral and mangrove",
           coral_core & !mangrove_core & !seagrass_core ~ "core in coral only",
           mangrove_core & !coral_core & !seagrass_core ~ "core in mangrove only",
           seagrass_core & !coral_core & !mangrove_core ~ "core in seagrass only",
           TRUE ~ "not classified")) %>%
  left_join(tax, by = "asv") %>% arrange(core_pattern, phylum, class, order, family, genus)
write_csv(core_overlap_asvs, file.path(out_dir, "within_study_habitat_core_overlap_asvs.csv"))

# summarize core membership patterns
core_overlap_pattern_summary <- core_overlap_asvs %>% count(core_pattern, name = "n_asvs") %>%
  arrange(desc(n_asvs))
write_csv(core_overlap_pattern_summary, file.path(out_dir, "within_study_habitat_core_overlap_pattern_summary.csv"))

message("script 15 done. outputs in: ", out_dir)