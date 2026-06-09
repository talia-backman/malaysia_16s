# 14_core_microbiome_comparison.R
#
# compare taxa from the published core mangrove microbiome to taxa identified in this study
# focus on overlap with habitat-associated, coast-associated, and combined significant taxa

set.seed(666)

library(tidyverse)
library(janitor)

# set paths
core_path <- "./data/core_microbiome/mangrove_core_dataframe.csv"
important_taxa_dir <- "./output/important_taxa"
out_dir <- "./output/core_microbiome_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# load mangrove core microbiome dataframe and important taxa outputs
mangrove_core_raw <- read_csv(core_path, show_col_types = FALSE) %>% janitor::clean_names()
habitat_sig <- read_csv(file.path(important_taxa_dir, "habitat_significant_taxa.csv"), show_col_types = FALSE) %>% 
  janitor::clean_names()
coast_sig <- read_csv(file.path(important_taxa_dir, "eastwest_significant_taxa.csv"), show_col_types = FALSE) %>% 
  janitor::clean_names()
model_overlap <- read_csv(file.path(important_taxa_dir, "significant_taxa_model_overlap.csv"), show_col_types = FALSE) %>% 
  janitor::clean_names()

# helper to standardize genus names
clean_genus <- function(x) {
  x %>% as.character() %>% str_trim() %>% na_if("") %>% str_replace("^Candidatus\\s+", "Candidatus ") %>% 
    str_squish()
}

# extract unique core mangrove genera
mangrove_core_genera <- mangrove_core_raw %>% mutate(genus_clean = clean_genus(genus)) %>%
  filter(!is.na(genus_clean), genus_clean != "NA", genus_clean != "uncultured", genus_clean != "unclassified") %>%
  distinct(genus_clean) %>% arrange(genus_clean)

write_csv(mangrove_core_genera, file.path(out_dir, "mangrove_core_genera.csv"))

# standardize genus names in malaysia significant taxa tables
habitat_sig_clean <- habitat_sig %>% mutate(genus_clean = clean_genus(genus), source_table = "habitat_significant_taxa")
coast_sig_clean <- coast_sig %>% mutate(genus_clean = clean_genus(genus), source_table = "eastwest_significant_taxa")
model_overlap_clean <- model_overlap %>% mutate(genus_clean = clean_genus(genus), source_table = "significant_taxa_model_overlap")

# compare core mangrove genera to significant taxa
habitat_overlap <- habitat_sig_clean %>% filter(!is.na(genus_clean), genus_clean %in% mangrove_core_genera$genus_clean) %>%
  mutate(overlap_set = "habitat-associated taxa")

coast_overlap <- coast_sig_clean %>% filter(!is.na(genus_clean), genus_clean %in% mangrove_core_genera$genus_clean) %>%
  mutate(overlap_set = "coast-associated taxa")

combined_overlap <- model_overlap_clean %>% filter(!is.na(genus_clean), genus_clean %in% mangrove_core_genera$genus_clean) %>%
  mutate(overlap_set = "any significant taxa")

all_overlap <- bind_rows(habitat_overlap %>% dplyr::select(overlap_set, taxon_id, phylum, class, order, 
                                                           family, genus, genus_clean, taxonomy, everything()),
  coast_overlap %>% dplyr::select(overlap_set, taxon_id, phylum, class, order, family, genus, genus_clean, 
                                  taxonomy, everything()),
  combined_overlap %>% dplyr::select(overlap_set, taxon_id, phylum, class, order, family, genus, genus_clean, 
                                     taxonomy, category, sig_habitat, sig_coast, sig_combined, everything()))

write_csv(all_overlap, file.path(out_dir, "mangrove_core_overlap_with_significant_taxa.csv"))

# summarize overlap counts
overlap_summary <- tibble(comparison = c("mangrove core genera", "habitat-associated significant taxa", 
                                         "coast-associated significant taxa", "all significant taxa"),
  n_total = c(n_distinct(mangrove_core_genera$genus_clean), n_distinct(habitat_sig_clean$genus_clean, na.rm = TRUE),
              n_distinct(coast_sig_clean$genus_clean, na.rm = TRUE), n_distinct(model_overlap_clean$genus_clean, na.rm = TRUE)),
  n_overlap_with_mangrove_core = c(NA_integer_, n_distinct(habitat_overlap$genus_clean), n_distinct(coast_overlap$genus_clean),
                                   n_distinct(combined_overlap$genus_clean)))

write_csv(overlap_summary, file.path(out_dir, "mangrove_core_overlap_summary.csv"))

# make publication overlap table
overlap_table <- combined_overlap %>% dplyr::select(taxon_id, phylum, class, order, family, genus, taxonomy, 
                                                    category, sig_habitat, sig_coast, sig_combined) %>%
  arrange(category, genus)

write_csv(overlap_table, file.path(out_dir, "mangrove_core_overlap_publication_table.csv"))