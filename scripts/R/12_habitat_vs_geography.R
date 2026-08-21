# 12_habitat_vs_geography.R
#
# directly compare habitat identity and coast identity as predictors of community similarity
# repeat comparisons after excluding within-site pairs and at the site level
# provide sample- and site-level comparisons for tables s4a and s4b
# save figure s4

set.seed(666)

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/habitat_vs_geography"
fig_dir <- "./figures/supplemental/figure_s4_habitat_vs_geography"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# load phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# extract metadata
meta <- sample_data(ps) %>% data.frame() %>% rownames_to_column("physeq_sample_id") %>% as_tibble()

meta_small <- meta %>% dplyr::select(physeq_sample_id, sample_id, habitat, east_or_west, site_code, latitude, longitude) %>%
  mutate(habitat = as.factor(habitat), east_or_west = as.factor(east_or_west), site_code = as.factor(site_code),
         latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# calculate unweighted UniFrac distances
dist_unifrac <- phyloseq::distance(ps, method = "unifrac")
dist_mat <- as.matrix(dist_unifrac)

# convert distance matrix to unique pairwise table
pairwise_dist <- dist_mat %>% as.data.frame() %>% rownames_to_column("sample_1") %>%
  pivot_longer(cols = -sample_1, names_to = "sample_2", values_to = "distance") %>% filter(sample_1 != sample_2) %>%
  mutate(pair_id = map2_chr(sample_1, sample_2, ~ paste(sort(c(.x, .y)), collapse = "__"))) %>% distinct(pair_id, .keep_all = TRUE)

# add metadata for both samples in each pair
pairwise_df <- pairwise_dist %>% left_join(meta_small, by = c("sample_1" = "physeq_sample_id")) %>%
  rename(habitat_1 = habitat, coast_1 = east_or_west, site_1 = site_code, latitude_1 = latitude,
         longitude_1 = longitude, original_sample_id_1 = sample_id) %>%
  left_join(meta_small, by = c("sample_2" = "physeq_sample_id")) %>%
  rename(habitat_2 = habitat, coast_2 = east_or_west, site_2 = site_code, latitude_2 = latitude,
         longitude_2 = longitude, original_sample_id_2 = sample_id) %>%
  filter(!is.na(distance))

# classify pair types
pairwise_df <- pairwise_df %>% mutate(same_habitat = habitat_1 == habitat_2, same_coast = coast_1 == coast_2, same_site = site_1 == site_2,
                                      pair_category = case_when(same_habitat & same_coast ~ "Same habitat / same coast",
                                                                !same_habitat & same_coast ~ "Different habitat / same coast",
                                                                same_habitat & !same_coast ~ "Same habitat / different coast",
                                                                !same_habitat & !same_coast ~ "Different habitat / different coast", TRUE ~ NA_character_),
                                      pair_category = factor(pair_category, levels = c("Same habitat / same coast", "Different habitat / same coast",
                                                                                       "Same habitat / different coast", "Different habitat / different coast")))

pairwise_summary <- pairwise_df %>% group_by(pair_category) %>%
  summarise(n_pairs = n(), mean_distance = mean(distance, na.rm = TRUE), median_distance = median(distance, na.rm = TRUE),
            sd_distance = sd(distance, na.rm = TRUE), .groups = "drop")

write_csv(pairwise_df, file.path(out_dir, "pairwise_habitat_geography_distances.csv"))
write_csv(pairwise_summary, file.path(out_dir, "pairwise_habitat_geography_summary.csv"))

# direct contrast between habitat-sharing and coast-sharing pairs
direct_contrast_df <- pairwise_df %>% filter(pair_category %in% c("Different habitat / same coast", "Same habitat / different coast")) %>%
  mutate(contrast_group = case_when(pair_category == "Different habitat / same coast" ~ "Different habitat\nsame coast",
                                    pair_category == "Same habitat / different coast" ~ "Same habitat\ndifferent coast"),
         contrast_group = factor(contrast_group, levels = c("Different habitat\nsame coast", "Same habitat\ndifferent coast")))

direct_contrast_summary <- direct_contrast_df %>% group_by(contrast_group) %>%
  summarise(n_pairs = n(), mean_distance = mean(distance, na.rm = TRUE), median_distance = median(distance, na.rm = TRUE),
            sd_distance = sd(distance, na.rm = TRUE), .groups = "drop")

wilcox_direct_contrast <- wilcox.test(distance ~ contrast_group, data = direct_contrast_df)

write_csv(direct_contrast_df, file.path(out_dir, "direct_habitat_vs_geography_contrast_distances.csv"))
write_csv(direct_contrast_summary, file.path(out_dir, "direct_habitat_vs_geography_contrast_summary.csv"))
capture.output(wilcox_direct_contrast, file = file.path(out_dir, "direct_habitat_vs_geography_wilcox_test.txt"))

# repeat after removing within-site comparisons
pairwise_df_no_site <- pairwise_df %>% filter(!same_site)

pairwise_summary_no_site <- pairwise_df_no_site %>% group_by(pair_category) %>%
  summarise(n_pairs = n(), mean_distance = mean(distance, na.rm = TRUE), median_distance = median(distance, na.rm = TRUE),
            sd_distance = sd(distance, na.rm = TRUE), .groups = "drop")

direct_contrast_no_site <- pairwise_df_no_site %>% filter(pair_category %in% c("Different habitat / same coast", "Same habitat / different coast")) %>%
  mutate(contrast_group = case_when(pair_category == "Different habitat / same coast" ~ "Different habitat\nsame coast",
                                    pair_category == "Same habitat / different coast" ~ "Same habitat\ndifferent coast"),
         contrast_group = factor(contrast_group, levels = c("Different habitat\nsame coast", "Same habitat\ndifferent coast")))

direct_contrast_summary_no_site <- direct_contrast_no_site %>% group_by(contrast_group) %>%
  summarise(n_pairs = n(), mean_distance = mean(distance, na.rm = TRUE), median_distance = median(distance, na.rm = TRUE),
            sd_distance = sd(distance, na.rm = TRUE), .groups = "drop")

wilcox_direct_contrast_no_site <- wilcox.test(distance ~ contrast_group, data = direct_contrast_no_site)

write_csv(pairwise_summary_no_site, file.path(out_dir, "pairwise_habitat_geography_summary_no_site.csv"))
write_csv(direct_contrast_summary_no_site, file.path(out_dir, "direct_habitat_vs_geography_summary_no_site.csv"))
capture.output(wilcox_direct_contrast_no_site, file = file.path(out_dir, "direct_habitat_vs_geography_wilcox_no_site.txt"))

# site-level robustness check
ps_site <- merge_samples(ps, "site_code")
sample_data(ps_site)$site_code <- sample_names(ps_site)

site_meta <- meta_small %>% group_by(site_code) %>% summarise(habitat = first(habitat), east_or_west = first(east_or_west), .groups = "drop")
site_dist <- phyloseq::distance(ps_site, method = "unifrac")

site_df <- as.matrix(site_dist) %>% as.data.frame() %>% rownames_to_column("site_1") %>%
  pivot_longer(-site_1, names_to = "site_2", values_to = "distance") %>% filter(site_1 != site_2) %>%
  mutate(pair_id = map2_chr(site_1, site_2, ~ paste(sort(c(.x, .y)), collapse = "__"))) %>% distinct(pair_id, .keep_all = TRUE) %>%
  left_join(site_meta, by = c("site_1" = "site_code")) %>% rename(habitat_1 = habitat, coast_1 = east_or_west) %>%
  left_join(site_meta, by = c("site_2" = "site_code")) %>% rename(habitat_2 = habitat, coast_2 = east_or_west) %>%
  mutate(same_habitat = habitat_1 == habitat_2, same_coast = coast_1 == coast_2,
         pair_category = case_when(same_habitat & same_coast ~ "Same habitat / same coast",
                                   !same_habitat & same_coast ~ "Different habitat / same coast",
                                   same_habitat & !same_coast ~ "Same habitat / different coast",
                                   !same_habitat & !same_coast ~ "Different habitat / different coast", TRUE ~ NA_character_),
         pair_category = factor(pair_category, levels = c("Same habitat / same coast", "Different habitat / same coast",
                                                          "Same habitat / different coast", "Different habitat / different coast")))

site_summary <- site_df %>% group_by(pair_category) %>%
  summarise(n_pairs = n(), mean_distance = mean(distance), median_distance = median(distance), sd_distance = sd(distance), .groups = "drop")

site_contrast <- site_df %>% filter(pair_category %in% c("Different habitat / same coast", "Same habitat / different coast")) %>%
  mutate(contrast_group = ifelse(pair_category == "Different habitat / same coast", "Different habitat\nsame coast", "Same habitat\ndifferent coast"))

site_contrast_summary <- site_contrast %>% group_by(contrast_group) %>%
  summarise(n_pairs = n(), mean_distance = mean(distance), median_distance = median(distance), sd_distance = sd(distance), .groups = "drop")

site_wilcox <- wilcox.test(distance ~ contrast_group, data = site_contrast)

write_csv(site_df, file.path(out_dir, "site_level_habitat_geography_distances.csv"))
write_csv(site_summary, file.path(out_dir, "site_level_habitat_geography_summary.csv"))
write_csv(site_contrast_summary, file.path(out_dir, "site_level_habitat_vs_geography_contrast_summary.csv"))
capture.output(site_wilcox, file = file.path(out_dir, "site_level_habitat_vs_geography_wilcox.txt"))

# plot figure s4
p_pairwise_categories_no_site <- ggplot(pairwise_df_no_site, aes(x = pair_category, y = distance)) +
  geom_boxplot(outlier.alpha = 0.15) + labs(x = NULL, y = "Unweighted UniFrac distance", tag = "A") +
  theme_bw() + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
                     axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 13), plot.tag = element_text(face = "bold", size = 16),
                     panel.grid.minor = element_blank())

p_site_categories <- ggplot(site_df, aes(x = pair_category, y = distance)) +
  geom_boxplot(outlier.alpha = 0.15) + labs(x = NULL, y = "Unweighted UniFrac distance", tag = "B") +
  theme_bw() + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
                     axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 13), plot.tag = element_text(face = "bold", size = 16),
                     panel.grid.minor = element_blank())

# save figure s4
figure_s4_habitat_vs_geography <- p_pairwise_categories_no_site +
  labs(tag = NULL)

figure_s4_habitat_vs_geography

ggsave(file.path(fig_dir, "figure_s4_habitat_vs_geography.pdf"), figure_s4_habitat_vs_geography, width = 7, height = 5)
ggsave(file.path(fig_dir, "figure_s4_habitat_vs_geography.png"), figure_s4_habitat_vs_geography, width = 7, height = 5, dpi = 600)

message("script 12 done. outputs in: ", out_dir)
