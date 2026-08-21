# 11_environmental_analysis.R
#
# extract bio-oracle environmental variables and summarize environmental structure
# fit environmental vectors onto microbial nmds ordination
# summarize environmental envfit statistics for table s7
# save figure 5 and figure s3

set.seed(666)

library(tidyverse)
library(phyloseq)
library(janitor)
library(sdmpredictors)
library(raster)
library(ggplot2)
library(vegan)
library(ggrepel)

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
env_dir <- "./data/biooracle"
out_dir <- "./output/environmental_analysis"
fig_dir <- "./figures/figure5_environment"
fig_supp_dir <- "./figures/supplemental/figure_s3_environment_pca"

dir.create(env_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_supp_dir, showWarnings = FALSE, recursive = TRUE)

options(sdmpredictors_datadir = env_dir)

# load phyloseq object and metadata
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

meta <- sample_data(ps) %>% data.frame() %>% janitor::clean_names() %>%
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

coords <- meta %>% dplyr::select(sample_id, latitude, longitude, habitat, east_or_west, site_code, protected_area, control) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>% distinct()

write_csv(coords, file.path(out_dir, "biooracle_sample_coordinates.csv"))

# selected bio-oracle layers
selected_layer_codes <- c("BO22_tempmean_ss", "BO22_temprange_ss", "BO22_salinitymean_ss", "BO22_salinityrange_ss",
                          "BO22_chlomean_ss", "BO22_chlorange_ss", "BO22_nitratemean_ss", "BO22_phosphatemean_ss",
                          "BO22_dissoxmean_ss", "BO22_ph", "BO22_parmean")

env_layers <- load_layers(selected_layer_codes)

# extract environmental values
env_values <- raster::extract(env_layers, coords %>% dplyr::select(longitude, latitude)) %>% as_tibble()
env_data <- bind_cols(coords, env_values) %>% mutate(env_missing = if_any(starts_with("BO22"), is.na))

write_csv(env_data, file.path(out_dir, "biooracle_environmental_data_raw.csv"))

# recover missing coastal values from nearest nearby marine raster cell
extract_nearest_ocean <- function(r, lon, lat, search_radius = 0.2) {
  pt <- matrix(c(lon, lat), ncol = 2)
  val <- raster::extract(r, pt)
  if (!is.na(val)) { return(val) }
  offsets <- expand.grid(x = seq(-search_radius, search_radius, by = 0.01), y = seq(-search_radius, search_radius, by = 0.01))
  nearby_pts <- cbind(lon + offsets$x, lat + offsets$y)
  nearby_vals <- raster::extract(r, nearby_pts)
  first_valid <- nearby_vals[!is.na(nearby_vals)][1]
  return(first_valid)
}

env_data_fixed <- env_data
env_cols <- names(env_data_fixed)[str_detect(names(env_data_fixed), "^BO22")]

for (col in env_cols) {
  r <- env_layers[[col]]
  missing_idx <- which(is.na(env_data_fixed[[col]]))
  message("fixing ", col, " (", length(missing_idx), " missing)")
  for (i in missing_idx) {
    env_data_fixed[[col]][i] <- extract_nearest_ocean(r = r, lon = env_data_fixed$longitude[i], lat = env_data_fixed$latitude[i])
  }
}

write_csv(env_data_fixed, file.path(out_dir, "biooracle_environmental_data_fixed.csv"))

# summarize environmental variables
env_site_summary <- env_data_fixed %>% group_by(site_code, habitat, east_or_west, protected_area) %>%
  summarise(across(all_of(env_cols), mean), n_samples = n(), .groups = "drop")

env_coast_summary <- env_data_fixed %>% group_by(east_or_west) %>%
  summarise(across(all_of(env_cols), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"), n_samples = n(), .groups = "drop")

env_habitat_summary <- env_data_fixed %>% group_by(habitat) %>%
  summarise(across(all_of(env_cols), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"), n_samples = n(), .groups = "drop")

write_csv(env_site_summary, file.path(out_dir, "biooracle_site_summary.csv"))
write_csv(env_coast_summary, file.path(out_dir, "biooracle_coast_summary.csv"))
write_csv(env_habitat_summary, file.path(out_dir, "biooracle_habitat_summary.csv"))

# environmental correlations and pca
env_matrix <- env_data_fixed %>% dplyr::select(all_of(env_cols))
env_scaled <- scale(env_matrix)
env_cor <- cor(env_scaled)

write_csv(as.data.frame(env_cor) %>% rownames_to_column("variable"), file.path(out_dir, "biooracle_environment_correlations.csv"))

env_pca <- prcomp(env_scaled, center = TRUE, scale. = TRUE)
env_pca_scores <- bind_cols(env_data_fixed, as_tibble(env_pca$x))

capture.output(summary(env_pca), file = file.path(out_dir, "biooracle_environment_pca_summary.txt"))

habitat_cols <- c(Coral = "#440154FF", Mangrove = "#21908CFF", Seagrass = "#FDE725FF")
coast_cols <- c(East = "#39568CFF", West = "#95D840FF")

p_env_pca <- ggplot(env_pca_scores, aes(PC1, PC2)) +
  geom_point(aes(fill = east_or_west, shape = habitat), color = "black", size = 3, alpha = 0.9, stroke = 0.3,
             show.legend = c(fill = FALSE, shape = TRUE)) +
  scale_fill_manual(values = coast_cols) + scale_shape_manual(values = c(Coral = 21, Mangrove = 24,Seagrass = 22)) +
  theme_bw() + labs(x = "PC1", y = "PC2", fill = "coast", shape = "habitat")

ggsave(file.path(fig_supp_dir, "figure_s3_environment_pca.pdf"), p_env_pca, width = 7, height = 5)
ggsave(file.path(fig_supp_dir, "figure_s3_environment_pca.png"), p_env_pca, width = 7, height = 5, dpi = 300)

# fit environmental vectors to microbial nmds
envfit_meta <- env_data_fixed %>% mutate(sample_id = paste0("S", sample_id)) %>% slice(match(sample_names(ps), sample_id))
stopifnot(all(sample_names(ps) == envfit_meta$sample_id))

envfit_matrix <- envfit_meta %>% dplyr::select(all_of(env_cols))

dist_unif <- phyloseq::distance(ps, method = "unifrac")
ord_nmds <- ordinate(ps, method = "NMDS", distance = dist_unif)

envfit_res <- envfit(ord_nmds, envfit_matrix, permutations = 999)

envfit_stats <- data.frame(variable = names(envfit_res$vectors$r), r2 = envfit_res$vectors$r,
                           p_value = envfit_res$vectors$pvals) %>% arrange(p_value)

write_csv(envfit_stats, file.path(out_dir, "envfit_statistics.csv"))

envfit_table <- envfit_stats %>% mutate(variable = recode(variable,
                                                          BO22_tempmean_ss = "Sea surface temperature (mean)", BO22_temprange_ss = "Sea surface temperature (range)",
                                                          BO22_salinitymean_ss = "Sea surface salinity (mean)", BO22_salinityrange_ss = "Sea surface salinity (range)",
                                                          BO22_chlomean_ss = "Chlorophyll concentration (mean)", BO22_chlorange_ss = "Chlorophyll concentration (range)",
                                                          BO22_nitratemean_ss = "Nitrate concentration (mean)", BO22_phosphatemean_ss = "Phosphate concentration (mean)",
                                                          BO22_dissoxmean_ss = "Dissolved oxygen (mean)", BO22_ph = "pH", BO22_parmean = "Photosynthetically available radiation"),
                                        r2 = round(r2, 3), p_value = round(p_value, 3)) %>%
  rename(environmental_variable = variable, effect_size_r2 = r2)

write_csv(envfit_table, file.path(out_dir, "supplemental_envfit_table.csv"))

# build figure 5
envfit_vectors <- scores(envfit_res, display = "vectors") %>% as.data.frame() %>% rownames_to_column("variable")

envfit_vectors_plot <- envfit_vectors %>% mutate(variable_clean = recode(variable,
                                                                         BO22_tempmean_ss = "temperature mean", BO22_temprange_ss = "temperature range",
                                                                         BO22_salinitymean_ss = "salinity mean", BO22_salinityrange_ss = "salinity range",
                                                                         BO22_chlomean_ss = "chlorophyll mean", BO22_chlorange_ss = "chlorophyll range",
                                                                         BO22_nitratemean_ss = "nitrate mean", BO22_phosphatemean_ss = "phosphate mean",
                                                                         BO22_dissoxmean_ss = "dissolved oxygen mean", BO22_ph = "pH", BO22_parmean = "PAR mean")) %>%
  left_join(envfit_stats, by = "variable") %>% filter(p_value <= 0.05) %>%
  filter(variable_clean %in% c("chlorophyll mean", "chlorophyll range", "temperature mean", "temperature range",
                               "salinity mean", "dissolved oxygen mean", "nitrate mean"))

arrow_mult <- 0.9

envfit_vectors_plot <- envfit_vectors_plot %>% mutate(NMDS1_end = NMDS1 * arrow_mult, NMDS2_end = NMDS2 * arrow_mult)

scr <- as.data.frame(ord_nmds$points) %>% rownames_to_column("sample_id") %>% left_join(envfit_meta, by = "sample_id")

p_nmds_envfit <- ggplot(scr, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = east_or_west, group = east_or_west), linewidth = 0.8, level = 0.95) +
  geom_point(aes(fill = east_or_west, shape = habitat), color = "black", size = 3, alpha = 0.9, stroke = 0.3) +
  geom_segment(data = envfit_vectors_plot, aes(x = 0, y = 0, xend = NMDS1_end, yend = NMDS2_end),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.9) +
  geom_label_repel(data = envfit_vectors_plot, aes(x = NMDS1_end, y = NMDS2_end, label = variable_clean),
                   inherit.aes = FALSE, size = 3.5, fill = "white", label.size = 0.2) +
  scale_fill_manual(values = coast_cols, guide = "none") + scale_color_manual(values = coast_cols) +
  scale_shape_manual(values = c(Coral = 21, Mangrove = 24, Seagrass = 22)) + theme_bw() +
  labs(x = "NMDS1", y = "NMDS2", color = "coast", shape = "habitat")
p_nmds_envfit

ggsave(file.path(fig_dir, "figure5_environment_envfit.pdf"), p_nmds_envfit, width = 7, height = 5)
ggsave(file.path(fig_dir, "figure5_environment_envfit.png"), p_nmds_envfit, width = 7, height = 5, dpi = 300)

message("script 11 done. outputs in: ", out_dir)
