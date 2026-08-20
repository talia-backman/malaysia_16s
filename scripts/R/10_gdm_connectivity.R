# 10_gdm_connectivity.R
#
# fit generalized dissimilarity models using presence-absence community data
# compare euclidean geographic distance and ocean-connected distance with mrm
# save figure 4

set.seed(666)

library(tidyverse)
library(phyloseq)
library(gdm)
library(janitor)
library(vegan)
library(ecodist)
library(patchwork)

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
ocean_dist_path <- "./output/spatial_connectivity/ocean_leastcost_distance_km.rds"
out_dir <- "./output/gdm_connectivity"
fig_dir <- "./figures/figure4_gdm_connectivity"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# load phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# prepare metadata and presence-absence community matrix
meta <- data.frame(sample_data(ps)) %>% janitor::clean_names()
meta$sample_id <- rownames(meta)

otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_mat <- t(otu_mat) }
otu_mat <- otu_mat[meta$sample_id, , drop = FALSE]

pa_mat <- otu_mat
pa_mat[pa_mat > 0] <- 1

sample_sums_pa <- rowSums(pa_mat)
empty_samples <- names(sample_sums_pa[sample_sums_pa == 0])
message("empty samples removed from gdm: ", paste(empty_samples, collapse = ", "))

pa_mat_gdm <- pa_mat[sample_sums_pa > 0, , drop = FALSE]
meta_gdm <- meta %>% filter(sample_id %in% rownames(pa_mat_gdm))
pa_mat_gdm <- pa_mat_gdm[meta_gdm$sample_id, , drop = FALSE]

# build gdm site-pair input
gdm_comm <- pa_mat_gdm %>% as.data.frame() %>% rownames_to_column("sample_id")

gdm_site_data <- meta_gdm %>% dplyr::select(sample_id, longitude, latitude, habitat, east_or_west, protected_area) %>%
  mutate(habitat_coral = as.numeric(habitat == "Coral"), habitat_mangrove = as.numeric(habitat == "Mangrove"),
         habitat_seagrass = as.numeric(habitat == "Seagrass"), coast_east = as.numeric(east_or_west == "East"),
         protected_yes = as.numeric(protected_area == "Y")) %>%
  dplyr::select(sample_id, longitude, latitude, habitat_coral, habitat_mangrove, habitat_seagrass, coast_east, protected_yes)

gdm_input <- formatsitepair(bioData = gdm_comm, bioFormat = 1, XColumn = "longitude", YColumn = "latitude",
                            predData = gdm_site_data, siteColumn = "sample_id")

gdm_input <- gdm_input %>% mutate(distance = as.numeric(distance), weights = as.numeric(weights),
                                  s1.xCoord = as.numeric(s1.xCoord), s1.yCoord = as.numeric(s1.yCoord),
                                  s2.xCoord = as.numeric(s2.xCoord), s2.yCoord = as.numeric(s2.yCoord))

# fit gdm models
gdm_geo <- gdm(gdm_input[, 1:6], geo = TRUE)
gdm_full <- gdm(gdm_input, geo = TRUE)

gdm_model_summary <- tibble(model = c("Euclidean geography only", "Euclidean geography + habitat/coast/protected area"),
                            explained_deviance = c(gdm_geo$explained, gdm_full$explained))

write_csv(gdm_model_summary, file.path(out_dir, "gdm_presence_absence_model_summary.csv"))

capture.output(summary(gdm_geo), file = file.path(out_dir, "gdm_geography_only_summary.txt"))
capture.output(summary(gdm_full), file = file.path(out_dir, "gdm_full_model_summary.txt"))

# predictor importance
gdm_importance <- gdm.varImp(gdm_input, geo = TRUE, nPerm = 100, parallel = FALSE)

gdm_importance_tbl <- tibble(predictor = rownames(gdm_importance$`Predictor Importance`),
                             importance = gdm_importance$`Predictor Importance`[, 1],
                             p_value = gdm_importance$`Predictor p-values`[, 1]) %>% arrange(desc(importance))

write_csv(gdm_importance_tbl, file.path(out_dir, "gdm_presence_absence_importance.csv"))

# full gdm spline plot for supplement
png(file.path(fig_dir, "gdm_presence_absence_splines.png"), width = 2200, height = 1800, res = 300)
plot(gdm_full, plot.layout = c(3, 2), rug = TRUE, ylab = "compositional turnover", xlab = "predictor gradient")
dev.off()

pdf(file.path(fig_dir, "gdm_presence_absence_splines.pdf"), width = 7.3, height = 6)
plot(gdm_full, plot.layout = c(3, 2), rug = TRUE, ylab = "compositional turnover", xlab = "predictor gradient")
dev.off()

# compare euclidean geography and ocean-connected distance
ocean_dist <- readRDS(ocean_dist_path)
ocean_dist <- as.matrix(ocean_dist)
ocean_dist_gdm <- ocean_dist[meta_gdm$sample_id, meta_gdm$sample_id]

comm_dist_pa <- vegdist(pa_mat_gdm, method = "bray", binary = TRUE)
geo_dist <- dist(meta_gdm %>% dplyr::select(longitude, latitude))
ocean_dist_as_dist <- as.dist(ocean_dist_gdm)

mrm_geo <- MRM(comm_dist_pa ~ geo_dist, nperm = 999)
mrm_ocean <- MRM(comm_dist_pa ~ ocean_dist_as_dist, nperm = 999)
mrm_geo_ocean <- MRM(comm_dist_pa ~ geo_dist + ocean_dist_as_dist, nperm = 999)

mrm_compare_tbl <- tibble(model = c("Euclidean geographic distance", "Ocean-connected distance", "Both predictors"),
                          r2 = c(mrm_geo$r.squared["R2"], mrm_ocean$r.squared["R2"], mrm_geo_ocean$r.squared["R2"]),
                          p_value = c(mrm_geo$r.squared["pval"], mrm_ocean$r.squared["pval"], mrm_geo_ocean$r.squared["pval"])) %>%
  mutate(model = factor(model, levels = model))

mrm_coef_tbl <- bind_rows(as.data.frame(mrm_geo$coef) %>% rownames_to_column("term") %>% mutate(model = "Euclidean geographic distance"),
                          as.data.frame(mrm_ocean$coef) %>% rownames_to_column("term") %>% mutate(model = "Ocean-connected distance"),
                          as.data.frame(mrm_geo_ocean$coef) %>% rownames_to_column("term") %>% mutate(model = "Both predictors")) %>% relocate(model, term)

write_csv(mrm_compare_tbl, file.path(out_dir, "mrm_ocean_distance_comparison.csv"))
write_csv(mrm_coef_tbl, file.path(out_dir, "mrm_ocean_distance_coefficients.csv"))

capture.output(mrm_geo, file = file.path(out_dir, "mrm_euclidean_geography.txt"))
capture.output(mrm_ocean, file = file.path(out_dir, "mrm_ocean_least_cost.txt"))
capture.output(mrm_geo_ocean, file = file.path(out_dir, "mrm_both_distance_predictors.txt"))

# prepare ocean-connected distance relationship
ocean_plot_tbl <- tibble(ocean_distance_km = as.vector(ocean_dist_as_dist),
                         community_dissimilarity = as.vector(comm_dist_pa))

write_csv(ocean_plot_tbl, file.path(out_dir, "ocean_connected_distance_decay.csv"))

# build figure 4 panels
gdm_importance_plot_tbl <- gdm_importance_tbl %>% mutate(predictor = recode(predictor,
                                                                            Geographic = "Euclidean geographic distance", coast_east = "East-west coast",
                                                                            habitat_mangrove = "Mangrove habitat", habitat_coral = "Coral habitat", protected_yes = "Protected area"),
                                                         sig_label = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ "ns"),
                                                         predictor = factor(predictor, levels = predictor[order(importance)]))

p_gdm_importance <- ggplot(gdm_importance_plot_tbl, aes(x = predictor, y = importance)) +
  geom_col(width = 0.75) + geom_text(aes(label = sig_label), hjust = -0.2, size = 4) + coord_flip() +
  labs(x = NULL, y = "Permutation importance", title = "A") + xlim(gdm_importance_plot_tbl$predictor) +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 16))

p_ocean_distance <- ggplot(ocean_plot_tbl, aes(x = ocean_distance_km, y = community_dissimilarity)) +
  geom_point(alpha = 0.15, size = 1) + geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(x = "Ocean-connected distance (km)", y = "Community dissimilarity", title = "B") +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 16))

p_mrm_compare <- ggplot(mrm_compare_tbl, aes(x = model, y = r2)) + geom_col(width = 0.7) +
  geom_text(aes(label = paste0("R² = ", round(r2, 3))), vjust = -0.4, size = 4) +
  labs(x = NULL, y = expression(R^2), title = "C") + ylim(0, 0.115) + theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold", size = 16))

figure4_gdm <- wrap_plots(p_gdm_importance, p_ocean_distance, p_mrm_compare, ncol = 3, widths = c(1.2, 1, 1))

figure4_gdm

ggsave(file.path(fig_dir, "figure4_gdm_connectivity.pdf"), figure4_gdm, width = 15, height = 4.5)
ggsave(file.path(fig_dir, "figure4_gdm_connectivity.png"), figure4_gdm, width = 12, height = 4.5,
       dpi = 600, device = ragg::agg_png)

message("script 10 done. outputs in: ", out_dir)