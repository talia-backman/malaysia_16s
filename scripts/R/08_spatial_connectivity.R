# 08_spatial_connectivity.R
#
# calculate ocean least-cost distances
# test distance-decay with mantel tests
# save figure 2

set.seed(666)

library(tidyverse)
library(janitor)
library(sf)
library(raster)
library(gdistance)
library(rnaturalearth)
library(rnaturalearthdata)
library(phyloseq)
library(vegan)
library(patchwork)
library(ggspatial)
library(viridis)

# load helper functions
source("./scripts/R/functions.R")

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/spatial_connectivity"
fig_dir <- "./figures/figure2_spatial_connectivity"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# load phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# prepare metadata and coordinates
meta_ps <- as(sample_data(ps), "data.frame") %>% janitor::clean_names() %>% rownames_to_column("sample_name")

pts_df <- meta_ps %>% mutate(longitude = str_trim(as.character(longitude)), latitude = str_trim(as.character(latitude)),
                             longitude = str_replace_all(longitude, "[^0-9\\-\\.]+", ""), latitude = str_replace_all(latitude, "[^0-9\\-\\.]+", ""),
                             lon = suppressWarnings(as.numeric(longitude)), lat = suppressWarnings(as.numeric(latitude))) %>%
  filter(!is.na(lon), !is.na(lat)) %>% transmute(sample_name, east_or_west = str_to_lower(as.character(east_or_west)),
                                                 habitat = str_to_lower(as.character(habitat)), lon, lat) %>% distinct(sample_name, .keep_all = TRUE)

# build ocean least-cost distance matrix
pts_sf <- st_as_sf(pts_df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

bb <- st_bbox(pts_sf)
bb_buf <- bb
bb_buf[c("xmin", "ymin")] <- bb[c("xmin", "ymin")] - 3
bb_buf[c("xmax", "ymax")] <- bb[c("xmax", "ymax")] + 3
bb_poly <- st_as_sfc(bb_buf)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land <- st_intersection(world, st_transform(bb_poly, st_crs(world))) %>% st_make_valid()

crs_m <- 3857
pts_m <- st_transform(pts_sf, crs_m)
land_m <- st_transform(land, crs_m)
bb_m <- st_bbox(st_transform(bb_poly, crs_m))

res_m <- 1000
shrink_m <- 4000

land_m2 <- suppressWarnings(st_buffer(land_m, dist = -shrink_m)) %>% st_make_valid()

r <- raster(xmn = bb_m["xmin"], xmx = bb_m["xmax"], ymn = bb_m["ymin"], ymx = bb_m["ymax"], res = res_m, crs = st_crs(crs_m)$wkt)
land_sp <- as(land_m2, "Spatial")
land_r <- rasterize(land_sp, r, field = 1, background = 0)

cost_r <- land_r
cost_r[land_r[] == 1] <- NA
cost_r[land_r[] == 0] <- 1

# snap nearshore points to nearest water cell
pts_sp <- as(pts_m, "Spatial")
cell_ids <- raster::cellFromXY(cost_r, sp::coordinates(pts_sp))
on_na <- is.na(cost_r[cell_ids])
message("points on land/na before snapping: ", sum(on_na), " / ", length(on_na))

if (any(on_na)) {
  water_cells <- which(!is.na(cost_r[]))
  water_xy <- raster::xyFromCell(cost_r, water_cells)
  coords_mat <- sp::coordinates(pts_sp)
  rownames(coords_mat) <- pts_df$sample_name
  
  for (i in which(on_na)) {
    p <- coords_mat[i, , drop = FALSE]
    d2 <- (water_xy[, 1] - p[1, 1])^2 + (water_xy[, 2] - p[1, 2])^2
    coords_mat[i, ] <- water_xy[which.min(d2), ]
  }
  
  pts_sp <- sp::SpatialPoints(coords_mat, proj4string = sp::CRS(sp::proj4string(pts_sp)))
}

cell_ids2 <- raster::cellFromXY(cost_r, sp::coordinates(pts_sp))
on_na2 <- is.na(cost_r[cell_ids2])
message("points on land/na after snapping: ", sum(on_na2), " / ", length(on_na2))

tr <- transition(cost_r, transitionFunction = mean, directions = 8)
tr <- geoCorrection(tr, type = "c", multpl = FALSE)

ocean_dist_m <- costDistance(tr, pts_sp)
ocean_km <- as.matrix(ocean_dist_m / 1000)
rownames(ocean_km) <- pts_df$sample_name
colnames(ocean_km) <- pts_df$sample_name

message("inf entries in ocean_km: ", sum(is.infinite(ocean_km)), " / ", length(ocean_km))

saveRDS(ocean_km, file.path(out_dir, "ocean_leastcost_distance_km.rds"))
write.csv(ocean_km, file.path(out_dir, "ocean_leastcost_distance_km.csv"), row.names = TRUE)

# build community distance matrices
dist_unif <- as.matrix(phyloseq::distance(ps, method = "unifrac"))
dist_wunif <- as.matrix(phyloseq::distance(ps, method = "wunifrac"))

rownames(dist_unif) <- sample_names(ps)
colnames(dist_unif) <- sample_names(ps)
rownames(dist_wunif) <- sample_names(ps)
colnames(dist_wunif) <- sample_names(ps)

keep <- Reduce(intersect, list(rownames(dist_unif), rownames(ocean_km), meta_ps$sample_name))
dist_unif <- dist_unif[keep, keep, drop = FALSE]
dist_wunif <- dist_wunif[keep, keep, drop = FALSE]
ocean_km <- ocean_km[keep, keep, drop = FALSE]

meta_keep <- meta_ps %>% mutate(east_or_west = str_to_lower(as.character(east_or_west)), habitat = str_to_lower(as.character(habitat))) %>%
  filter(sample_name %in% keep) %>% slice(match(keep, sample_name))

# run mantel tests
mantel_results <- list()

group_list <- list(all = meta_keep$sample_name,
                   east_only = meta_keep$sample_name[meta_keep$east_or_west == "east"],
                   west_only = meta_keep$sample_name[meta_keep$east_or_west == "west"])

habitat_levels <- meta_keep %>% filter(!is.na(habitat), habitat != "") %>% distinct(habitat) %>% pull(habitat) %>% sort()

for (h in habitat_levels) { group_list[[paste0("habitat_", h)]] <- meta_keep$sample_name[meta_keep$habitat == h] }

for (nm in names(group_list)) {
  mantel_results[[paste0("unifrac_", nm)]] <- run_mantel(ids = group_list[[nm]], label = nm, dist_mat = dist_unif, geo_mat = ocean_km, dist_name = "unifrac")
  mantel_results[[paste0("wunifrac_", nm)]] <- run_mantel(ids = group_list[[nm]], label = nm, dist_mat = dist_wunif, geo_mat = ocean_km, dist_name = "wunifrac")
}

mantel_tbl <- bind_rows(mantel_results)
write_csv(mantel_tbl, file.path(out_dir, "mantel_unifrac_wunifrac_vs_ocean_summary.csv"))

# build pairwise distance-decay data
pair_unif <- dist_to_df(dist_unif, ocean_km, "unifrac", meta_keep)
pair_wunif <- dist_to_df(dist_wunif, ocean_km, "wunifrac", meta_keep)
pair_df <- bind_rows(pair_unif, pair_wunif)

write_csv(pair_df, file.path(out_dir, "distance_decay_pairwise_data.csv"))

# plot all-sample distance-decay
p_all <- ggplot(pair_df, aes(x = ocean_distance_km, y = community_distance)) +
  geom_point(alpha = 0.35, size = 1) + geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Ocean-connected distance (km)", y = "community distance", tag = "B") + theme_bw()

ggsave(file.path(fig_dir, "distance_decay_all_samples.pdf"), p_all, width = 9, height = 4.5)
ggsave(file.path(fig_dir, "distance_decay_all_samples.png"), p_all, width = 9, height = 4.5, dpi = 300)

# plot habitat-specific distance-decay
pair_df_same_habitat <- pair_df %>% filter(habitat_1 == habitat_2) %>% mutate(habitat = habitat_1)

p_habitat <- ggplot(pair_df_same_habitat, aes(x = ocean_distance_km, y = community_distance)) +
  geom_point(alpha = 0.35, size = 1) + geom_smooth(method = "loess", se = TRUE) +
  facet_grid(metric ~ habitat, scales = "free_y") +
  labs(x = "Ocean-connected distance (km)", y = "community distance", tag = "C") + theme_bw()

ggsave(file.path(fig_dir, "distance_decay_by_habitat.pdf"), p_habitat, width = 12, height = 6)
ggsave(file.path(fig_dir, "distance_decay_by_habitat.png"), p_habitat, width = 12, height = 6, dpi = 300)

# inspect plotted site coordinates before adjusting map display positions
coord_check <- pts_df %>% dplyr::select(east_or_west, habitat, lon, lat) %>%
  dplyr::distinct() %>% dplyr::arrange(east_or_west, lon, lat)
coord_check
coord_counts <- pts_df %>% dplyr::count(east_or_west, habitat, lon, lat, name = "n_samples") %>%
  dplyr::arrange(east_or_west, lon, lat)
coord_counts

# build map panel
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
malaysia_map <- world_sf %>% filter(admin %in% c("Malaysia", "Thailand", "Indonesia", "Singapore", "Brunei", "Philippines")) %>% st_make_valid()
habitat_cols <- c(coral = "#440154FF", mangrove = "#21908CFF", seagrass = "#FDE725FF")
# habitat_cols_map <- scales::alpha(habitat_cols, 0.45)
pts_df <- pts_df %>% dplyr::left_join(
  meta_keep %>% dplyr::select(sample_name, site_code, site_replicate),
  by = "sample_name"
)
# use one actual sampling coordinate per site for map display
pts_sites <- pts_df %>% dplyr::group_by(site_code, east_or_west, habitat) %>%
  dplyr::slice_min(site_replicate, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()
# map
p_map <- ggplot() + geom_sf(data = malaysia_map, linewidth = 0.2, fill = "grey90", color = "grey50") +
  geom_jitter(data = pts_sites, aes(x = lon, y = lat, fill = habitat), shape = 21, size = 3, alpha = 0.85,
              height = 0.1, width = 0.1, color = "black", stroke = 0.3) +
  coord_sf(xlim = range(pts_sites$lon) + c(-1.5, 1.5), ylim = range(pts_sites$lat) + c(-1.5, 1.5), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  scale_fill_manual(values = habitat_cols, guide = guide_legend(override.aes = list(
    fill = unname(habitat_cols), alpha = 1))) +
  labs(x = "longitude", y = "latitude", fill = "habitat", tag = "A") + theme_bw()
p_map
# combine figure 2
figure2_spatial <- p_map | (p_all / p_habitat)
figure2_spatial <- figure2_spatial &
  theme(text = element_text(size = 14),
        plot.tag = element_text(face = "bold", size = 16))

figure2_spatial

ggsave(file.path(fig_dir, "figure2_spatial_connectivity.pdf"), figure2_spatial, width = 11, height = 5)
ggsave(file.path(fig_dir, "figure2_spatial_connectivity.png"), figure2_spatial, width = 11, height = 5, dpi = 300)

message("script 8 done. outputs in: ", out_dir)