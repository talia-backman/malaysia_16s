# 07_beta_diversity.R
#
# test microbial community differences by habitat, coast, site, and protection
# assess whether permanova results reflect centroid shifts or dispersion differences
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

# load phyloseq object
ps <- readRDS(ps_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# prepare metadata
meta <- as(sample_data(ps), "data.frame") %>% janitor::clean_names() %>% rownames_to_column("sample_name") %>%
  mutate(east_or_west = factor(east_or_west), habitat = factor(habitat), site_code = factor(site_code), protected_area = factor(protected_area))

meta <- meta %>% slice(match(sample_names(ps), sample_name))
stopifnot(identical(meta$sample_name, sample_names(ps)))

# calculate unifrac distance and nmds
dist_unif <- phyloseq::distance(ps, method = "unifrac")
ord_nmds <- ordinate(ps, method = "NMDS", distance = dist_unif)

scr <- as.data.frame(ord_nmds$points) %>% rownames_to_column("sample_name") %>% left_join(meta, by = "sample_name")

habitat_cols <- c(Coral = "#440154FF", Mangrove = "#21908CFF", Seagrass = "#FDE725FF")
coast_cols <- c(East = "#39568CFF", West = "#95D840FF")
# plot nmds panels
p_habitat <- ggplot(scr, aes(x = MDS1, y = MDS2)) + geom_point(aes(fill = habitat), shape = 21, color = "black", size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(color = habitat), level = 0.75, linewidth = 0.7) + scale_fill_manual(values = habitat_cols) +
  scale_color_manual(values = habitat_cols) + theme_bw() + labs(title = "A", x = "NMDS1", y = "NMDS2", fill = "habitat", color = "habitat")

p_coast <- ggplot(scr, aes(x = MDS1, y = MDS2)) + geom_point(aes(fill = east_or_west), shape = 21, color = "black", size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(color = east_or_west), level = 0.75, linewidth = 0.7) + scale_fill_manual(values = coast_cols) +
  scale_color_manual(values = coast_cols) + theme_bw() + labs(title = "B", x = "NMDS1", y = "NMDS2", fill = "coast", color = "coast")

p_habitat_by_coast <- ggplot(scr, aes(x = MDS1, y = MDS2)) + geom_point(aes(fill = habitat), shape = 21, color = "black", size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(color = habitat, group = habitat), level = 0.75, linewidth = 0.7) + facet_wrap(~ east_or_west) +
  scale_fill_manual(values = habitat_cols) + scale_color_manual(values = habitat_cols) + theme_bw() +
  labs(title = "C", x = "NMDS1", y = "NMDS2", fill = "habitat", color = "habitat")

figure1_beta <- (p_habitat | p_coast) / p_habitat_by_coast + plot_layout(heights = c(1, 1.1))
figure1_beta

ggsave(file.path(fig_dir, "figure1_beta_diversity.pdf"), figure1_beta, width = 10, height = 8)
ggsave(file.path(fig_dir, "figure1_beta_diversity.png"), figure1_beta, width = 10, height = 8, dpi = 300)

# permanova models
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

message("done. outputs in: ", out_dir)