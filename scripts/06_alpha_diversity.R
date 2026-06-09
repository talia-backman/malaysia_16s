# 06_alpha_diversity.R
#
# calculate alpha diversity
# test habitat, coast, and habitat-by-coast effects using site as a random effect
# save figure s2

set.seed(666)

library(tidyverse)
library(phyloseq)
library(lmerTest)
library(broom)
library(emmeans)
library(patchwork)

# load helper functions
source("./scripts/R/functions.R")

# set paths
ps_path <- "./output/16s_physeq_cleaned_w_tree.RDS"
out_dir <- "./output/alpha_diversity"
fig_dir <- "./figures/supplemental/figure_s2_alpha_diversity"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# load phyloseq object
ps <- readRDS(ps_path)

# prepare metadata
meta <- as(sample_data(ps), "data.frame") %>% tibble::rownames_to_column("sample_name") %>% janitor::clean_names() %>%
  mutate(east_west = str_to_lower(trimws(east_or_west)), habitat = str_to_lower(trimws(habitat)), site_name = trimws(site_name))

message("samples in ps: ", nsamples(ps), " | asvs: ", ntaxa(ps))
print(meta %>% count(east_west, habitat, name = "n") %>% arrange(desc(n)))

# calculate alpha diversity
alpha <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson")) %>% rownames_to_column("sample_name") %>%
  mutate(sample_name = stringr::str_remove(sample_name, "^X")) %>% left_join(meta, by = "sample_name")

write_csv(alpha, file.path(out_dir, "alpha_diversity_16s.csv"))

# plot alpha diversity
habitat_cols <- c(coral = "#440154FF", mangrove = "#21908CFF", seagrass = "#FDE725FF")

p_rich <- alpha %>% ggplot(aes(x = east_west, y = Observed)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
  geom_jitter(aes(fill = habitat), shape = 21, color = "black", width = 0.15, height = 0, alpha = 0.75, size = 2.3, stroke = 0.3) +
  facet_wrap(~ habitat, scales = "free_y") + scale_fill_manual(values = habitat_cols) +
  theme_bw() + labs(x = "coast", y = "asv richness", fill = "habitat")

p_shan <- alpha %>% ggplot(aes(x = east_west, y = Shannon)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
  geom_jitter(aes(fill = habitat), shape = 21, color = "black", width = 0.15, height = 0, alpha = 0.75, size = 2.3, stroke = 0.3) +
  facet_wrap(~ habitat, scales = "free_y") + scale_fill_manual(values = habitat_cols) +
  theme_bw() + labs(x = "coast", y = "shannon diversity", fill = "habitat")

p_alpha <- (p_rich / p_shan) + plot_layout(guides = "collect") &
  theme(text = element_text(size = 14), plot.title = element_text(face = "bold"), legend.position = "none")

p_alpha

ggsave(file.path(fig_dir, "figure_s2_alpha_diversity.pdf"), p_alpha, width = 10, height = 8)
ggsave(file.path(fig_dir, "figure_s2_alpha_diversity.png"), p_alpha, width = 10, height = 8, dpi = 300)

# fit mixed models with site as random intercept
alpha2 <- alpha %>% filter(!is.na(east_west), !is.na(habitat), !is.na(site_name)) %>%
  mutate(east_west = factor(east_west), habitat = factor(habitat), site_name = factor(site_name))

mod_rich <- lmer(Observed ~ east_west * habitat + (1 | site_name), data = alpha2)
mod_shan <- lmer(Shannon ~ east_west * habitat + (1 | site_name), data = alpha2)

sink(file.path(out_dir, "alpha_models_summary.txt"))
cat("alpha diversity mixed models\n\n")
cat("\nrichness observed\n")
print(summary(mod_rich))
cat("\nanova type iii\n")
print(anova(mod_rich, type = 3))
cat("\nshannon diversity\n")
print(summary(mod_shan))
cat("\nanova type iii\n")
print(anova(mod_shan, type = 3))
sink()

message("done. outputs in: ", out_dir)