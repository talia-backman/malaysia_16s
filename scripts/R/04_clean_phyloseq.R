# 04_clean_phyloseq.R
#
# clean phyloseq object
# remove negative controls
# remove organelles
# tidy taxonomy labels
# save cleaned phyloseq objects

library(phyloseq)
library(tidyverse)
library(janitor)
library(vegan)
library(patchwork)

# set paths
out_dir <- "./output"
fig_dir <- "./figures/supplemental/figure_s1_dataset_characteristics"
ps_path <- file.path(out_dir, "phyloseq_not_cleaned_with_controls.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# load phyloseq object
ps <- readRDS(ps_path)

message("starting samples: ", nsamples(ps))
message("starting ASVs: ", ntaxa(ps))

# remove negative controls
sd_df <- as(sample_data(ps), "data.frame")
stopifnot("control" %in% colnames(sd_df))

samples_no_neg <- rownames(sd_df)[is.na(sd_df$control) | sd_df$control != "Neg"]
ps_no_neg <- prune_samples(samples_no_neg, ps)
ps_no_neg <- prune_taxa(taxa_sums(ps_no_neg) > 0, ps_no_neg)

# remove organelles and taxa without phylum assignment
tax_df <- as.data.frame(tax_table(ps_no_neg))

keep_fam <- if ("Family" %in% colnames(tax_df)) { is.na(tax_df$Family) | tax_df$Family != "Mitochondria" } else { rep(TRUE, nrow(tax_df)) }
keep_ord <- if ("Order" %in% colnames(tax_df)) { is.na(tax_df$Order) | tax_df$Order != "Chloroplast" } else { rep(TRUE, nrow(tax_df)) }

if ("Phylum" %in% colnames(tax_df)) {
  keep_phy <- !is.na(tax_df$Phylum)
} else { stop("no Phylum column found in tax_table(ps_no_neg). check taxonomy ranks.") }

keep_taxa <- keep_phy & keep_fam & keep_ord
ps_no_org <- prune_taxa(keep_taxa, ps_no_neg)
ps_no_org <- prune_taxa(taxa_sums(ps_no_org) > 0, ps_no_org)

# build table S1 sample metadata and post-filtering read counts
# apply the final retained bacterial taxon set to the original object so that
# environmental samples and negative controls are summarized consistently

final_bacterial_taxa <- taxa_names(ps_no_org)

ps_qc_all <- prune_taxa(final_bacterial_taxa, ps)

# retained bacterial reads after quality control and taxonomic filtering
retained_reads <- sample_sums(ps_qc_all)

table_s1 <- as(sample_data(ps_qc_all), "data.frame") %>%
  tibble::rownames_to_column("phyloseq_sample_name") %>%
  janitor::clean_names() %>%
  mutate(
    retained_reads = as.numeric(retained_reads[phyloseq_sample_name]),
    
    sample_type = case_when(
      control == "Neg" ~ "Negative control",
      environmental_sample == "water" ~ "Environmental sample",
      TRUE ~ as.character(environmental_sample)
    ),
    
    final_analysis_status = case_when(
      control == "Neg" ~ "Negative control",
      environmental_sample == "water" & retained_reads == 0 ~
        "Excluded: zero bacterial reads after QC/taxonomic filtering",
      environmental_sample == "water" & retained_reads > 0 ~ "Retained",
      TRUE ~ NA_character_
    ),
    
    sample = case_when(
      environmental_sample == "water" ~ paste0("S", sample_id),
      TRUE ~ phyloseq_sample_name
    )
  ) %>%
  transmute(
    Sample = sample,
    `Sample type` = sample_type,
    `Site name` = site_name,
    `Site code` = site_code,
    `Site replicate` = site_replicate,
    Habitat = habitat,
    Coast = east_or_west,
    `Protected area` = protected_area,
    Latitude = latitude,
    Longitude = longitude,
    `Sampling date` = as.Date(sampling_date),
    `Volume (mL)` = volume_ml,
    `Retained reads` = retained_reads,
    `Final analysis status` = final_analysis_status
  )

write_csv(
  table_s1,
  file.path(out_dir, "table_s1_sample_metadata_qc.csv")
)

# print checks used in manuscript
env_qc <- table_s1 %>%
  filter(`Sample type` == "Environmental sample")

neg_qc <- table_s1 %>%
  filter(`Sample type` == "Negative control")

message("environmental samples total: ", nrow(env_qc))
message("environmental samples retained: ",
        sum(env_qc$`Final analysis status` == "Retained", na.rm = TRUE))
message("environmental samples excluded: ",
        sum(str_detect(env_qc$`Final analysis status`, "^Excluded"), na.rm = TRUE))

message(
  "retained environmental sequencing depth min/median/max: ",
  min(env_qc$`Retained reads`[env_qc$`Retained reads` > 0]), " / ",
  median(env_qc$`Retained reads`[env_qc$`Retained reads` > 0]), " / ",
  max(env_qc$`Retained reads`[env_qc$`Retained reads` > 0])
)

message("negative controls: ", nrow(neg_qc))
message(
  "negative-control retained reads: ",
  paste(neg_qc$`Retained reads`, collapse = ", ")
)

# tidy species labels
tax_df_clean <- as.data.frame(tax_table(ps_no_org))

if ("Species" %in% colnames(tax_df_clean)) {
  tax_df_clean$Species <- as.character(tax_df_clean$Species)
  tax_df_clean$Species[tax_df_clean$Species == "unclassified"] <- NA
  tax_df_clean$Species[is.na(tax_df_clean$Species)] <- "sp."
  tax_table(ps_no_org) <- tax_table(as.matrix(tax_df_clean))
}

# save cleaned objects
ps_clean <- ps_no_org

saveRDS(ps_no_neg, file.path(out_dir, "phyloseq_noNeg.rds"))
saveRDS(ps_clean, file.path(out_dir, "phyloseq_cleaned_16s.rds"))

# prepare non-empty environmental samples for sequencing-depth plots
ps_plot <- prune_samples(sample_sums(ps_clean) > 0, ps_clean)
ps_plot <- prune_taxa(taxa_sums(ps_plot) > 0, ps_plot)

# summarize sequencing depth
sequencing_depth <- sample_sums(ps_plot)
message("sequencing depth min: ", min(sequencing_depth))
message("sequencing depth median: ", median(sequencing_depth))
message("sequencing depth max: ", max(sequencing_depth))

# plot retained asv length distribution
asv_lengths <- nchar(taxa_names(ps_clean))

p_length <- ggplot(tibble(length_bp = asv_lengths), aes(length_bp)) + geom_histogram(binwidth = 25) + theme_bw() +
  labs(x = "ASV length (bp)", y = "Number of ASVs")

ggsave(file.path(fig_dir, "asv_length_distribution.pdf"), p_length, width = 7, height = 5)
ggsave(file.path(fig_dir, "asv_length_distribution.png"), p_length, width = 7, height = 5, dpi = 300)

# plot sample rarefaction curves
otu_mat <- as(otu_table(ps_plot), "matrix")
if (taxa_are_rows(ps_plot)) otu_mat <- t(otu_mat)

rare_df <- vegan::rarecurve(otu_mat, step = 100, label = FALSE, tidy = TRUE) %>%
  dplyr::rename(sample_name = Site, sequencing_depth = Sample, observed_asvs = Species)

rare_meta <- as(sample_data(ps_plot), "data.frame") %>% rownames_to_column("sample_name") %>% janitor::clean_names()

p_rarefaction <- ggplot(rare_df, aes(x = sequencing_depth, y = observed_asvs, group = sample_name)) +
  geom_line(alpha = 0.6, linewidth = 0.5) + theme_bw() +
  labs(x = "Sequencing depth", y = "Observed ASVs")

ggsave(file.path(fig_dir, "sample_rarefaction_curves.pdf"), p_rarefaction, width = 7, height = 5)
ggsave(file.path(fig_dir, "sample_rarefaction_curves.png"), p_rarefaction, width = 7, height = 5, dpi = 300)

# plot phylum-level composition
# plot phylum-level composition
ps_phylum_clean <- tax_glom(ps_clean, taxrank = "Phylum")

phylum_cols <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                 "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100",
                 "#6699CC", "#AA4466", "#4477AA", "#228833", "#EE6677",
                 "#EE8866", "#BBBBBB")

p_phylum <- plot_bar(ps_phylum_clean, fill = "Phylum") + theme_bw() +
  scale_fill_manual(values = phylum_cols) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Read abundance")

ggsave(file.path(fig_dir, "phylum_composition_cleaned.pdf"), p_phylum, width = 10, height = 5)
ggsave(file.path(fig_dir, "phylum_composition_cleaned.png"), p_phylum, width = 10, height = 5, dpi = 300)

# combine figure S1
figure_s1 <- p_length + p_rarefaction + p_phylum + plot_annotation(tag_levels = "A")

ggsave(file.path(fig_dir, "figure_s1_dataset_characteristics.pdf"), figure_s1, width = 16, height = 5)
ggsave(file.path(fig_dir, "figure_s1_dataset_characteristics.png"), figure_s1, width = 16, height = 5, dpi = 300)

# write cleaning summary
summary_file <- file.path(out_dir, "physeq_cleaning_summary.txt")
con <- file(summary_file, open = "wt")
on.exit(close(con), add = TRUE)

writeLines(c(
  "cleaning summary",
  paste0("original ps:        ", nsamples(ps), " samples; ", ntaxa(ps), " ASVs"),
  paste0("after drop neg:     ", nsamples(ps_no_neg), " samples; ", ntaxa(ps_no_neg), " ASVs"),
  paste0("after organelle rm: ", nsamples(ps_no_org), " samples; ", ntaxa(ps_no_org), " ASVs"),
  paste0("final cleaned:      ", nsamples(ps_clean), " samples; ", ntaxa(ps_clean), " ASVs")
), con = con)

message("script 4 done. outputs in: ", out_dir)