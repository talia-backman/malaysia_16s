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

# plot retained asv length distribution
asv_lengths <- nchar(taxa_names(ps_clean))

p_length <- ggplot(tibble(length_bp = asv_lengths), aes(length_bp)) + geom_histogram(binwidth = 25) + theme_bw() +
  labs(x = "ASV length (bp)", y = "Number of ASVs")

ggsave(file.path(fig_dir, "asv_length_distribution.pdf"), p_length, width = 7, height = 5)
ggsave(file.path(fig_dir, "asv_length_distribution.png"), p_length, width = 7, height = 5, dpi = 300)

# plot phylum-level composition
ps_phylum_clean <- tax_glom(ps_clean, taxrank = "Phylum")

p_phylum <- plot_bar(ps_phylum_clean, fill = "Phylum") + theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Read abundance")

ggsave(file.path(fig_dir, "phylum_composition_cleaned.pdf"), p_phylum, width = 10, height = 5)
ggsave(file.path(fig_dir, "phylum_composition_cleaned.png"), p_phylum, width = 10, height = 5, dpi = 300)

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

message("done. outputs in: ", out_dir)