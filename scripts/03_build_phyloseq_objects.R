# 03_build_phyloseq_objects.R
#
# build phyloseq objects from dada2 asv table, taxonomy, and metadata
# rename numeric sample ids to safe phyloseq sample names

set.seed(666)

library(tidyverse)
library(janitor)
library(phyloseq)

# set paths
in_meta <- "./MetaData_eDNA_Malaysia_withQC.tsv"
out_dir <- "./output"
seqtab_path <- file.path(out_dir, "seqtab_nochim.rds")
tax_path <- file.path(out_dir, "taxonomy.rds")

# load metadata, asv table, and taxonomy
meta <- read_tsv(in_meta, show_col_types = FALSE) %>% janitor::clean_names() %>%
  mutate(sample_name = as.character(sample_id),
         control = if_else(sample_id == "Blank" | site_code == "Blank" | environmental_sample == "Blank", "Neg", NA_character_))

seqtab_nochim <- readRDS(seqtab_path)
tax <- readRDS(tax_path)

# create safe sample names
rownames(seqtab_nochim) <- as.character(rownames(seqtab_nochim))

sample_name_orig <- rownames(seqtab_nochim)
sample_name_safe <- ifelse(grepl("^[0-9]+$", sample_name_orig), paste0("S", sample_name_orig), sample_name_orig)
stopifnot(anyDuplicated(sample_name_safe) == 0)

rownames(seqtab_nochim) <- sample_name_safe

# build metadata ordered to otu table
meta_ps_all <- meta %>% mutate(sample_name = as.character(sample_name)) %>% filter(sample_name %in% sample_name_orig) %>%
  distinct(sample_name, .keep_all = TRUE) %>%
  mutate(sample_name_orig = sample_name, sample_name = if_else(grepl("^[0-9]+$", sample_name_orig), paste0("S", sample_name_orig), sample_name_orig)) %>%
  slice(match(sample_name_safe, sample_name))

stopifnot(!anyNA(meta_ps_all$sample_name))
stopifnot(identical(meta_ps_all$sample_name, rownames(seqtab_nochim)))

# build phyloseq object with controls
sam_df <- as.data.frame(meta_ps_all, stringsAsFactors = FALSE)
rownames(sam_df) <- sam_df$sample_name
sam_df$sample_name <- NULL

sam <- sample_data(sam_df)
otu <- otu_table(seqtab_nochim, taxa_are_rows = FALSE)
taxmat <- tax_table(tax)

sn_otu <- as.character(sample_names(otu))
sn_sam <- as.character(sample_names(sam))

message("n OTU samples = ", length(sn_otu), "; n sample_data samples = ", length(sn_sam))
message("first 10 OTU sample names: ", paste(head(sn_otu, 10), collapse = ", "))
message("first 10 sample_data names: ", paste(head(sn_sam, 10), collapse = ", "))
message("identical(as.character) = ", identical(sn_otu, sn_sam))

stopifnot(identical(sn_otu, sn_sam))

ps_all <- phyloseq(otu, sam, taxmat)
saveRDS(ps_all, file.path(out_dir, "phyloseq_not_cleaned_with_controls.rds"))

# build phyloseq object without negative controls
meta_ps_no_neg <- meta_ps_all %>% filter(is.na(control) | control != "Neg")
seqtab_no_neg <- seqtab_nochim[meta_ps_no_neg$sample_name, , drop = FALSE]

sam_no_neg_df <- as.data.frame(meta_ps_no_neg, stringsAsFactors = FALSE)
rownames(sam_no_neg_df) <- sam_no_neg_df$sample_name
sam_no_neg_df$sample_name <- NULL

ps_no_neg <- phyloseq(otu_table(seqtab_no_neg, taxa_are_rows = FALSE), sample_data(sam_no_neg_df), taxmat)
saveRDS(ps_no_neg, file.path(out_dir, "phyloseq_not_cleaned_noNeg.rds"))

message("done. outputs in: ", out_dir)