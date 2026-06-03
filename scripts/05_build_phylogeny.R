# 05_build_phylogeny.R
#
# align retained 16s asvs
# build neighbor-joining tree
# add phylogenetic tree to cleaned phyloseq object

set.seed(666)

library(tidyverse)
library(phyloseq)
library(DECIPHER)
library(phylogram)

# load helper functions
source("./scripts/R/functions.R")

# set paths
out_dir <- "./output"
ps_path <- file.path(out_dir, "phyloseq_cleaned_16s.rds")

# load cleaned phyloseq object
ps <- readRDS(ps_path)

# extract asv sequences
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_", seq_along(seqs))

# align sequences
aln <- DECIPHER::AlignSeqs(DNAStringSet(seqs), processors = NULL)
saveRDS(aln, file.path(out_dir, "16s_dna_alignment.RDS"))

# build neighbor-joining tree
dist <- DistanceMatrix(aln, type = "dist", correction = "Jukes-Cantor", verbose = FALSE, processors = NULL)
nj <- TreeLine(myDistMatrix = dist, method = "NJ", cutoff = 0.05, showPlot = FALSE, verbose = FALSE)
saveRDS(nj, file.path(out_dir, "16s_NJtree.RDS"))

# convert tree and restore asv labels
tree <- phylogram::as.phylo(nj)
asv_order <- tree$tip.label %>% str_remove("ASV_") %>% as.numeric()
tree$tip.label <- taxa_names(ps)[asv_order]

# add tree to phyloseq object
ps_w_tree <- phyloseq(otu_table(ps), tax_table(ps), sample_data(ps), phy_tree(tree))
saveRDS(ps_w_tree, file.path(out_dir, "16s_physeq_cleaned_w_tree.RDS"))

message("done. outputs in: ", out_dir)