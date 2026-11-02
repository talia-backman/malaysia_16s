# 04b_rarefy.R
#
# rarefy cleaned phyloseq object to even sequencing depth
# used for unweighted UniFrac and presence/absence analyses
# saves rarefied phyloseq object for downstream use

library(phyloseq)
library(tidyverse)

out_dir <- "./output"
ps_clean <- readRDS(file.path(out_dir, "phyloseq_cleaned_16s.rds"))

# explore rarefaction depth tradeoffs
depths_to_try <- c(500, 525, 750, 1000, 1500, 2000)

message("samples retained at each rarefaction depth:")
for (d in depths_to_try) {
  n <- sum(sample_sums(ps_clean) >= d)
  message("  ", d, " reads: ", n, " samples retained")
}

# rarefy to minimum depth -- 525 retains all 74 filtered samples
rarefy_depth <- 525

set.seed(666) # recorded for reproducibility
ps_rare <- rarefy_even_depth(
  ps_clean,
  sample.size = rarefy_depth,
  rngseed = FALSE,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = TRUE
)

message("rarefied object: ", nsamples(ps_rare), " samples, ", ntaxa(ps_rare), " ASVs")
message("rarefaction depth: ", rarefy_depth)

# confirm all samples are at target depth
stopifnot(all(sample_sums(ps_rare) == rarefy_depth))
message("all samples confirmed at ", rarefy_depth, " reads")

# check coast balance after rarefaction
coast_check <- data.frame(
  coast = sample_data(ps_rare)$east_or_west,
  habitat = sample_data(ps_rare)$habitat
) %>%
  group_by(coast, habitat) %>%
  summarise(n = n(), .groups = "drop")

message("sample distribution after rarefaction:")
print(coast_check)

saveRDS(ps_rare, file.path(out_dir, "phyloseq_rarefied_16s.rds"))
message("rarefied phyloseq saved. script 04b done.")
