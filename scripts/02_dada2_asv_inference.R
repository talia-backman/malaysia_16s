# 02_dada2_asv_inference.R
#
# infer asvs from filtered pacbio hifi 16s reads
# remove chimeras
# assign taxonomy with silva

set.seed(666)

library(tidyverse)
library(dada2)
library(janitor)

# set threads and avoid nested parallelism
nthreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
message("using nthreads = ", nthreads)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

# set paths
in_meta <- "./MetaData_eDNA_Malaysia_withQC.tsv"
filt_dir <- "./data/16S/filtered"
out_dir <- "./output"
tax_db <- "./taxonomy/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
species_db <- "./taxonomy/silva_v138.2_assignSpecies.fa.gz"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# set dada2 parameters
min_len <- 1200
max_len <- 1700
pooling <- FALSE # manuscript analysis

# load metadata and filtered read paths
meta <- read_tsv(in_meta, show_col_types = FALSE) %>% janitor::clean_names() %>%
  mutate(sample_name = as.character(sample_id),
         filtered_path = file.path(filt_dir, paste0(sample_name, "_filtered.fastq.gz")))

if (any(!file.exists(meta$filtered_path))) {
  warning("some filtered FASTQs are missing; subsetting to those that exist.")
  meta <- meta[file.exists(meta$filtered_path), ]
}

stopifnot(nrow(meta) > 0)

if (!file.exists(tax_db)) {
  stop("tax_db not found: ", tax_db, "\n", "you need a dada2-formatted 16s training set at this path.")
}

fn <- meta$filtered_path
names(fn) <- meta$sample_name
message("number of FASTQs in analysis: ", length(fn))

# learn error model
err_rds <- file.path(out_dir, "err_model.rds")

if (!file.exists(err_rds)) {
  message("learning error model from filtered reads...")
  learn_n <- min(12, length(fn))
  learn_on <- sample(fn, learn_n)
  message("learning errors on ", learn_n, " FASTQs; nbases = 2e7; randomize = TRUE")
  err <- learnErrors(learn_on, multithread = max(1, nthreads - 1), verbose = TRUE, nbases = 2e7, randomize = TRUE)
  saveRDS(err, err_rds)
} else { err <- readRDS(err_rds) }

pdf(file.path(out_dir, "error_model.pdf"))
plotErrors(err, nominalQ = TRUE)
dev.off()

# infer asvs sample-by-sample with checkpoints
message("dereplication...")
message("running dada...")

dada_dir <- file.path(out_dir, "dada_per_sample")
dir.create(dada_dir, showWarnings = FALSE, recursive = TRUE)

getN <- function(x) sum(getUniques(x))
input_reads <- setNames(rep(NA_integer_, length(fn)), names(fn))
denoised_reads <- setNames(rep(NA_integer_, length(fn)), names(fn))

dada_list <- vector("list", length(fn))
names(dada_list) <- names(fn)

for (s in names(fn)) {
  message("sample: ", s)
  
  dada_s_rds <- file.path(dada_dir, paste0(s, "_dada.rds"))
  track_s_rds <- file.path(dada_dir, paste0(s, "_track.rds"))
  
  if (file.exists(dada_s_rds) && file.exists(track_s_rds)) {
    dada_s <- readRDS(dada_s_rds)
    track_s <- readRDS(track_s_rds)
    input_reads[s] <- as.integer(track_s$input)
    denoised_reads[s] <- as.integer(track_s$denoised)
    dada_list[[s]] <- dada_s
    rm(dada_s, track_s)
    gc()
    next
  }
  
  derep_s <- derepFastq(fn[[s]], verbose = TRUE)
  input_reads[s] <- sum(derep_s$uniques)
  dada_s <- dada(derep_s, err = err, multithread = max(1, nthreads - 1), pool = pooling)
  denoised_reads[s] <- getN(dada_s)
  
  saveRDS(dada_s, dada_s_rds)
  saveRDS(list(input = input_reads[s], denoised = denoised_reads[s]), track_s_rds)
  
  dada_list[[s]] <- dada_s
  rm(derep_s, dada_s)
  gc()
}

# build sequence table and remove chimeras
message("building sequence table...")
seqtab <- makeSequenceTable(dada_list)

seqlens <- nchar(colnames(seqtab))
keep <- seqlens >= min_len & seqlens <= max_len
seqtab_lenf <- seqtab[, keep, drop = FALSE]
message("length filter kept ", sum(keep), " ASVs out of ", length(keep))

message("removing chimeras...")
seqtab_nochim <- removeBimeraDenovo(seqtab_lenf, method = "consensus", multithread = max(1, nthreads - 1), verbose = TRUE)
nonchim_reads <- rowSums(seqtab_nochim)

saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_nochim.rds"))

track <- tibble(sample_name = meta$sample_name, input = as.integer(input_reads[meta$sample_name]),
                denoised = as.integer(denoised_reads[meta$sample_name]),
                nonchim = as.integer(nonchim_reads[meta$sample_name]))

write_csv(track, file.path(out_dir, "read_tracking_dada2.csv"))

# assign taxonomy
tax_rds <- file.path(out_dir, "taxonomy.rds")

if (!file.exists(tax_rds)) {
  message("assigning taxonomy...")
  tax <- assignTaxonomy(seqtab_nochim, tax_db, multithread = max(1, nthreads - 1), minBoot = 50)
  
  if (file.exists(species_db)) {
    tax <- addSpecies(tax, species_db, allowMultiple = TRUE)
  } else { warning("species_db not found; skipping addSpecies(): ", species_db) }
  
  saveRDS(tax, tax_rds)
} else { tax <- readRDS(tax_rds) }

message("done. outputs in: ", out_dir)