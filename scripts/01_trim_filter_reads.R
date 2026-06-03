# 01_trim_filter_reads.R
#
# process pacbio hifi full-length 16s reads
# remove ambiguous reads
# remove primers with cutadapt
# quality filter reads
# retain full-length amplicons (1200-1700 bp)

library(tidyverse)
library(dada2)
library(janitor)

# set threads and avoid nested parallelism
nthreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
message("using nthreads = ", nthreads)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

# load helper functions
source("./scripts/R/functions.R")

# set filtering parameters
bact_27f <- "GAGAGTTTGATCCTGGCTCAG"
bact_1541r <- "AAGGAGGTGATCCAGCCGCA"

min_len <- 1200
max_len <- 1700
max_ee <- 4
trunc_q <- 20
max_n <- 0

# set paths
meta_path <- "./MetaData_eDNA_Malaysia_withQC.tsv"
in_dir <- "fastq"
out_filt_n <- "data/16S/filtN"
out_cutadapt <- "data/16S/cutadapt"
out_filtered <- "data/16S/filtered"

dir.create(out_filt_n, showWarnings = FALSE, recursive = TRUE)
dir.create(out_cutadapt, showWarnings = FALSE, recursive = TRUE)
dir.create(out_filtered, showWarnings = FALSE, recursive = TRUE)

# load metadata
meta <- read_tsv(meta_path, show_col_types = FALSE) %>% janitor::clean_names() %>%
  mutate(sample_name = as.character(sample_id), region = "16S",
         filepaths = file.path(in_dir, paste0(sample_name, "_HiFi.fastq.gz")),
         filt_n_paths = file.path(out_filt_n, paste0(sample_name, "_filtN.fastq.gz")),
         cutadapt_paths = file.path(out_cutadapt, paste0(sample_name, "_cutadapt.fastq.gz")),
         filtered_paths = file.path(out_filtered, paste0(sample_name, "_filtered.fastq.gz")))

if (any(!file.exists(meta$filepaths))) {
  warning("some input FASTQs are missing; subsetting to those that exist.")
  meta <- meta[file.exists(meta$filepaths), ]
}

stopifnot(nrow(meta) > 0)

# remove reads with ambiguous bases
todo <- meta %>% filter(!file.exists(filt_n_paths))

if (nrow(todo) > 0) {
  message("prefiltering ", nrow(todo), " files (maxN = 0)...")
  for (i in seq_len(nrow(todo))) {
    message("[filtN] ", todo$sample_name[i], " : ", basename(todo$filepaths[i]))
    filterAndTrim(fwd = todo$filepaths[i], filt = todo$filt_n_paths[i], maxN = max_n,
                  compress = TRUE, multithread = nthreads, verbose = TRUE)
  }
} else { message("all filtN files already exist; skipping maxN step.") }

# remove primers with cutadapt
if (any(!file.exists(meta$cutadapt_paths))) {
  remove_primers(metadata = meta, amplicon.colname = "region", amplicon = "16S",
                 sampleid.colname = "sample_name", filtN.colname = "filt_n_paths",
                 cutadapt.colname = "cutadapt_paths", fwd_primer = bact_27f,
                 rev_primer = bact_1541r, multithread = max(1, nthreads - 1))
} else { message("all cutadapt outputs already exist; skipping primer removal.") }

# quality filter and retain full-length amplicons
todo <- meta %>% filter(!file.exists(filtered_paths))

if (nrow(todo) > 0) {
  message("qc filtering ", nrow(todo), " files with length bounds [", min_len, ", ", max_len, "] ...")
  
  for (i in seq_len(nrow(todo))) {
    sid <- todo$sample_name[i]
    inf <- todo$cutadapt_paths[i]
    outf <- todo$filtered_paths[i]
    
    message("[filter] ", sid, " : ", basename(inf))
    
    filt_out <- filterAndTrim(fwd = inf, filt = outf, maxN = max_n, maxEE = max_ee,
                              truncQ = trunc_q, minLen = min_len, maxLen = max_len,
                              rm.phix = FALSE, compress = TRUE, multithread = nthreads,
                              verbose = TRUE)
    
    message("[filter] ", sid, " reads.in=", filt_out[1, "reads.in"],
            " reads.out=", filt_out[1, "reads.out"])
  }
} else { message("all qc filtered files already exist; skipping filter step.") }

message("done. outputs in: ", out_filtered)