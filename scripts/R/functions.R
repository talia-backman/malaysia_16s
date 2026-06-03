# scripts/R/functions.R
# Minimal functions used by this repo's pipeline

remove_primers <- function(metadata,
                           amplicon.colname = "region",
                           amplicon = "16S",
                           sampleid.colname = "sample_name",
                           filtN.colname = "filt_n_paths",
                           cutadapt.colname = "cutadapt_paths",
                           fwd_primer = "GTGCCAGCMGCCGCGGTAA",
                           rev_primer = "GGACTACHVGGGTWTCTAAT",
                           multithread = parallel::detectCores() - 1) {

  x <- metadata[metadata[[amplicon.colname]] == amplicon, ]

  fnFs.filtN    <- x[[filtN.colname]]
  fnFs.cutadapt <- x[[cutadapt.colname]]

  FWD <- fwd_primer
  REV <- rev_primer

  REV.RC <- dada2:::rc(REV)

  R1.flags <- paste("-g", FWD, "-a", REV.RC)

  for (i in seq_along(fnFs.filtN)) {
    message("[cutadapt] ", basename(fnFs.filtN[i]))
    system2(
      "cutadapt",
      args = c(
        R1.flags,
        "-n", 2,
        "--minimum-length", "500",
        "--cores", as.character(max(1, multithread)),
        "-o", fnFs.cutadapt[i],
        fnFs.filtN[i]
      )
    )
  }

  invisible(TRUE)
}


qc_filter <- function(metadata,
                      amplicon.colname = "region",
                      amplicon = "16S",
                      input.colname = "cutadapt_paths",
                      output.colname = "filtered_paths",
                      max.ee = 1,
                      trucq = 2,
                      threads = parallel::detectCores() - 1) {

  x <- metadata[metadata[[amplicon.colname]] == amplicon, ]

  input.files  <- x[[input.colname]]
  output.files <- x[[output.colname]]

  out <- dada2::filterAndTrim(
    fwd = input.files,
    filt = output.files,
    maxN = 0,
    maxEE = max.ee,
    truncQ = trucq,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = threads
  )

  saveRDS(out, paste0("./data/", amplicon, "/filtration_stats.RDS"))
  invisible(out)
}

# tidy permanova table
tidy_adonis <- function(x, model_name) {
  broom::tidy(x) %>%
    dplyr::rename(
      sum_of_sqs = dplyr::any_of(c("sumsqs", "SumOfSqs")),
      r2 = dplyr::any_of(c("r2", "R2"))
    ) %>%
    dplyr::select(term, df, sum_of_sqs, r2, statistic, p.value) %>%
    dplyr::mutate(model = model_name)
}

# betadisperser
run_betadisp <- function(dist_obj, group_vec, model_name, permutations = 999) {
  if (length(group_vec) != attr(dist_obj, "Size")) {
    stop("length of group_vec must match the number of samples in dist_obj")
  }
  
  group_vec <- as.factor(group_vec)
  
  bd <- vegan::betadisper(dist_obj, group = group_vec)
  pt <- vegan::permutest(bd, permutations = permutations)
  
  list(
    bd = bd,
    pt = pt,
    summary = tibble::tibble(
      model = model_name,
      F = unname(pt$tab[1, "F"]),
      p = unname(pt$tab[1, "Pr(>F)"])
    ),
    distances = tibble::tibble(
      model = model_name,
      sample_name = names(bd$distances),
      group = as.character(bd$group),
      dispersion = as.numeric(bd$distances)
    )
  )
}

# mantel tests
run_mantel <- function(ids, label, dist_mat, geo_mat, dist_name,
                       permutations = 999, out_dir = NULL) {
  
  ids <- intersect(ids, rownames(dist_mat))
  
  if (length(ids) < 3) {
    warning("fewer than 3 samples for ", label, " in ", dist_name, "; skipping")
    return(NULL)
  }
  
  res <- vegan::mantel(
    stats::as.dist(dist_mat[ids, ids, drop = FALSE]),
    stats::as.dist(geo_mat[ids, ids, drop = FALSE]),
    method = "spearman",
    permutations = permutations
  )
  
  if (!is.null(out_dir)) {
    capture.output(
      res,
      file = file.path(out_dir,
                       paste0("mantel_", dist_name, "_vs_ocean_", label, ".txt"))
    )
  }
  
  print(res)
  
  tibble::tibble(
    metric = dist_name,
    test = label,
    n_samples = length(ids),
    statistic_r = unname(res$statistic),
    p = res$signif
  )
}

# turn two distance matrices into a pairwise table for plotting
dist_to_df <- function(comm_mat, geo_mat, metric_name, meta_df) {
  # check matrix dimensions and names
  if (!all(rownames(comm_mat) %in% rownames(geo_mat))) {
    stop("not all rownames in comm_mat are present in geo_mat")
  }
  
  if (!all(colnames(comm_mat) %in% colnames(geo_mat))) {
    stop("not all colnames in comm_mat are present in geo_mat")
  }
  
  required_cols <- c("sample_name", "east_or_west", "habitat")
  if (!all(required_cols %in% colnames(meta_df))) {
    stop("meta_df must contain: sample_name, east_or_west, habitat")
  }
  
  comm_long <- as.data.frame(as.table(comm_mat)) %>%
    dplyr::rename(sample_1 = Var1, sample_2 = Var2, community_distance = Freq)
  
  geo_long <- as.data.frame(as.table(geo_mat)) %>%
    dplyr::rename(sample_1 = Var1, sample_2 = Var2, ocean_distance_km = Freq)
  
  pair_df <- comm_long %>%
    dplyr::left_join(geo_long, by = c("sample_1", "sample_2")) %>%
    dplyr::filter(sample_1 != sample_2) %>%
    dplyr::mutate(
      pair_id = purrr::map2_chr(sample_1, sample_2, ~ paste(sort(c(.x, .y)), collapse = "__"))
    ) %>%
    dplyr::distinct(pair_id, .keep_all = TRUE) %>%
    dplyr::select(-pair_id)
  
  meta_small <- meta_df %>%
    tibble::as_tibble() %>%
    dplyr::select(sample_name, east_or_west, habitat)
  
  pair_df <- pair_df %>%
    dplyr::left_join(meta_small, by = c("sample_1" = "sample_name")) %>%
    dplyr::rename(east_or_west_1 = east_or_west, habitat_1 = habitat) %>%
    dplyr::left_join(meta_small, by = c("sample_2" = "sample_name")) %>%
    dplyr::rename(east_or_west_2 = east_or_west, habitat_2 = habitat) %>%
    dplyr::mutate(metric = metric_name)
  
  pair_df
}