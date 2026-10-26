# Shared helpers for FEMA analytic sample tables (iid/eid -> model_sample).

bmi_designmat_model_dirs <- function(design_tag) {
  c(
    sprintf("designmat1f_%sSageSpdsIncEducPCsScanSoft_nodemean_females", design_tag),
    sprintf("designmat1m_%sSageSpdsIncEducPCsScanSoft_nodemean_males", design_tag)
  )
}

bmi_discover_designmat_models <- function(results_dir) {
  if (!dir.exists(results_dir)) {
    return(character())
  }
  models <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
  models[grepl("^designmat1[fm]_", models)]
}

bmi_model_sample_paths <- function(sample_tables_dir, design_tag, modality = "rni") {
  models <- bmi_designmat_model_dirs(design_tag)
  file.path(sample_tables_dir, models, paste0("model_sample_", modality, ".csv"))
}

bmi_read_iid_eid <- function(path) {
  stopifnot(file.exists(path))
  idlist <- read.delim(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  stopifnot(ncol(idlist) >= 2)
  names(idlist)[1:2] <- c("participant_id", "session_id")
  idlist$participant_id <- as.character(idlist$participant_id)
  idlist$session_id <- as.character(idlist$session_id)
  idlist[, c("participant_id", "session_id"), drop = FALSE]
}

bmi_visit_keys <- function(df) {
  df %>%
    dplyr::transmute(
      participant_id = as.character(participant_id),
      session_id = as.character(session_id)
    ) %>%
    dplyr::distinct()
}

bmi_load_model_sample_visits <- function(sample_tables_dir, design_tag, modality = "rni") {
  paths <- bmi_model_sample_paths(sample_tables_dir, design_tag, modality)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Missing model_sample file(s):\n  ",
      paste(missing, collapse = "\n  "),
      call. = FALSE
    )
  }
  dfs <- lapply(paths, function(p) {
    readr::read_csv(
      p,
      show_col_types = FALSE,
      col_select = c("participant_id", "session_id"),
      col_types = readr::cols(
        participant_id = readr::col_character(),
        session_id = readr::col_character()
      )
    )
  })
  bmi_visit_keys(dplyr::bind_rows(dfs))
}

bmi_compare_visit_sets <- function(keys_a, keys_b, label_a = "A", label_b = "B") {
  only_a <- dplyr::anti_join(keys_a, keys_b, by = c("participant_id", "session_id"))
  only_b <- dplyr::anti_join(keys_b, keys_a, by = c("participant_id", "session_id"))
  list(
    n_a = nrow(keys_a),
    n_b = nrow(keys_b),
    n_intersection = nrow(dplyr::inner_join(keys_a, keys_b, by = c("participant_id", "session_id"))),
    n_only_a = nrow(only_a),
    n_only_b = nrow(only_b),
    only_a = only_a,
    only_b = only_b,
    match = nrow(only_a) == 0 && nrow(only_b) == 0 && nrow(keys_a) == nrow(keys_b)
  )
}

bmi_print_visit_set_comparison <- function(cmp, label_a, label_b) {
  cat(sprintf(
    "\n%s vs %s visit comparison:\n  %s: %s visits\n  %s: %s visits\n  intersection: %s\n",
    label_a, label_b,
    label_a, format(cmp$n_a, big.mark = ","),
    label_b, format(cmp$n_b, big.mark = ","),
    format(cmp$n_intersection, big.mark = ",")
  ))
  if (cmp$match) {
    cat("  MATCH: identical participant_id + session_id sets.\n")
  } else {
    cat(sprintf(
      "  MISMATCH: only in %s: %s | only in %s: %s\n",
      label_a, format(cmp$n_only_a, big.mark = ","),
      label_b, format(cmp$n_only_b, big.mark = ",")
    ))
    if (cmp$n_only_a > 0 && cmp$n_only_a <= 10) {
      cat("  Only in", label_a, ":\n")
      print(cmp$only_a)
    }
    if (cmp$n_only_b > 0 && cmp$n_only_b <= 10) {
      cat("  Only in", label_b, ":\n")
      print(cmp$only_b)
    }
  }
  invisible(cmp)
}

bmi_paths_for_pheno <- function(pheno_var, repo_root = Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")) {
  old_pheno <- Sys.getenv("BMI_PHENO_VARIABLE", unset = NA)
  old_analysis <- Sys.getenv("BMI_IMAGING_ANALYSIS_DIR", unset = NA)
  on.exit({
    if (!is.na(old_pheno)) Sys.setenv(BMI_PHENO_VARIABLE = old_pheno) else Sys.unsetenv("BMI_PHENO_VARIABLE")
    if (!is.na(old_analysis)) Sys.setenv(BMI_IMAGING_ANALYSIS_DIR = old_analysis) else Sys.unsetenv("BMI_IMAGING_ANALYSIS_DIR")
  }, add = TRUE)

  Sys.setenv(BMI_PHENO_VARIABLE = pheno_var)
  source(file.path(normalizePath(repo_root, mustWork = FALSE), "scripts", "config_project.R"), local = FALSE)
  profile <- bmi_pheno_profile(pheno_var)
  Sys.setenv(BMI_IMAGING_ANALYSIS_DIR = profile$analysis_dir)
  bmi_project_paths()
}
