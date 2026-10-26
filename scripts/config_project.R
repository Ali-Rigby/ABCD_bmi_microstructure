# Unified BMI imaging pipeline: shared s0/s1 + separate bmip / perc_median analyses.
#
# Environment variables (BMI_IMAGING_* prefix):
#   BMI_IMAGING_REPO_ROOT       - this repository
#   BMI_IMAGING_CEPH_ROOT       - .../2026_bmi_imaging_analysis
#   BMI_PHENO_VARIABLE          - bmip | perc_median (required for profile scripts s2-s9)
#   BMI_IMAGING_ANALYSIS_DIR    - bmip | perc_median (default from pheno profile)
#   BMI_IMAGING_SHARED_DATAFILES_DIR - optional override for shared/datafiles
#   BMI_IMAGING_DATAFILES_DIR   - optional override for <profile>/datafiles
#   BMI_IMAGING_DESIGNMATS_DIR / DESIGNMATS_SUBDIR
#   BMI_IMAGING_RESULTS_ROOT / RESULTS_RUN
#   BMI_IMAGING_RUN_TAG            - e.g. 2026_06_29 (auto from date via env/run_date.csh)
#   BMI_IMAGING_AUTO_RUN_DATE      - 1 = fill subdirs from RUN_TAG when unset (default)
#   ABCD_TABULATED_ROOT

bmi_env_datafile_name <- function(var, default) {
  setting <- trimws(Sys.getenv(var, unset = ""))
  if (nzchar(setting)) {
    return(setting)
  }
  default
}

bmi_default_repo_root <- function() {
  "/home/arigby/projects/ABCD_bmi_microstructure"
}

bmi_default_ceph_root <- function() {
  "/space/ceph/1/ABCD/users/arigby/2026_bmi_imaging_analysis"
}

bmi_default_abcd_tabulated <- function() {
  "/space/ceph/1/ABCD/data/release_archive/6_1/tabulated"
}

if (!nzchar(Sys.getenv("BMI_IMAGING_CEPH_ROOT", unset = ""))) {
  Sys.setenv(BMI_IMAGING_CEPH_ROOT = bmi_default_ceph_root())
}
if (!nzchar(Sys.getenv("ABCD_TABULATED_ROOT", unset = ""))) {
  Sys.setenv(ABCD_TABULATED_ROOT = bmi_default_abcd_tabulated())
}

bmi_resolve_under_ceph <- function(ceph_root, setting) {
  if (!nzchar(setting)) {
    return("")
  }
  if (grepl("^/", setting)) {
    return(normalizePath(setting, mustWork = FALSE))
  }
  normalizePath(file.path(ceph_root, setting), mustWork = FALSE)
}

bmi_shared_datafiles_dir <- function(ceph_root) {
  setting <- Sys.getenv("BMI_IMAGING_SHARED_DATAFILES_DIR", unset = "")
  if (!nzchar(setting)) {
    return(normalizePath(file.path(ceph_root, "shared", "datafiles"), mustWork = FALSE))
  }
  bmi_resolve_under_ceph(ceph_root, setting)
}

bmi_resolve_analysis_base <- function(ceph_root) {
  setting <- Sys.getenv("BMI_IMAGING_ANALYSIS_DIR", unset = "")
  if (!nzchar(setting)) {
    stop(
      "Set BMI_IMAGING_ANALYSIS_DIR or BMI_PHENO_VARIABLE (via env/run_<pheno>.csh) before sourcing config for profile scripts.",
      call. = FALSE
    )
  }
  bmi_resolve_under_ceph(ceph_root, setting)
}

bmi_resolve_profile_datafiles_dir <- function(analysis_base, ceph_root) {
  setting <- Sys.getenv("BMI_IMAGING_DATAFILES_DIR", unset = "")
  if (!nzchar(setting)) {
    return(normalizePath(file.path(analysis_base, "datafiles"), mustWork = FALSE))
  }
  bmi_resolve_under_ceph(ceph_root, setting)
}

bmi_resolve_designmats_dir <- function(analysis_base, ceph_root) {
  setting <- Sys.getenv("BMI_IMAGING_DESIGNMATS_DIR", unset = "")
  if (nzchar(setting)) {
    return(bmi_resolve_under_ceph(ceph_root, setting))
  }
  base_dm <- file.path(analysis_base, "designmats")
  sub <- Sys.getenv("BMI_IMAGING_DESIGNMATS_SUBDIR", unset = "")
  if (nzchar(sub)) {
    return(normalizePath(file.path(base_dm, sub), mustWork = FALSE))
  }
  normalizePath(base_dm, mustWork = FALSE)
}

bmi_resolve_results_root <- function(analysis_base, ceph_root) {
  setting <- Sys.getenv("BMI_IMAGING_RESULTS_ROOT", unset = "")
  if (!nzchar(setting)) {
    return(normalizePath(file.path(analysis_base, "results"), mustWork = FALSE))
  }
  bmi_resolve_under_ceph(ceph_root, setting)
}

bmi_auto_run_date_enabled <- function() {
  setting <- tolower(trimws(Sys.getenv("BMI_IMAGING_AUTO_RUN_DATE", unset = "1")))
  setting %in% c("1", "true", "yes", "on")
}

bmi_run_tag <- function() {
  tag <- trimws(Sys.getenv("BMI_IMAGING_RUN_TAG", unset = ""))
  if (!nzchar(tag)) {
    tag <- format(Sys.Date(), "%Y_%m_%d")
    Sys.setenv(BMI_IMAGING_RUN_TAG = tag)
  }
  tag
}

bmi_ensure_run_tag_env <- function() {
  if (!bmi_auto_run_date_enabled()) {
    return(invisible(NULL))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_DESIGNMATS_DIR", unset = ""))) {
    return(invisible(NULL))
  }
  tag <- bmi_run_tag()
  if (!nzchar(Sys.getenv("BMI_IMAGING_DESIGNMATS_SUBDIR", unset = ""))) {
    Sys.setenv(BMI_IMAGING_DESIGNMATS_SUBDIR = paste0("designmats_", tag))
  }
  if (!nzchar(Sys.getenv("BMI_IMAGING_RESULTS_RUN", unset = ""))) {
    Sys.setenv(BMI_IMAGING_RESULTS_RUN = paste0("results_", tag))
  }
  invisible(NULL)
}

bmi_project_paths_core <- function() {
  repo <- Sys.getenv("BMI_IMAGING_REPO_ROOT", unset = bmi_default_repo_root())
  ceph <- Sys.getenv("BMI_IMAGING_CEPH_ROOT", unset = bmi_default_ceph_root())
  tab <- Sys.getenv("ABCD_TABULATED_ROOT", unset = bmi_default_abcd_tabulated())
  repo <- normalizePath(repo, mustWork = FALSE)
  shared_datafiles <- bmi_shared_datafiles_dir(ceph)

  list(
    repo_root = repo,
    abcd_tabulated = tab,
    ceph_project = ceph,
    shared_datafiles = shared_datafiles,
    base_datafile_pheno_qc = bmi_env_datafile_name(
      "BMI_IMAGING_BASE_DATAFILE_PHENO_QC",
      "datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_removed.txt"
    ),
    base_datafile_imaging = bmi_env_datafile_name(
      "BMI_IMAGING_BASE_DATAFILE_IMAGING",
      "datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_dmri_incl.txt"
    ),
    cmig_tools_r = "/usr/pubsw/packages/cmig/cmig_tools/cmig_tools_utils/r",
    cdc_ref_rda = file.path(repo, "data", "cdc_ref_data.rda"),
    cdcanthro_r = file.path(repo, "utilities", "r_utils", "cdcanthro.R")
  )
}

bmi_shared_paths <- function() {
  paths <- bmi_project_paths_core()
  paths$final_datafiles <- paths$shared_datafiles
  paths
}

bmi_age_discrepancy_csv_path <- function(paths) {
  f <- "age_scan_anthro_discrepancy_discrepant_only.csv"
  shared <- paths$shared_datafiles %||% paths$final_datafiles
  p <- file.path(shared, f)
  if (file.exists(p)) {
    return(p)
  }
  stop("Age discrepancy CSV not found under shared datafiles: ", p, call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

bmi_assert_file <- function(path, what = "File") {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", what, path), call. = FALSE)
  }
  invisible(path)
}

bmi_assert_dir <- function(path, what = "Directory") {
  if (!dir.exists(path)) {
    stop(sprintf("%s not found or not a directory: %s", what, path), call. = FALSE)
  }
  invisible(path)
}

bmi_ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  invisible(path)
}

bmi_assert_cols <- function(df, cols, what = "data.frame") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(sprintf("%s missing required columns: %s", what, paste(miss, collapse = ", ")), call. = FALSE)
  }
  invisible(df)
}

bmi_assert_positive_nrow <- function(df, what = "data.frame") {
  if (!nrow(df)) {
    stop(sprintf("%s has zero rows", what), call. = FALSE)
  }
  invisible(df)
}

source(file.path(
  normalizePath(bmi_default_repo_root(), mustWork = FALSE),
  "scripts",
  "bmi_pheno_config.R"
))

if (!nzchar(Sys.getenv("BMI_IMAGING_ANALYSIS_DIR", unset = ""))) {
  if (nzchar(Sys.getenv("BMI_PHENO_VARIABLE", unset = ""))) {
    Sys.setenv(BMI_IMAGING_ANALYSIS_DIR = bmi_pheno_profile()$analysis_dir)
  }
}

bmi_ensure_run_tag_env()

bmi_project_paths <- function() {
  paths <- bmi_project_paths_core()
  if (!nzchar(Sys.getenv("BMI_PHENO_VARIABLE", unset = ""))) {
    stop("Set BMI_PHENO_VARIABLE (bmip or perc_median) before calling bmi_project_paths().", call. = FALSE)
  }
  bmi_ensure_run_tag_env()
  pheno <- bmi_pheno_profile()
  merged_override <- trimws(Sys.getenv("BMI_IMAGING_MERGED_DATAFILE", unset = ""))
  if (nzchar(merged_override)) {
    pheno$merged_datafile <- merged_override
  }
  ceph <- paths$ceph_project
  analysis_base <- bmi_resolve_analysis_base(ceph)
  final_datafiles <- bmi_resolve_profile_datafiles_dir(analysis_base, ceph)
  designmats_dir <- bmi_resolve_designmats_dir(analysis_base, ceph)
  results_root <- bmi_resolve_results_root(analysis_base, ceph)
  results_run <- Sys.getenv("BMI_IMAGING_RESULTS_RUN", unset = "")
  if (!nzchar(results_run)) {
    results_run <- paste0("results_", bmi_run_tag())
    Sys.setenv(BMI_IMAGING_RESULTS_RUN = results_run)
  }

  paths$pheno <- pheno
  paths$analysis_base <- analysis_base
  paths$final_analysis <- analysis_base
  paths$final_datafiles <- final_datafiles
  paths$profile_datafiles <- final_datafiles
  paths$designmats_dir <- designmats_dir
  paths$supplementary_designmats <- designmats_dir
  paths$results_root <- results_root
  paths$results_run <- results_run
  paths$results_dir <- file.path(results_root, results_run)
  paths$run_tag <- Sys.getenv("BMI_IMAGING_RUN_TAG", unset = bmi_run_tag())
  paths$explore_figures_dir <- file.path(paths$ceph_project, "shared", "bmip_metric_figures")
  paths
}
