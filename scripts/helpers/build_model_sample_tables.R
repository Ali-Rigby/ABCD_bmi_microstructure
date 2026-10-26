################################
# build_model_sample_tables.R
# 1) extract_iid_eid.m for bmip and perc_median (unless BMI_SKIP_FEMA_EXTRACT=1)
# 2) intersecting_iid_eid.R for each profile
# 3) sanity-check that visit sets match across phenotypes
################################

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_model_sample_utils.R"))

phenos <- c("bmip", "perc_median")
skip_extract <- identical(tolower(Sys.getenv("BMI_SKIP_FEMA_EXTRACT", "0")), "1")

run_matlab_extract <- function(pheno) {
  profile <- bmi_pheno_profile(pheno)
  env <- c(
    sprintf("BMI_PHENO_VARIABLE=%s", pheno),
    sprintf("BMI_IMAGING_ANALYSIS_DIR=%s", profile$analysis_dir),
    sprintf("BMI_IMAGING_REPO_ROOT=%s", imaging_repo_root)
  )
  if (nzchar(Sys.getenv("BMI_IMAGING_CEPH_ROOT", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_CEPH_ROOT=%s", Sys.getenv("BMI_IMAGING_CEPH_ROOT")))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_RESULTS_RUN", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_RESULTS_RUN=%s", Sys.getenv("BMI_IMAGING_RESULTS_RUN")))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_RUN_TAG", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_RUN_TAG=%s", Sys.getenv("BMI_IMAGING_RUN_TAG")))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_AUTO_RUN_DATE", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_AUTO_RUN_DATE=%s", Sys.getenv("BMI_IMAGING_AUTO_RUN_DATE")))
  }

  matlab_cmd <- Sys.getenv("BMI_IMAGING_MATLAB_CMD", "matlab")
  batch <- sprintf(
    "addpath('%s'); addpath('%s'); run('%s'); run('%s');",
    file.path(imaging_repo_root, "scripts"),
    file.path(imaging_repo_root, "utilities", "matlab_utils"),
    file.path(imaging_repo_root, "scripts", "matlab_add_paths.m"),
    file.path(imaging_repo_root, "scripts", "helpers", "extract_iid_eid.m")
  )
  cat("\n=== MATLAB extract_iid_eid (", pheno, ") ===\n", sep = "")
  status <- system2(
    matlab_cmd,
    c("-batch", batch),
    env = c(Sys.getenv(), env)
  )
  if (status != 0) {
    stop("MATLAB extract_iid_eid failed for ", pheno, " (exit ", status, ")", call. = FALSE)
  }
}

run_intersect <- function(pheno) {
  profile <- bmi_pheno_profile(pheno)
  env <- c(
    sprintf("BMI_PHENO_VARIABLE=%s", pheno),
    sprintf("BMI_IMAGING_ANALYSIS_DIR=%s", profile$analysis_dir),
    sprintf("BMI_IMAGING_REPO_ROOT=%s", imaging_repo_root)
  )
  if (nzchar(Sys.getenv("BMI_IMAGING_CEPH_ROOT", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_CEPH_ROOT=%s", Sys.getenv("BMI_IMAGING_CEPH_ROOT")))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_RESULTS_RUN", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_RESULTS_RUN=%s", Sys.getenv("BMI_IMAGING_RESULTS_RUN")))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_RUN_TAG", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_RUN_TAG=%s", Sys.getenv("BMI_IMAGING_RUN_TAG")))
  }
  if (nzchar(Sys.getenv("BMI_IMAGING_AUTO_RUN_DATE", unset = ""))) {
    env <- c(env, sprintf("BMI_IMAGING_AUTO_RUN_DATE=%s", Sys.getenv("BMI_IMAGING_AUTO_RUN_DATE")))
  }

  cat("\n=== intersecting_iid_eid (", pheno, ") ===\n", sep = "")
  status <- system2(
    "Rscript",
    file.path(imaging_repo_root, "scripts", "helpers", "intersecting_iid_eid.R"),
    env = c(Sys.getenv(), env)
  )
  if (status != 0) {
    stop("intersecting_iid_eid failed for ", pheno, " (exit ", status, ")", call. = FALSE)
  }
}

for (pheno in phenos) {
  if (!skip_extract) {
    run_matlab_extract(pheno)
  } else {
    cat("\n=== Skipping MATLAB extract (BMI_SKIP_FEMA_EXTRACT=1) for ", pheno, " ===\n", sep = "")
  }
  run_intersect(pheno)
}

cat("\n=== Cross-pheno sanity check (bmip vs perc_median) ===\n")
visit_sets <- list()
for (pheno in phenos) {
  PATHS <- bmi_paths_for_pheno(pheno, imaging_repo_root)
  sample_dir <- file.path(PATHS$analysis_base, "sample_tables", "voxelwise")
  visit_sets[[pheno]] <- bmi_load_model_sample_visits(sample_dir, PATHS$pheno$design_tag, "rni")
  cat(sprintf("  %s: %s visits\n", pheno, format(nrow(visit_sets[[pheno]]), big.mark = ",")))
}

cmp <- bmi_compare_visit_sets(visit_sets$bmip, visit_sets$perc_median, "bmip", "perc_median")
bmi_print_visit_set_comparison(cmp, "bmip", "perc_median")

summary_path <- file.path(
  Sys.getenv("BMI_IMAGING_CEPH_ROOT", bmi_default_ceph_root()),
  "shared", "demographics", "model_sample_pheno_comparison.txt"
)
bmi_ensure_dir(dirname(summary_path))
writeLines(c(
  sprintf("bmip visits: %d", cmp$n_a),
  sprintf("perc_median visits: %d", cmp$n_b),
  sprintf("intersection: %d", cmp$n_intersection),
  sprintf("only bmip: %d", cmp$n_only_a),
  sprintf("only perc_median: %d", cmp$n_only_b),
  sprintf("match: %s", cmp$match)
), summary_path)
cat("\nWrote comparison summary:", summary_path, "\n")

if (!cmp$match) {
  stop("bmip and perc_median FEMA analytic samples do not match.", call. = FALSE)
}

cat("\nbuild_model_sample_tables: all profiles match.\n")
