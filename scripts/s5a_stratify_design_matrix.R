################################
# Stratify design matrix by sex
# Expects combined designmat from s5_make_design.R, e.g.:
#   designmat1_SbmiSageSpdsSexIncEducPCsScanSoft_nodemean.txt  (bmip)
#   designmat1_SpercSageSpdsSexIncEducPCsScanSoft_nodemean.txt  (perc_median)
################################

suppressPackageStartupMessages(library(dplyr))
# plyr (loaded by s5 / makeDesign) masks dplyr::n() — use dplyr:: explicitly below

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))

PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

designmat_dir <- PATHS$supplementary_designmats
bmi_ensure_dir(designmat_dir)
cat("Design matrix directory (created if needed):", designmat_dir, "\n")

cat("Analysis bundle:", PATHS$analysis_base, "\n")
if (!is.null(PHENO)) {
  cat("Pheno:", PHENO$pheno_var, "| expected tag:", PHENO$design_tag, "\n")
}

# Combined (non-stratified) design matrices: designmat<number>_<model>.txt
# Excludes files already stratified (designmat<num>m_ or designmat<num>f_)
all_files <- list.files(designmat_dir, pattern = "^designmat\\d+_.*\\.txt$", full.names = TRUE)
all_files <- all_files[!grepl("designmat\\d+[mf]_", basename(all_files))]

if (!length(all_files)) {
  stop(
    "No combined design matrices found in ", designmat_dir,
    "\nRun s5_make_design.R first.",
    call. = FALSE
  )
}

cat(sprintf("Found %d combined design matrix(ices) to stratify\n", length(all_files)))

for (input_file in all_files) {
  fname <- basename(input_file)

  num <- sub("^designmat(\\d+)_.*", "\\1", fname)
  model <- sub("^designmat\\d+_(.*)\\.txt$", "\\1", fname)

  model_no_sex <- gsub("Sex", "", model)
  output_males <- file.path(designmat_dir, paste0("designmat", num, "m_", model_no_sex, "_males.txt"))
  output_females <- file.path(designmat_dir, paste0("designmat", num, "f_", model_no_sex, "_females.txt"))

  cat(sprintf("\n###############################################################\n"))
  cat(sprintf("Design Matrix %s: %s\n", num, fname))
  cat(sprintf("###############################################################\n"))

  design_mat <- read.delim(input_file, sep = "\t")
  bmi_assert_positive_nrow(design_mat, paste("Design matrix", fname))
  bmi_assert_cols(design_mat, c("session_id", "ab_g_stc__cohort_sexF"), paste("Design matrix", fname))
  cat(sprintf("Loaded design matrix with %d rows and %d columns\n", nrow(design_mat), ncol(design_mat)))

  design_mat_males <- design_mat %>%
    dplyr::filter(ab_g_stc__cohort_sexF == 0) %>%
    dplyr::select(-ab_g_stc__cohort_sexF)

  design_mat_females <- design_mat %>%
    dplyr::filter(ab_g_stc__cohort_sexF == 1) %>%
    dplyr::select(-ab_g_stc__cohort_sexF)

  cat(sprintf("\nFemales: %d rows\n", nrow(design_mat_females)))
  cat(sprintf("Males: %d rows\n", nrow(design_mat_males)))

  # Save stratified matrices first (so reporting errors cannot block outputs)
  write.table(design_mat_females, file = output_females, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(design_mat_males, file = output_males, sep = "\t", row.names = FALSE, quote = FALSE)
  bmi_assert_file(output_females, "Stratified female design matrix")
  bmi_assert_file(output_males, "Stratified male design matrix")

  cat(sprintf("\nDesign matrices saved:\n"))
  cat(sprintf("  Females: %s\n", output_females))
  cat(sprintf("  Males: %s\n", output_males))

  file.remove(input_file)
  cat(sprintf("  Deleted original: %s\n", fname))

  cat("\n=== Observations by Session (optional report) ===\n")
  tryCatch({
    cat("\nCombined (both sexes):\n")
    session_counts_combined <- design_mat %>%
      dplyr::group_by(session_id) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(session_id)
    print(session_counts_combined)
    cat(sprintf("Total: %d\n", sum(session_counts_combined$n)))

    cat("\nFemales:\n")
    session_counts_females <- design_mat_females %>%
      dplyr::group_by(session_id) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(session_id)
    print(session_counts_females)
    cat(sprintf("Total: %d\n", sum(session_counts_females$n)))

    cat("\nMales:\n")
    session_counts_males <- design_mat_males %>%
      dplyr::group_by(session_id) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(session_id)
    print(session_counts_males)
    cat(sprintf("Total: %d\n", sum(session_counts_males$n)))
  }, error = function(e) {
    warning("Session count tables skipped: ", conditionMessage(e), call. = FALSE)
  })

  cat("\n================================\n")
}

cat("\nDone.\n")
