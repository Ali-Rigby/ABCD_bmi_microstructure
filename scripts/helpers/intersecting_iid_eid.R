################################
# intersecting_iid_eid.R
# Merge FEMA iid/eid visit lists with profile merged datafile -> model_sample_<modality>.csv
################################

suppressPackageStartupMessages(library(dplyr))

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
if (!nzchar(Sys.getenv("BMI_PHENO_VARIABLE", unset = ""))) {
  stop("Set BMI_PHENO_VARIABLE (bmip or perc_median) before running.", call. = FALSE)
}

source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_model_sample_utils.R"))

PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

modalities <- trimws(strsplit(
  Sys.getenv("BMI_IMAGING_EXTRACT_MODALITIES", "rni"),
  ",", fixed = TRUE
)[[1]])
modalities <- modalities[nzchar(modalities)]

datafile <- file.path(PATHS$final_datafiles, PHENO$merged_datafile)
bmi_assert_file(datafile, "Merged datafile for model sample")

results_base_dir <- PATHS$results_dir
output_base_dir <- file.path(PATHS$analysis_base, "sample_tables", "voxelwise")
bmi_ensure_dir(output_base_dir)

models <- bmi_discover_designmat_models(results_base_dir)
if (!length(models)) {
  models <- bmi_designmat_model_dirs(PHENO$design_tag)
}

cat("\nintersecting_iid_eid\n")
cat("  pheno:          ", PHENO$pheno_var, "\n")
cat("  merged datafile:", datafile, "\n")
cat("  results_dir:    ", results_base_dir, "\n")
cat("  output_dir:     ", output_base_dir, "\n")
cat("  models:         ", paste(models, collapse = ", "), "\n\n")

df_all <- read.delim(datafile, check.names = FALSE)

for (modality in modalities) {
  for (model in models) {
    iid_eid_file <- file.path(results_base_dir, model, paste0("iid_eid_", modality, ".txt"))
    sample_out <- file.path(output_base_dir, model, paste0("model_sample_", modality, ".csv"))

    bmi_assert_file(iid_eid_file, paste("iid/eid list for", model))

    dir.create(dirname(sample_out), recursive = TRUE, showWarnings = FALSE)

    idlist <- bmi_read_iid_eid(iid_eid_file)
    df_combined <- merge(idlist, df_all, by = c("participant_id", "session_id"), all.x = TRUE)

    n_missing <- sum(is.na(df_combined$ab_g_stc__cohort_sex))
    if (n_missing > 0) {
      warning(
        n_missing, " visit(s) in ", model, " iid/eid not found in merged datafile.",
        call. = FALSE
      )
    }

    write.csv(df_combined, file = sample_out, row.names = FALSE)
    cat(sprintf(
      "Wrote %s (%s rows)\n",
      sample_out, format(nrow(df_combined), big.mark = ",")
    ))
  }
}

cat("\nintersecting_iid_eid: done.\n")
