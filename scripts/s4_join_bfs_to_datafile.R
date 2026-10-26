################################
# Merge pheno + age + PDS basis functions into profile datafile
# Base table from shared/datafiles; age/PDS BFs from shared; pheno BFs from profile.
################################

for (p in c("tidyverse", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

shared_dir <- PATHS$shared_datafiles
profile_dir <- PATHS$profile_datafiles
bmi_ensure_dir(profile_dir)

base_name <- PATHS$base_datafile_imaging
base_file <- file.path(shared_dir, base_name)
bmi_assert_file(base_file, "Imaging base datafile (run s1 first)")
cat("Reading imaging base datafile from:", base_file, "\n")

pheno_bf_file <- file.path(profile_dir, PHENO$bf_table_interp)
age_bf_file <- file.path(shared_dir, "table_bf_age_interp.txt")
pds_bf_file <- file.path(shared_dir, "table_bf_pds_interp.txt")
output_file <- file.path(profile_dir, PHENO$merged_datafile)

merge_basis_functions <- function(df_clean, basis_file, prefix, drop_col = NULL) {
  bfs <- read.delim(basis_file, header = FALSE)
  colnames(bfs) <- paste0(prefix, seq_len(ncol(bfs)))
  df_merged <- cbind(df_clean, bfs)
  if (!is.null(drop_col)) {
    df_merged <- dplyr::select(df_merged, -dplyr::all_of(drop_col))
  }
  df_merged
}

df <- read.delim(base_file)
final_df <- df

pv <- PHENO$pheno_var
df_clean_pheno <- df[!is.na(df[[pv]]), c("participant_id", "session_id", pv)]
df_basis_pheno <- merge_basis_functions(df_clean_pheno, pheno_bf_file, PHENO$bf_col_prefix, pv)
final_df <- merge(final_df, df_basis_pheno, by = c("participant_id", "session_id"), all.x = TRUE)

df_clean_age <- df[!is.na(df$ab_g_dyn__visit_age), c("participant_id", "session_id", "ab_g_dyn__visit_age")]
df_basis_age <- merge_basis_functions(df_clean_age, age_bf_file, "bf_age_", "ab_g_dyn__visit_age")
final_df <- merge(final_df, df_basis_age, by = c("participant_id", "session_id"), all.x = TRUE)

df_clean_pds <- df[!is.na(df$pds_categ_y_p_all), c("participant_id", "session_id", "pds_categ_y_p_all")]
df_basis_pds <- merge_basis_functions(df_clean_pds, pds_bf_file, "bf_pds_", "pds_categ_y_p_all")
final_df <- merge(final_df, df_basis_pds, by = c("participant_id", "session_id"), all.x = TRUE)

write.table(final_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote", output_file, "\n")
