# Generate pheno / age / PDS vectors for basis-function modeling.
# Age/PDS vectors -> shared/datafiles; pheno vector -> <profile>/datafiles.

required_pkgs <- c("tidyverse", "dplyr")
lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
})

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

clean_and_write_vector <- function(df, var, out_prefix, out_dir) {
  clean <- df[!is.na(df[[var]]), c("participant_id", "session_id", var)]
  write.table(clean[[var]], file.path(out_dir, paste0(out_prefix, ".txt")),
              row.names = FALSE, col.names = FALSE)
}

base_name <- PATHS$base_datafile_imaging
base_file <- file.path(PATHS$shared_datafiles, base_name)
if (!file.exists(base_file)) {
  stop("Imaging base datafile not found: ", base_file, "\nRun s0/s1 after sourcing env/run_shared.csh.", call. = FALSE)
}

bmi_ensure_dir(PATHS$shared_datafiles)
bmi_ensure_dir(PATHS$profile_datafiles)
cat("Reading imaging base datafile from:", base_file, "\n")

df <- read.delim(base_file)

cat("Pheno variable:", PHENO$pheno_var, "(", PHENO$label, ")\n", sep = "")
cat("Shared vectors ->", PATHS$shared_datafiles, "\n")
cat("Pheno vector ->", PATHS$profile_datafiles, "\n")

clean_and_write_vector(df, "ab_g_dyn__visit_age", "age_vector_from_data", PATHS$shared_datafiles)
clean_and_write_vector(df, "pds_categ_y_p_all", "pds_vector_from_data", PATHS$shared_datafiles)
clean_and_write_vector(df, PHENO$pheno_var, PHENO$vector_prefix, PATHS$profile_datafiles)

v <- df[[PHENO$pheno_var]]
v <- v[!is.na(v)]
cat("\n", PHENO$pheno_var, " summary:\n", sep = "")
print(summary(v))
