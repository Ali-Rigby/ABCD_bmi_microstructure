###############################################################
### Design matrix 1: s(pheno) + s(age) + s(pds) + standard covariates
# Set BMI_PHENO_VARIABLE (bmip | perc_median) before running.
###############################################################

for (p in c("plyr", "dplyr")) {
  if (!eval(parse(text = paste0("require(", p, ")")))) {
    install.packages(p, repos = "https://cloud.r-project.org")
    lapply(p, library, character.only = TRUE)
  }
}

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

source(file.path(PATHS$cmig_tools_r, "makeDesign.R"))

outpath <- PATHS$supplementary_designmats
bmi_ensure_dir(outpath)

ndafile <- file.path(PATHS$final_datafiles, PHENO$merged_datafile)
bmi_assert_file(ndafile, "Merged datafile with pheno BFs")
nda <- read.delim(ndafile)
bmi_assert_positive_nrow(nda, "NDA table")
bmi_assert_cols(nda, c("participant_id", "session_id", "mr_y_qc__incl__dmri_indicator"), "NDA table")

gen_pc_cols <- grep("^ab_g_stc__gen_pc__", names(nda), value = TRUE)
nda[gen_pc_cols] <- lapply(nda[gen_pc_cols], function(x) as.numeric(x))

idx_dmri_inc <- which(nda$mr_y_qc__incl__dmri_indicator == 1)
nda_dmri_inc <- nda[idx_dmri_inc, ]

vars_interest <- c(
  "participant_id", "session_id", "ab_g_stc__design_id__fam", "ab_g_dyn__visit_age",
  paste0("ab_g_stc__gen_pc__", sprintf("%02d", 1:10)), "ab_g_stc__cohort_sex",
  "ab_g_dyn__cohort_edu__cgs", "ab_g_dyn__cohort_income__hhold__6lvl",
  "ab_g_dyn__design_mr__serial_software",
  names(nda)[grep("bf_", names(nda))]
)

nda_dmri_inc_completecases <- stats::na.omit(nda_dmri_inc[, vars_interest, drop = FALSE])
bmi_assert_positive_nrow(nda_dmri_inc_completecases, "NDA complete cases")

nda_dmri_inc_completecases$ab_g_stc__cohort_sex <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_stc__cohort_sex), ref = "M"
)
nda_dmri_inc_completecases$ab_g_dyn__design_mr__serial_software <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_dyn__design_mr__serial_software),
  ref = "HASH3935c89e_syngo MR E11 VE11C_LATEST_20160120"
)
nda_dmri_inc_completecases$ab_g_dyn__cohort_edu__cgs <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_dyn__cohort_edu__cgs), ref = "5"
)
nda_dmri_inc_completecases$ab_g_dyn__cohort_income__hhold__6lvl <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_dyn__cohort_income__hhold__6lvl), ref = "5"
)

bf_pat <- paste0("^", PHENO$bf_grep)
contvar <- c(
  names(nda)[grep(bf_pat, names(nda))],
  names(nda)[grep("bf_age", names(nda))],
  names(nda)[grep("bf_pds", names(nda))],
  paste0("ab_g_stc__gen_pc__", sprintf("%02d", 1:10))
)
catvar <- c(
  "ab_g_stc__cohort_sex",
  "ab_g_dyn__cohort_edu__cgs",
  "ab_g_dyn__cohort_income__hhold__6lvl",
  "ab_g_dyn__design_mr__serial_software"
)
time <- c("ses-00A", "ses-02A", "ses-04A", "ses-06A")

fname <- sprintf("designmat1_%sSageSpdsSexIncEducPCsScanSoft_nodemean.txt", PHENO$design_tag)
outfile <- file.path(outpath, fname)

cat("Pheno:", PHENO$pheno_var, "| Writing", outfile, "\n")

makeDesign(
  nda_dmri_inc_completecases, outfile, time,
  contvar = contvar, catvar = catvar,
  delta = NULL, interact = NULL, subjs = NULL,
  demean = FALSE, quadratic = NULL,
  familyID = "ab_g_stc__design_id__fam"
)

cat("Done.\n")
