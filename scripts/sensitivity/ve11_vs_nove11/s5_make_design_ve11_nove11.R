###############################################################
### VE11 vs noVE11 (6.0 matched): s(pheno)+s(age)+s(pds)+covariates
# Combined-sex (no s5a). Arm-specific scanner nuisance:
#   noVE11 -> ScanSoft (serial_software), modal ref in arm
#   VE11   -> scanner serial only, modal ref in arm
#
# Required env (see env/run_ve11_vs_nove11_arm.env):
#   BMI_PHENO_VARIABLE   bmip | perc_median
#   BMI_SCANNER_ARM      VE11 | noVE11
###############################################################

for (p in c("plyr", "dplyr")) {
  if (!eval(parse(text = paste0("require(", p, ")")))) {
    install.packages(p, repos = "https://cloud.r-project.org")
    lapply(p, library, character.only = TRUE)
  }
}

imaging_repo_root <- Sys.getenv(
  "BMI_IMAGING_REPO_ROOT",
  "/home/arigby/projects/ABCD_bmi_microstructure"
)
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

arm <- trimws(Sys.getenv("BMI_SCANNER_ARM", unset = ""))
if (!nzchar(arm) || !(arm %in% c("VE11", "noVE11"))) {
  stop("Set BMI_SCANNER_ARM to VE11 or noVE11.", call. = FALSE)
}

source(file.path(PATHS$cmig_tools_r, "makeDesign.R"))

outpath <- PATHS$supplementary_designmats
bmi_ensure_dir(outpath)

ndafile <- file.path(PATHS$final_datafiles, PHENO$merged_datafile)
bmi_assert_file(ndafile, "Merged matched datafile with pheno BFs")
nda <- read.delim(ndafile, stringsAsFactors = FALSE)
bmi_assert_positive_nrow(nda, "NDA table")
bmi_assert_cols(
  nda,
  c("participant_id", "session_id", "mr_y_qc__incl__dmri_indicator", "bmi_ve11_arm"),
  "NDA table"
)

analysis_root <- trimws(Sys.getenv("BMI_VE11_NOVE11_ANALYSIS_ROOT", unset = "ve11_vs_noVE11_6_0"))
visit_file <- file.path(
  PATHS$ceph_project, analysis_root, "matching",
  if (arm == "VE11") "matched_ve11_visits.txt" else "matched_nove11_visits.txt"
)
bmi_assert_file(visit_file, paste("Matched visit list for", arm))
visits <- read.delim(visit_file, stringsAsFactors = FALSE)
nda <- nda %>%
  inner_join(visits, by = c("participant_id", "session_id")) %>%
  filter(bmi_ve11_arm == arm)

cat(sprintf("Arm %s: %d matched visits in merged table\n", arm, nrow(nda)))

gen_pc_cols <- grep("^ab_g_stc__gen_pc__", names(nda), value = TRUE)
nda[gen_pc_cols] <- lapply(nda[gen_pc_cols], function(x) as.numeric(x))

idx_dmri_inc <- which(nda$mr_y_qc__incl__dmri_indicator == 1)
nda_dmri_inc <- nda[idx_dmri_inc, ]

scanner_col <- if (arm == "VE11") {
  "ab_g_dyn__design_mr__serial"
} else {
  "ab_g_dyn__design_mr__serial_software"
}
bmi_assert_cols(nda_dmri_inc, scanner_col, "Arm scanner covariate")

vars_interest <- c(
  "participant_id", "session_id", "ab_g_stc__design_id__fam", "ab_g_dyn__visit_age",
  paste0("ab_g_stc__gen_pc__", sprintf("%02d", 1:10)), "ab_g_stc__cohort_sex",
  "ab_g_dyn__cohort_edu__cgs", "ab_g_dyn__cohort_income__hhold__6lvl",
  scanner_col,
  names(nda)[grep("bf_", names(nda))]
)

nda_dmri_inc_completecases <- stats::na.omit(nda_dmri_inc[, vars_interest, drop = FALSE])
bmi_assert_positive_nrow(nda_dmri_inc_completecases, "NDA complete cases")

nda_dmri_inc_completecases$ab_g_stc__cohort_sex <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_stc__cohort_sex), ref = "M"
)
nda_dmri_inc_completecases$ab_g_dyn__cohort_edu__cgs <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_dyn__cohort_edu__cgs), ref = "5"
)
nda_dmri_inc_completecases$ab_g_dyn__cohort_income__hhold__6lvl <- stats::relevel(
  as.factor(nda_dmri_inc_completecases$ab_g_dyn__cohort_income__hhold__6lvl), ref = "5"
)

# Modal (or env) reference for arm-specific scanner covariate
ss_tab <- sort(table(nda_dmri_inc_completecases[[scanner_col]]), decreasing = TRUE)
cat(sprintf("%s level counts (%s):\n", scanner_col, arm))
print(ss_tab)

if (arm == "VE11") {
  ref_env <- trimws(Sys.getenv("BMI_SCANNER_SERIAL_REF", unset = ""))
  modal_flag <- tolower(trimws(Sys.getenv("BMI_SCANNER_SERIAL_MODAL_REF", unset = "1")))
  design_suffix <- "Scanner"
} else {
  ref_env <- trimws(Sys.getenv("BMI_SCANSOFT_REF", unset = ""))
  modal_flag <- tolower(trimws(Sys.getenv("BMI_SCANSOFT_MODAL_REF", unset = "1")))
  design_suffix <- "ScanSoft"
}
use_modal <- modal_flag %in% c("1", "true", "yes", "on")
if (nzchar(ref_env)) {
  scanner_ref <- ref_env
} else if (use_modal) {
  scanner_ref <- names(ss_tab)[1]
} else {
  stop("Set modal ref flag or explicit scanner reference for arm ", arm, call. = FALSE)
}
if (!(scanner_ref %in% names(ss_tab))) {
  stop("Scanner reference not present in arm sample: ", scanner_ref, call. = FALSE)
}
cat("Scanner reference level:", scanner_ref, "\n")
nda_dmri_inc_completecases[[scanner_col]] <- stats::relevel(
  as.factor(nda_dmri_inc_completecases[[scanner_col]]),
  ref = scanner_ref
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
  scanner_col
)
time <- c("ses-00A", "ses-02A", "ses-04A", "ses-06A")

fname <- sprintf(
  "designmat1_%sSageSpdsSexIncEducPCs%s_nodemean.txt",
  PHENO$design_tag, design_suffix
)
outfile <- file.path(outpath, fname)

cat("Pheno:", PHENO$pheno_var, "| Arm:", arm, "| Writing", outfile, "\n")

makeDesign(
  nda_dmri_inc_completecases, outfile, time,
  contvar = contvar, catvar = catvar,
  delta = NULL, interact = NULL, subjs = NULL,
  demean = FALSE, quadratic = NULL,
  familyID = "ab_g_stc__design_id__fam"
)

cat("Done.\n")
