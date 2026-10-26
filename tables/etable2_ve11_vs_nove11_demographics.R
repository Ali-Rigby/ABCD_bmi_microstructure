################################
# etable2_ve11_vs_nove11_demographics.R
# Demographics Table 1 for VE11 vs noVE11 (6.0) analytic sample.
#
# Visit-level characteristics: Overall | VE11 | noVE11 | p-value
# Not sex-stratified; sex composition is a table row.
#
# Sample definition:
#   FEMA analytic visits (iid/eid) for bmip ∩ perc_median within each arm,
#   then pooled for Overall. Requires
#   scripts/sensitivity/ve11_vs_nove11/extract_iid_eid_ve11_vs_nove11.m
#   (auto-run if iid_eid_rni.txt files are missing).
#
# Between-arm tests (visit-level):
#   continuous  -> Wilcoxon rank-sum (Mann–Whitney U)
#   categorical -> chi-squared
# Supplemental CSV also reports two-sample Kolmogorov–Smirnov for continuous vars.
#
# Usage (gtsummary/gt live in the sci conda env on this cluster):
#   bash -lc 'source env/run_ve11_vs_noVE11_6_0_shared.env
#             export BMI_PHENO_VARIABLE=bmip
#             /home/arigby/programs/anaconda3/envs/sci/bin/Rscript tables/etable2_ve11_vs_nove11_demographics.R'
#
# Optional env:
#   BMI_VE11_NOVE11_ANALYSIS_ROOT       default ve11_vs_noVE11_6_0
#   BMI_PHENO_VARIABLE                  bmip|perc_median (BMI percentile / CDC class rows)
#   BMI_IMAGING_RESULTS_RUN_BMIP        default results_2026_09_02
#   BMI_IMAGING_RESULTS_RUN_PERC_MEDIAN default results_2026_10_02
#   BMI_IMAGING_REPO_ROOT               this repo
################################

required_pkgs <- c("tidyverse", "gtsummary", "gt", "rlang")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(
      "Package '", p, "' is required. On this cluster use the sci conda env:\n",
      "  /home/arigby/programs/anaconda3/envs/sci/bin/Rscript tables/etable2_ve11_vs_nove11_demographics.R\n",
      "Or: conda install -c conda-forge r-gtsummary r-gt r-tidyverse r-arrow",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

imaging_repo_root <- Sys.getenv(
  "BMI_IMAGING_REPO_ROOT",
  "/home/arigby/projects/ABCD_bmi_microstructure"
)
Sys.setenv(BMI_IMAGING_REPO_ROOT = imaging_repo_root)
if (!nzchar(Sys.getenv("BMI_PHENO_VARIABLE", unset = ""))) {
  Sys.setenv(BMI_PHENO_VARIABLE = "bmip")
}

analysis_root <- trimws(Sys.getenv("BMI_VE11_NOVE11_ANALYSIS_ROOT", unset = "ve11_vs_noVE11_6_0"))
Sys.setenv(BMI_IMAGING_SHARED_DATAFILES_DIR = file.path(analysis_root, "shared", "datafiles"))
# Profile paths only needed so config_project resolves; demos load the matched combined file.
Sys.setenv(BMI_IMAGING_ANALYSIS_DIR = file.path(analysis_root, Sys.getenv("BMI_PHENO_VARIABLE")))

source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_demos_data.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_model_sample_utils.R"))

PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno

matched_file <- file.path(PATHS$shared_datafiles, "datafile_matched_combined.txt")
bmi_assert_file(matched_file, "Matched combined datafile (run s1e matching first)")

out_dir <- file.path(PATHS$ceph_project, analysis_root, "demographics")
bmi_ensure_dir(out_dir)

run_bmip <- trimws(Sys.getenv("BMI_IMAGING_RESULTS_RUN_BMIP", unset = "results_2026_09_02"))
run_perc <- trimws(Sys.getenv("BMI_IMAGING_RESULTS_RUN_PERC_MEDIAN", unset = "results_2026_10_02"))
modality <- "rni"

ve11_iid_eid_path <- function(pheno, arm, results_run) {
  results_dir <- file.path(PATHS$ceph_project, analysis_root, pheno, arm, "results", results_run)
  models <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
  models <- models[grepl("/designmat1_[^/]+$", models)]
  files <- file.path(models, paste0("iid_eid_", modality, ".txt"))
  files <- files[file.exists(files)]
  if (!length(files)) {
    return(NA_character_)
  }
  if (length(files) > 1) {
    warning(
      "Multiple iid_eid files for ", pheno, "/", arm, "; using ", files[[1]],
      call. = FALSE
    )
  }
  files[[1]]
}

needed <- expand.grid(
  pheno = c("bmip", "perc_median"),
  arm = c("VE11", "noVE11"),
  stringsAsFactors = FALSE
)
needed$results_run <- ifelse(needed$pheno == "bmip", run_bmip, run_perc)
needed$path <- mapply(
  ve11_iid_eid_path, needed$pheno, needed$arm, needed$results_run,
  USE.NAMES = FALSE
)

if (anyNA(needed$path)) {
  cat("Missing iid_eid_rni.txt for one or more arms; running extract_iid_eid_ve11_vs_nove11.m\n")
  cmd <- sprintf(
    "cd %s && matlab -batch %s",
    shQuote(imaging_repo_root),
    shQuote("run('scripts/matlab_add_paths.m'); run('scripts/sensitivity/ve11_vs_nove11/extract_iid_eid_ve11_vs_nove11.m')")
  )
  status <- system(cmd)
  if (!identical(status, 0L)) {
    stop("Failed to extract iid/eid lists (matlab exit ", status, ").", call. = FALSE)
  }
  needed$path <- mapply(
    ve11_iid_eid_path, needed$pheno, needed$arm, needed$results_run,
    USE.NAMES = FALSE
  )
}
if (anyNA(needed$path)) {
  stop(
    "Still missing iid_eid_rni.txt after extract. Expected under:\n  ",
    paste(
      file.path(
        PATHS$ceph_project, analysis_root,
        needed$pheno[is.na(needed$path)],
        needed$arm[is.na(needed$path)],
        "results", needed$results_run[is.na(needed$path)], "designmat1_*/"
      ),
      collapse = "\n  "
    ),
    call. = FALSE
  )
}

# Force re-extract if lists still look like pre-intersect matched N (7931/arm)
pre_intersect_n <- 7931L
n_by_path <- vapply(needed$path, function(p) {
  max(0L, length(readLines(p, warn = FALSE)) - 1L)
}, integer(1))
if (all(n_by_path == pre_intersect_n)) {
  cat("iid_eid lists look pre-intersect (n=7931); re-running extract for post-FEMA QC lists\n")
  cmd <- sprintf(
    "cd %s && matlab -batch %s",
    shQuote(imaging_repo_root),
    shQuote("run('scripts/matlab_add_paths.m'); run('scripts/sensitivity/ve11_vs_nove11/extract_iid_eid_ve11_vs_nove11.m')")
  )
  status <- system(cmd)
  if (!identical(status, 0L)) {
    stop("Failed to re-extract post-intersect iid/eid (matlab exit ", status, ").", call. = FALSE)
  }
}

# Analytic visits = bmip ∩ perc_median within each arm, then pool for Overall
analytic_by_arm <- lapply(c("VE11", "noVE11"), function(arm) {
  p_bmip <- needed$path[needed$pheno == "bmip" & needed$arm == arm]
  p_perc <- needed$path[needed$pheno == "perc_median" & needed$arm == arm]
  keys_bmip <- bmi_visit_keys(bmi_read_iid_eid(p_bmip))
  keys_perc <- bmi_visit_keys(bmi_read_iid_eid(p_perc))
  keys <- dplyr::inner_join(keys_bmip, keys_perc, by = c("participant_id", "session_id"))
  cat(sprintf(
    "  %s analytic visits: bmip=%s, perc_median=%s, intersection=%s\n",
    arm,
    format(nrow(keys_bmip), big.mark = ","),
    format(nrow(keys_perc), big.mark = ","),
    format(nrow(keys), big.mark = ",")
  ))
  keys %>% dplyr::mutate(bmi_ve11_arm = arm)
})
analytic_keys <- dplyr::bind_rows(analytic_by_arm)

df <- read.delim(matched_file, check.names = FALSE, stringsAsFactors = FALSE)
bmi_assert_cols(df, c(
  "participant_id", "session_id", "bmi_ve11_arm", "ab_g_stc__cohort_sex",
  "ab_g_dyn__visit_age", "ab_g_dyn__cohort_edu__cgs",
  "ab_g_dyn__cohort_income__hhold__6lvl", "pds_categ_y_p_all",
  "anthro_bmi_corr", "perc_median"
), "Matched combined datafile")

df <- df %>%
  dplyr::mutate(
    participant_id = as.character(participant_id),
    session_id = as.character(session_id),
    bmi_ve11_arm = factor(bmi_ve11_arm, levels = c("VE11", "noVE11"))
  ) %>%
  dplyr::filter(!is.na(bmi_ve11_arm)) %>%
  dplyr::inner_join(analytic_keys, by = c("participant_id", "session_id", "bmi_ve11_arm")) %>%
  dplyr::mutate(bmi_ve11_arm = factor(bmi_ve11_arm, levels = c("VE11", "noVE11")))

n_matched <- length(readLines(matched_file)) - 1L
cat(sprintf(
  "Restricted matched combined (%s visits) -> analytic sample (%s visits)\n",
  format(n_matched, big.mark = ","),
  format(nrow(df), big.mark = ",")
))

df <- bmi_demos_merge_from_parquet(
  df,
  c("ab_g_stc__cohort_race__nih", "ab_g_stc__cohort_ethn"),
  PATHS$abcd_tabulated
)

df_visits <- bmi_demos_add_study_wave(df)
label_maps <- bmi_demos_label_maps()

df_main <- df_visits %>%
  bmi_demos_mutate_characteristics(label_maps) %>%
  dplyr::mutate(
    scanner_arm = factor(as.character(bmi_ve11_arm), levels = c("VE11", "noVE11"))
  )

build_table_vars <- function() {
  vars <- c(
    "sex_group",
    "age_years",
    "bmi_kgm2"
  )
  if ("bmip" %in% names(df_main) && any(is.finite(df_main$bmip_score))) {
    vars <- c(vars, "bmip_score")
  }
  vars <- c(vars, "perc_from_median", "pubertal_stage")
  if ("cdc_weight_class" %in% names(df_main) && !all(is.na(df_main$cdc_weight_class))) {
    vars <- c(vars, "cdc_weight_class")
  }
  c(
    vars,
    "study_wave",
    "race_nih",
    "ethnicity",
    "household_income",
    "parental_education"
  )
}

# study_wave label is not in the shared map yet
local_var_labels <- c(
  bmi_demos_table_var_labels(),
  study_wave = "Study wave, No. (%)"
)

main_vars <- build_table_vars()
size_summary <- bmi_demos_sample_size_summary(df_main)
arm_ns <- df_main %>%
  dplyr::count(scanner_arm, name = "n") %>%
  dplyr::arrange(scanner_arm)
arm_ns_vec <- stats::setNames(arm_ns$n, as.character(arm_ns$scanner_arm))

n_part_by_arm <- df_main %>%
  dplyr::group_by(scanner_arm) %>%
  dplyr::summarise(
    n_participants = dplyr::n_distinct(participant_id),
    n_female = dplyr::n_distinct(participant_id[sex_group == "Female"]),
    n_male = dplyr::n_distinct(participant_id[sex_group == "Male"]),
    .groups = "drop"
  )

cat("\nVE11 vs noVE11 Table 1\n")
cat("  Sample source: post-FEMA_intersect_design visits (bmip ∩ perc_median) within each arm\n")
cat("  Demographics covariates from: ", matched_file, "\n", sep = "")
cat("  Results runs: bmip=", run_bmip, "; perc_median=", run_perc, "\n", sep = "")
cat("  Pheno profile (for BMI percentile / CDC class rows): ", PHENO$pheno_var, "\n", sep = "")
cat(sprintf(
  "  Observations: %s overall (VE11: %s, noVE11: %s)\n",
  bmi_demos_format_count(size_summary$n_obs),
  bmi_demos_format_count(arm_ns_vec[["VE11"]]),
  bmi_demos_format_count(arm_ns_vec[["noVE11"]])
))
cat(sprintf(
  "  Unique participants: %s (Female: %s, Male: %s)\n",
  bmi_demos_format_count(size_summary$n_participants),
  bmi_demos_format_count(size_summary$n_participants_female),
  bmi_demos_format_count(size_summary$n_participants_male)
))
print(as.data.frame(n_part_by_arm), row.names = FALSE)

df_tbl <- df_main %>%
  dplyr::select(scanner_arm, dplyr::all_of(main_vars))

labels <- as.list(local_var_labels[names(local_var_labels) %in% main_vars])
specs <- bmi_demos_jama_table_specs(main_vars)

tbl_main <- df_tbl %>%
  gtsummary::tbl_summary(
    by = scanner_arm,
    include = dplyr::all_of(main_vars),
    label = labels,
    type = specs$var_types,
    statistic = specs$statistic,
    digits = specs$digits,
    missing = "no",
    percent = "column"
  ) %>%
  gtsummary::add_overall(last = FALSE, col_label = "**Overall**") %>%
  gtsummary::add_p(
    test = list(
      gtsummary::all_continuous() ~ "wilcox.test",
      gtsummary::all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = gtsummary::label_style_pvalue(digits = 3)
  ) %>%
  gtsummary::bold_labels()

stat_cols <- grep("^stat_", names(tbl_main$table_body), value = TRUE)
tbl_main$table_body <- tbl_main$table_body %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(stat_cols),
      ~ bmi_demos_jama_fmt_iqr(.x)
    )
  )

format_header_count <- function(label, n) {
  paste0("**", label, "**\n(n = ", bmi_demos_format_count(n), " observations)")
}

header_args <- list(
  label = "**Characteristic**",
  stat_0 = format_header_count("Overall", size_summary$n_obs)
)
level_names <- names(arm_ns_vec)
for (i in seq_along(arm_ns_vec)) {
  header_args[[paste0("stat_", i)]] <- format_header_count(level_names[[i]], arm_ns_vec[[i]])
}
tbl_main <- rlang::exec(gtsummary::modify_header, tbl_main, !!!header_args)
tbl_main <- gtsummary::modify_header(tbl_main, p.value ~ "**p-value**")

arm_part_text <- paste(
  vapply(
    seq_len(nrow(n_part_by_arm)),
    function(i) {
      paste0(
        as.character(n_part_by_arm$scanner_arm[[i]]), ": ",
        bmi_demos_format_count(n_part_by_arm$n_participants[[i]]),
        " participants (",
        bmi_demos_format_count(n_part_by_arm$n_female[[i]]), " female; ",
        bmi_demos_format_count(n_part_by_arm$n_male[[i]]), " male)"
      )
    },
    character(1)
  ),
  collapse = "; "
)

footnote_cfg <- bmi_demos_jama_footnote_config(
  include_bmip = "bmip_score" %in% main_vars,
  extra_entries = list(
    list(
      letter = if ("bmip_score" %in% main_vars) "e" else "d",
      text = paste0(
        "Visit-level characteristics for the VE11 vs noVE11 voxelwise analytic sample ",
        "(visits retained after FEMA imaging QC / design–imaging intersect in both bmip and ",
        "perc_median models within each scanner arm, following 1:1 matching on ",
        "percent-from-median BMI within ±10 and age within ±0.25 years). ",
        "A total of ", bmi_demos_format_count(size_summary$n_obs), " observations were included (",
        bmi_demos_format_count(size_summary$n_obs_female), " female; ",
        bmi_demos_format_count(size_summary$n_obs_male), " male). ",
        "Unique participants overall: ", bmi_demos_format_count(size_summary$n_participants),
        " (", bmi_demos_format_count(size_summary$n_participants_female), " female; ",
        bmi_demos_format_count(size_summary$n_participants_male), " male). ",
        "By scanner arm — ", arm_part_text, ". ",
        if (nzchar(bmi_demos_format_wave_counts(size_summary$wave_counts))) {
          paste0(
            "Observations by study wave: ",
            bmi_demos_format_wave_counts(size_summary$wave_counts), ". "
          )
        } else {
          ""
        },
        "Sex composition is summarized as No. (%) of visits; the table is not sex-stratified. ",
        "Participants may contribute more than one observation."
      ),
      attach = NULL
    ),
    list(
      letter = if ("bmip_score" %in% main_vars) "f" else "e",
      text = paste(
        "P values compare VE11 vs noVE11 at the visit level using Wilcoxon rank-sum tests",
        "(Mann–Whitney U) for continuous variables and chi-squared tests for categorical variables.",
        "Because matching used percent-from-median BMI and age, those variables are expected to be similar by design.",
        "A supplemental CSV also reports two-sample Kolmogorov–Smirnov tests for continuous variables."
      ),
      attach = NULL
    )
  )
)

bmi_demos_jama_save_table(
  tbl_main,
  table_number = "1",
  table_title = "Demographic Data for Analytic Observations, Stratified by Scanner Software Arm (VE11 vs noVE11)",
  footnote_config = footnote_cfg,
  out_dir = out_dir,
  file_prefix = "table1_ve11_vs_noVE11_6_0"
)

################################
# Supplemental continuous arm comparisons: Wilcoxon + KS
################################

cont_vars <- intersect(
  c("age_years", "bmi_kgm2", "bmip_score", "perc_from_median", "pubertal_stage"),
  main_vars
)

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.4f", p)
}

cmp_rows <- lapply(cont_vars, function(v) {
  x_a <- df_main[[v]][df_main$scanner_arm == "VE11"]
  x_b <- df_main[[v]][df_main$scanner_arm == "noVE11"]
  x_a <- x_a[is.finite(x_a)]
  x_b <- x_b[is.finite(x_b)]

  wt <- tryCatch(
    stats::wilcox.test(x_a, x_b, exact = FALSE),
    error = function(e) NULL
  )
  ks <- tryCatch(
    suppressWarnings(stats::ks.test(x_a, x_b, exact = FALSE)),
    error = function(e) NULL
  )

  data.frame(
    variable = v,
    label = unname(local_var_labels[[v]]),
    n_VE11 = length(x_a),
    n_noVE11 = length(x_b),
    mean_VE11 = if (length(x_a)) mean(x_a) else NA_real_,
    sd_VE11 = if (length(x_a)) stats::sd(x_a) else NA_real_,
    mean_noVE11 = if (length(x_b)) mean(x_b) else NA_real_,
    sd_noVE11 = if (length(x_b)) stats::sd(x_b) else NA_real_,
    wilcoxon_W = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
    wilcoxon_p = if (!is.null(wt)) wt$p.value else NA_real_,
    wilcoxon_p_fmt = if (!is.null(wt)) fmt_p(wt$p.value) else NA_character_,
    ks_D = if (!is.null(ks)) unname(ks$statistic) else NA_real_,
    ks_p = if (!is.null(ks)) ks$p.value else NA_real_,
    ks_p_fmt = if (!is.null(ks)) fmt_p(ks$p.value) else NA_character_,
    stringsAsFactors = FALSE
  )
})

cmp_tbl <- dplyr::bind_rows(cmp_rows)
cmp_path <- file.path(out_dir, "table1_ve11_vs_noVE11_6_0_continuous_tests.csv")
readr::write_csv(cmp_tbl, cmp_path)

cat("\nContinuous VE11 vs noVE11 tests (Wilcoxon + KS):\n")
print(cmp_tbl %>% dplyr::select(label, wilcoxon_p_fmt, ks_D, ks_p_fmt), row.names = FALSE)
cat("\nWrote:\n  ", cmp_path, "\n", sep = "")
cat("Main Table 1 outputs under:\n  ", out_dir, "\n", sep = "")
