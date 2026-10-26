################################
# tables/table1_demographics.R
#
# Title: Table 1. Demographic Data for All Observations, Stratified by Sex and Study Timepoint
# Structure: Female block stacked above Male block; columns = Overall + study timepoints
# Unit: imaging visits/observations
################################

required_pkgs <- c("tidyverse", "gtsummary", "gt")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(
      "Package '", p, "' is required. Example:\n",
      "  conda install -c conda-forge r-gtsummary r-gt r-tidyverse r-arrow",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
if (!nzchar(Sys.getenv("BMI_PHENO_VARIABLE", unset = ""))) {
  Sys.setenv(BMI_PHENO_VARIABLE = "bmip")
}
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_model_sample_utils.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_demos_data.R"))

PATHS <- bmi_project_paths()
PHENO <- PATHS$pheno
out_dir <- file.path(PATHS$ceph_project, "shared", "demographics")
bmi_ensure_dir(out_dir)

index_session <- tolower(Sys.getenv("BMI_IMAGING_TABLE1_INDEX_SESSION", unset = "baseline"))
skip_supplement <- tolower(Sys.getenv("BMI_IMAGING_TABLE1_SKIP_SUPPLEMENT", unset = "1")) %in%
  c("1", "true", "yes", "on")
legacy_participant <- tolower(Sys.getenv("BMI_IMAGING_TABLE1_LEGACY_PARTICIPANT", unset = "0")) %in%
  c("1", "true", "yes", "on")

sample_tables_dir <- Sys.getenv("BMI_IMAGING_SAMPLE_TABLES_DIR", unset = "")
loaded <- bmi_demos_load_model_sample(PATHS, PHENO, sample_tables_dir)
df <- loaded$df
sample_source <- loaded$sample_source

bmi_assert_cols(df, c(
  "participant_id", "session_id", "ab_g_stc__cohort_sex", "ab_g_dyn__visit_age",
  "ab_g_dyn__cohort_edu__cgs", "ab_g_dyn__cohort_income__hhold__6lvl",
  "pds_categ_y_p_all", "anthro_bmi_corr", "perc_median"
), "Model sample")

df <- bmi_demos_merge_from_parquet(
  df,
  c("ab_g_stc__cohort_race__nih", "ab_g_stc__cohort_ethn"),
  PATHS$abcd_tabulated
)

cat("\nTable 1 sample source:", sample_source, "\n")
cat("Pheno profile:", PHENO$pheno_var, "\n")

df_visits <- bmi_demos_add_study_wave(df)
label_maps <- bmi_demos_label_maps()

build_table_vars <- function(include_sex = FALSE, sex_var = "sex_group") {
  vars <- c("age_years")
  if (include_sex) {
    vars <- c(vars, sex_var)
  }
  vars <- c(vars, "bmi_kgm2")
  if ("bmip" %in% names(df_visits)) {
    vars <- c(vars, "bmip_score")
  }
  vars <- c(vars, "perc_from_median", "pubertal_stage")
  if ("bmip" %in% names(df_visits)) {
    vars <- c(vars, "cdc_weight_class")
  }
  c(
    vars,
    "race_nih",
    "ethnicity",
    "household_income",
    "parental_education"
  )
}

build_observation_footnote <- function(
  size_summary,
  include_bmip = TRUE,
  next_letter = "d",
  stacked_by_sex = FALSE
) {
  wave_text <- bmi_demos_format_wave_counts(size_summary$wave_counts)
  sex_wave_text <- bmi_demos_format_sex_wave_counts(size_summary$wave_counts_by_sex)
  participant_text <- paste0(
    "The sample comprised ",
    bmi_demos_format_count(size_summary$n_participants),
    " unique participants (",
    bmi_demos_format_count(size_summary$n_participants_female),
    " female; ",
    bmi_demos_format_count(size_summary$n_participants_male),
    " male)."
  )
  observation_total_text <- paste0(
    "A total of ",
    bmi_demos_format_count(size_summary$n_obs),
    " observations were included (",
    bmi_demos_format_count(size_summary$n_obs_female),
    " female; ",
    bmi_demos_format_count(size_summary$n_obs_male),
    " male)."
  )
  structure_text <- if (isTRUE(stacked_by_sex)) {
    paste(
      "Female and male results are presented in stacked blocks.",
      "Within each block, the Overall column summarizes all observations for that sex,",
      "and wave columns summarize observations within that sex and study wave.",
      "Column percentages are computed within each column of the corresponding sex block."
    )
  } else {
    paste(
      "Participants may contribute more than one observation; visit-level summaries therefore weight",
      "participants with repeated imaging more heavily than those with a single visit."
    )
  }
 
  observation_text <- paste0(
    observation_total_text,
    " ",
    participant_text,
    if (nzchar(wave_text)) {
      paste0(" Observations by study timepoint (all participants): ", wave_text, ".")
    } else {
      ""
    },
    " Characteristics are summarized across all eligible imaging visits in the voxelwise model sample.",
    " ",
    structure_text
  )

  bmi_demos_footnote_config(
    include_bmip = include_bmip,
    stacked_sex_wave_obs_n = isTRUE(stacked_by_sex),
    extra_entries = list(
      list(
        letter = next_letter,
        text = observation_text,
        attach = NULL
      )
    )
  )
}

################################
# Main Table 1: visit-level, stacked Female/Male blocks by study wave
################################

df_main <- df_visits %>%
  bmi_demos_mutate_characteristics(label_maps)

main_vars <- build_table_vars(include_sex = FALSE)
size_summary <- bmi_demos_sample_size_summary(df_main)

cat(sprintf(
  "\nMain Table 1 observations: %s overall (Female: %s, Male: %s)\n",
  bmi_demos_format_count(size_summary$n_obs),
  bmi_demos_format_count(size_summary$n_obs_female),
  bmi_demos_format_count(size_summary$n_obs_male)
))
cat(sprintf(
  "  Unique participants: %s (Female: %s, Male: %s)\n",
  bmi_demos_format_count(size_summary$n_participants),
  bmi_demos_format_count(size_summary$n_participants_female),
  bmi_demos_format_count(size_summary$n_participants_male)
))
if (!is.null(size_summary$wave_counts) && nrow(size_summary$wave_counts)) {
  cat("  By study wave:", bmi_demos_format_wave_counts(size_summary$wave_counts), "\n")
}

df_main_tbl <- df_main %>%
  dplyr::select(sex_group, study_wave, dplyr::all_of(main_vars))

tbl_main <- bmi_demos_build_stacked_sex_wave_table(
  df_tbl = df_main_tbl,
  table_vars = main_vars,
  header_count_descriptor = "observations"
)

main_footnotes <- build_observation_footnote(
  size_summary,
  include_bmip = "bmip_score" %in% main_vars,
  next_letter = if ("bmip_score" %in% main_vars) "e" else "d",
  stacked_by_sex = TRUE
)

bmi_demos_save_table(
  tbl_main,
  table_number = "1",
  table_title = "Demographic Data for All Observations, Stratified by Sex and Study Timepoint",
  footnote_config = main_footnotes,
  out_dir = out_dir,
  file_prefix = "table1_demos"
)

################################
# Optional legacy Table 1: participant-level by sex (baseline-indexed)
################################

if (isTRUE(legacy_participant)) {
  df_legacy <- bmi_demos_to_participant_level(df_visits, index_session = index_session) %>%
    bmi_demos_mutate_characteristics(label_maps) %>%
    dplyr::mutate(sex_display = sex_group)

  legacy_vars <- build_table_vars(include_sex = TRUE, sex_var = "sex_display")
  n_participants <- dplyr::n_distinct(df_legacy$participant_id)
  n_female <- sum(df_legacy$sex_group == "Female", na.rm = TRUE)
  n_male <- sum(df_legacy$sex_group == "Male", na.rm = TRUE)

  cat(sprintf(
    "\nLegacy participant Table 1: %s overall (Female: %s, Male: %s); index session = %s\n",
    bmi_demos_format_count(n_participants),
    bmi_demos_format_count(n_female),
    bmi_demos_format_count(n_male),
    index_session
  ))

  df_legacy_tbl <- df_legacy %>%
    dplyr::select(sex_group, dplyr::all_of(legacy_vars))

  tbl_legacy <- bmi_demos_build_table(
    df_tbl = df_legacy_tbl,
    by_var = "sex_group",
    table_vars = legacy_vars,
    header_overall_n = n_participants,
    header_level_ns = c(Female = n_female, Male = n_male),
    header_overall_label = "Overall",
    header_level_prefix = "N"
  )

  legacy_footnotes <- bmi_demos_footnote_config(
    include_bmip = "bmip_score" %in% legacy_vars,
    extra_entries = list(
      list(
        letter = if ("bmip_score" %in% legacy_vars) "e" else "d",
        text = paste0(
          "Participant-level characteristics were indexed at ",
          if (index_session == "baseline") {
            "baseline (session ses-00A), with earliest available visit used when baseline was absent."
          } else {
            "each participant's earliest eligible imaging visit."
          }
        ),
        attach = NULL
      )
    )
  )

  bmi_demos_save_table(
    tbl_legacy,
    table_number = "1",
    table_title = "Demographic Data for All Participants, Stratified by Sex",
    footnote_config = legacy_footnotes,
    out_dir = out_dir,
    file_prefix = "table1_demos_participant_level"
  )
}

################################
# Supplement Table S1: visit-level by study wave
################################

if (!skip_supplement) {
  df_supp <- df_visits %>%
    bmi_demos_mutate_characteristics(label_maps)

  supp_vars <- build_table_vars(include_sex = TRUE)
  supp_size <- bmi_demos_sample_size_summary(df_supp)
  wave_ns <- stats::setNames(supp_size$wave_counts$n, as.character(supp_size$wave_counts$study_wave))

  cat(sprintf(
    "\nSupplement Table S1 observations: %s overall (%s)\n",
    bmi_demos_format_count(supp_size$n_obs),
    bmi_demos_format_wave_counts(supp_size$wave_counts)
  ))

  df_supp_tbl <- df_supp %>%
    dplyr::select(study_wave, dplyr::all_of(supp_vars))

  tbl_supp <- bmi_demos_build_table(
    df_tbl = df_supp_tbl,
    by_var = "study_wave",
    table_vars = supp_vars,
    header_overall_n = supp_size$n_obs,
    header_level_ns = wave_ns,
    header_overall_label = "Overall",
    header_level_prefix = "N",
    header_count_descriptor = "observations"
  )

  supp_footnotes <- build_observation_footnote(
    supp_size,
    include_bmip = "bmip_score" %in% supp_vars,
    next_letter = if ("bmip_score" %in% supp_vars) "e" else "d"
  )

  bmi_demos_save_table(
    tbl_supp,
    table_number = "S1",
    table_title = "Demographic Data for All Observations, Stratified by Study Wave",
    footnote_config = supp_footnotes,
    out_dir = out_dir,
    file_prefix = "table_s1_demos_by_wave"
  )
}
