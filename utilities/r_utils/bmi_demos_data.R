# Shared helpers for demographics table.

bmi_demos_load_model_sample <- function(PATHS, PHENO, sample_tables_dir = "") {
  if (!nzchar(sample_tables_dir)) {
    sample_tables_dir <- file.path(PATHS$analysis_base, "sample_tables", "voxelwise")
  }
  sample_tables_dir <- normalizePath(sample_tables_dir, mustWork = FALSE)

  design_female <- sprintf(
    "designmat1f_%sSageSpdsIncEducPCsScanSoft_nodemean_females",
    PHENO$design_tag
  )
  design_male <- sprintf(
    "designmat1m_%sSageSpdsIncEducPCsScanSoft_nodemean_males",
    PHENO$design_tag
  )

  female_file <- file.path(sample_tables_dir, design_female, "model_sample_rni.csv")
  male_file <- file.path(sample_tables_dir, design_male, "model_sample_rni.csv")

  load_pair <- function(path_female, path_male) {
    if (!file.exists(path_female) || !file.exists(path_male)) {
      return(NULL)
    }
    df_female <- readr::read_csv(
      path_female,
      show_col_types = FALSE,
      col_types = readr::cols(ab_g_stc__cohort_sex = readr::col_character())
    )
    df_male <- readr::read_csv(
      path_male,
      show_col_types = FALSE,
      col_types = readr::cols(ab_g_stc__cohort_sex = readr::col_character())
    )
    dplyr::bind_rows(df_female, df_male)
  }

  df <- load_pair(female_file, male_file)
  sample_source <- "profile sample_tables"

  if (is.null(df)) {
    merged_datafile <- file.path(PATHS$final_datafiles, PHENO$merged_datafile)
    bmi_assert_file(merged_datafile, "Merged datafile for model sample rebuild")
    df_all <- read.delim(merged_datafile, check.names = FALSE)

    results_dir <- PATHS$results_dir
    iid_female <- file.path(results_dir, design_female, "iid_eid_rni.txt")
    iid_male <- file.path(results_dir, design_male, "iid_eid_rni.txt")
    if (!file.exists(iid_female) || !file.exists(iid_male)) {
      stop(
        "Could not find model_sample_rni.csv or iid_eid_rni.txt files.\n",
        "  Tried female sample: ", female_file, "\n",
        "  Tried male sample:   ", male_file, "\n",
        "  Run scripts/helpers/intersecting_iid_eid.R (or build_model_sample_tables.R) after voxelwise results exist.",
        call. = FALSE
      )
    }

    read_iid_eid <- function(path) {
      idlist <- read.delim(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
      stopifnot(ncol(idlist) >= 2)
      names(idlist)[1:2] <- c("participant_id", "session_id")
      idlist[, c("participant_id", "session_id"), drop = FALSE]
    }

    idlist <- dplyr::bind_rows(read_iid_eid(iid_female), read_iid_eid(iid_male))
    df <- merge(idlist, df_all, by = c("participant_id", "session_id"), all.x = TRUE)
    sample_source <- "rebuilt from iid_eid + merged datafile"
  }

  if (!"perc_median" %in% names(df)) {
    merged_datafile <- file.path(PATHS$final_datafiles, PHENO$merged_datafile)
    if (!file.exists(merged_datafile)) {
      merged_datafile <- file.path(PATHS$shared_datafiles, PATHS$base_datafile_imaging)
    }
    if (file.exists(merged_datafile)) {
      pheno_cols <- read.delim(merged_datafile, check.names = FALSE, nrows = 1)
      join_cols <- intersect(
        c("participant_id", "session_id", "perc_median"),
        names(pheno_cols)
      )
      if (length(join_cols) == 3) {
        pheno_df <- read.delim(merged_datafile, check.names = FALSE)[, join_cols, drop = FALSE]
        pheno_df$participant_id <- as.character(pheno_df$participant_id)
        pheno_df$session_id <- as.character(pheno_df$session_id)
        pheno_df$perc_median <- as.numeric(pheno_df$perc_median)
        df <- dplyr::left_join(df, pheno_df, by = c("participant_id", "session_id"))
      }
    }
  }

  list(
    df = df,
    sample_source = sample_source,
    sample_tables_dir = sample_tables_dir,
    design_female = design_female,
    design_male = design_male
  )
}

bmi_demos_merge_from_parquet <- function(df, cols, abcd_tabulated) {
  need <- setdiff(cols, names(df))
  if (!length(need)) {
    return(df)
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(
      "Package 'arrow' is required to merge columns from ABCD parquet: ",
      paste(need, collapse = ", "),
      call. = FALSE
    )
  }

  stc_cols <- intersect(
    need,
    c(
      "ab_g_stc__cohort_race__nih",
      "ab_g_stc__cohort_ethn",
      "ab_g_stc__design_id__group"
    )
  )
  dyn_cols <- setdiff(need, stc_cols)

  if (length(stc_cols)) {
    stc_path <- file.path(abcd_tabulated, "ab_g_stc.parquet")
    bmi_assert_file(stc_path, "ab_g_stc.parquet")
    stc <- arrow::read_parquet(
      stc_path,
      col_select = c("participant_id", tidyselect::all_of(stc_cols))
    )
    df <- dplyr::left_join(df, stc, by = "participant_id")
  }

  if (length(dyn_cols)) {
    dyn_path <- file.path(abcd_tabulated, "ab_g_dyn.parquet")
    bmi_assert_file(dyn_path, "ab_g_dyn.parquet")
    dyn <- arrow::read_parquet(
      dyn_path,
      col_select = c("participant_id", "session_id", tidyselect::all_of(dyn_cols))
    )
    df <- dplyr::left_join(df, dyn, by = c("participant_id", "session_id"))
  }

  df
}

bmi_demos_session_levels <- function() {
  c("ses-00A", "ses-02A", "ses-04A", "ses-06A")
}

bmi_demos_session_labels <- function() {
  c("Baseline", "2-Year FU", "4-Year FU", "6-Year FU")
}

bmi_demos_add_study_wave <- function(df) {
  session_levels <- bmi_demos_session_levels()
  session_labels <- bmi_demos_session_labels()
  df %>%
    dplyr::mutate(
      session_id = as.character(session_id),
      study_wave = factor(session_id, levels = session_levels, labels = session_labels)
    ) %>%
    dplyr::filter(!is.na(study_wave))
}

bmi_demos_label_maps <- function() {
  list(
    edu = c(
      "1" = "Up to high school (no diploma)",
      "2" = "High school diploma/GED",
      "3" = "Some college",
      "4" = "Bachelor's degree",
      "5" = "Graduate or professional degree"
    ),
    income = c(
      "1" = "<$25 000",
      "2" = "$25 000-$49 999",
      "3" = "$50 000-$74 999",
      "4" = "$75 000-$99 999",
      "5" = "$100 000-$199 999",
      "6" = "≥$200 000",
      "777" = "Declined to answer",
      "999" = "Don't know"
    ),
    race = c(
      "2" = "White",
      "3" = "Black or African American",
      "4" = "Asian",
      "5" = "American Indian or Alaska Native",
      "6" = "Native Hawaiian or Other Pacific Islander",
      "8" = "More than 1 race",
      "13" = "Other or unknown"
    ),
    ethn = c(
      "1" = "Hispanic",
      "2" = "Non-Hispanic"
    ),
    visit_count = c("1" = "1", "2" = "2", "3" = "3", "4" = "4")
  )
}

bmi_demos_cdc_weight_class_levels <- function() {
  c("Underweight", "Healthy weight", "Overweight", "Obesity")
}

bmi_demos_cdc_weight_class_from_bmip <- function(bmip) {
  # Continuous CDC extended BMI percentile cuts (no gaps between classes):
  # underweight <5; healthy 5 to <85; overweight 85 to <95; obesity >=95.
  # Integer-style cuts (<=4, 5-84, ...) leave (4,5), (84,85), (94,95) unclassified.
  bmip <- as.numeric(bmip)
  dplyr::case_when(
    is.na(bmip) ~ NA_character_,
    bmip < 5 ~ "Underweight",
    bmip < 85 ~ "Healthy weight",
    bmip < 95 ~ "Overweight",
    TRUE ~ "Obesity"
  )
}

bmi_demos_mutate_cdc_weight_class <- function(bmip) {
  raw <- bmi_demos_cdc_weight_class_from_bmip(bmip)
  is_missing <- is.na(bmip) | is.na(raw)
  raw[is_missing] <- "Missing"
  levels <- bmi_demos_cdc_weight_class_levels()
  if (any(is_missing)) {
    levels <- c(levels, "Missing")
  }
  factor(raw, levels = levels)
}

bmi_demos_factor_with_missing <- function(x, map, missing_label = "Missing") {
  raw <- as.character(x)
  labeled <- unname(map[raw])
  is_missing <- is.na(raw) | is.na(labeled)
  labeled[is_missing] <- missing_label
  levels <- if (any(is_missing)) c(unname(map), missing_label) else unname(map)
  factor(labeled, levels = levels)
}

bmi_demos_mutate_characteristics <- function(df, maps = bmi_demos_label_maps(), wave_col = "") {
  out <- df %>%
    dplyr::mutate(
      sex_group = factor(
        ab_g_stc__cohort_sex,
        levels = c("F", "M"),
        labels = c("Female", "Male")
      ),
      age_years = as.numeric(ab_g_dyn__visit_age),
      parental_education = bmi_demos_factor_with_missing(
        ab_g_dyn__cohort_edu__cgs,
        maps$edu
      ),
      household_income = bmi_demos_factor_with_missing(
        ab_g_dyn__cohort_income__hhold__6lvl,
        maps$income
      ),
      race_nih = bmi_demos_factor_with_missing(
        ab_g_stc__cohort_race__nih,
        maps$race
      ),
      ethnicity = bmi_demos_factor_with_missing(
        ab_g_stc__cohort_ethn,
        maps$ethn
      ),
      pubertal_stage = as.numeric(pds_categ_y_p_all),
      bmi_kgm2 = as.numeric(anthro_bmi_corr),
      perc_from_median = as.numeric(perc_median),
      bmip_score = if ("bmip" %in% names(df)) as.numeric(bmip) else NA_real_,
      cdc_weight_class = if ("bmip" %in% names(df)) {
        bmi_demos_mutate_cdc_weight_class(bmip)
      } else {
        NA
      }
    )

  if ("n_visits" %in% names(out)) {
    out$n_study_waves <- factor(
      pmin(as.integer(out$n_visits), 4L),
      levels = 1:4,
      labels = maps$visit_count
    )
  }

  if (nzchar(wave_col) && wave_col %in% names(out)) {
    out$recruitment_wave <- as.factor(out[[wave_col]])
  }

  out
}

bmi_demos_table_var_labels <- function() {
  c(
    age_years = "Age, y",
    sex_group = "Sex, No. (%)",
    sex_display = "Sex, No. (%)",
    race_nih = "Race, No. (%)",
    ethnicity = "Ethnicity, No. (%)",
    household_income = "Annual household income, No. (%)",
    parental_education = "Highest parental education, No. (%)",
    n_study_waves = "No. of study waves with imaging, No. (%)",
    recruitment_wave = "Recruitment wave, No. (%)",
    pubertal_stage = "Puberty",
    bmi_kgm2 = "BMI, kg/m^2",
    perc_from_median = "Percent from median BMI, age/sex",
    bmip_score = "CDC extended BMI percentile",
    cdc_weight_class = "CDC weight class, No. (%)",
    rni_mean = "Global mean RNI"
  )
}

bmi_demos_fmt_iqr <- function(x) {
  dplyr::if_else(
    grepl("\\([^)]+, [^)]+\\)", x),
    gsub("\\(([^,]+), ([^)]+)\\)", "(\\1-\\2)", x, perl = TRUE),
    x
  )
}

bmi_demos_table_specs <- function(table_vars) {
  cont_type_vars <- intersect(
    c("age_years", "pubertal_stage", "bmi_kgm2", "perc_from_median", "bmip_score", "rni_mean"),
    table_vars
  )
  var_types <- as.list(stats::setNames(rep("continuous", length(cont_type_vars)), cont_type_vars))

  statistic <- as.list(c(
    age_years = "{mean} ({sd})",
    pubertal_stage = "{mean} ({sd})",
    bmi_kgm2 = "{mean} ({sd})",
    perc_from_median = "{mean} ({sd})",
    bmip_score = "{mean} ({sd})",
    rni_mean = "{mean} ({sd})"
  )[names(c(
    age_years = "{mean} ({sd})",
    pubertal_stage = "{mean} ({sd})",
    bmi_kgm2 = "{mean} ({sd})",
    perc_from_median = "{mean} ({sd})",
    bmip_score = "{mean} ({sd})",
    rni_mean = "{mean} ({sd})"
  )) %in% table_vars])

  cont_digits <- c(
    age_years = 1,
    pubertal_stage = 1,
    bmi_kgm2 = 1,
    perc_from_median = 1,
    bmip_score = 1,
    rni_mean = 2
  )
  digits <- c(
    as.list(cont_digits[names(cont_digits) %in% table_vars]),
    list(all_categorical() ~ c(0, 1))
  )

  list(
    var_types = var_types,
    statistic = statistic,
    digits = digits
  )
}

bmi_demos_prepend_observation_n_row <- function(
  tbl,
  overall_n,
  level_ns,
  pooled_n = overall_n,
  row_label = "No. of observations, No. (%)"
) {
  stat_cols <- grep("^stat_", names(tbl$table_body), value = TRUE)
  if (!length(stat_cols)) {
    return(tbl)
  }

  n_vals <- rep(NA_character_, length(stat_cols))
  n_vals[[1]] <- bmi_demos_format_count_pct(overall_n, pooled_n)
  if (length(stat_cols) > 1L && length(level_ns)) {
    wave_ns <- unname(level_ns)[seq_len(min(length(level_ns), length(stat_cols) - 1L))]
    for (i in seq_along(wave_ns)) {
      n_vals[[i + 1L]] <- bmi_demos_format_count_pct(wave_ns[[i]], overall_n)
    }
  }

  template <- tbl$table_body[1, , drop = FALSE]
  n_row <- template
  n_row$variable <- "__obs_n__"
  n_row$row_type <- "label"
  n_row$var_type <- NA_character_
  n_row$var_label <- row_label
  n_row$label <- row_label
  for (i in seq_along(stat_cols)) {
    n_row[[stat_cols[[i]]]] <- n_vals[[i]]
  }

  tbl$table_body <- dplyr::bind_rows(n_row, tbl$table_body)
  tbl
}

bmi_demos_build_table <- function(
  df_tbl,
  by_var,
  table_vars,
  header_overall_n,
  header_level_ns,
  header_overall_label = "Overall",
  header_level_prefix = "n",
  header_count_descriptor = NULL,
  include_header_counts = TRUE,
  prepend_observation_n_row = FALSE,
  observation_n_pooled_n = NULL
) {
  format_header_count <- function(label, n, prefix) {
    # Plain text labels (no markdown **); bolding comes from gtsummary/gt styling.
    if (!isTRUE(include_header_counts)) {
      return(as.character(label))
    }
    count_text <- bmi_demos_format_count(n)
    suffix <- if (!is.null(header_count_descriptor) && nzchar(header_count_descriptor)) {
      paste0(" ", header_count_descriptor)
    } else {
      ""
    }
    if (prefix == "N") {
      paste0(label, "\n(N = ", count_text, suffix, ")")
    } else {
      paste0(label, "\n(n = ", count_text, suffix, ")")
    }
  }

  labels <- as.list(bmi_demos_table_var_labels()[names(bmi_demos_table_var_labels()) %in% table_vars])
  specs <- bmi_demos_table_specs(table_vars)

  tbl <- df_tbl %>%
    gtsummary::tbl_summary(
      by = !!rlang::sym(by_var),
      include = dplyr::all_of(table_vars),
      label = labels,
      type = specs$var_types,
      statistic = specs$statistic,
      digits = specs$digits,
      missing = "no",
      percent = "column"
    ) %>%
    gtsummary::add_overall(last = FALSE, col_label = "Overall")

  tbl <- tbl %>%
    gtsummary::bold_labels()

  stat_cols <- grep("^stat_", names(tbl$table_body), value = TRUE)
  tbl$table_body <- tbl$table_body %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(stat_cols),
        ~ bmi_demos_fmt_iqr(.x)
      )
    )

  overall_header <- format_header_count(header_overall_label, header_overall_n, header_level_prefix)

  header_args <- list(label = "Characteristic")
  header_args$stat_0 <- overall_header
  level_names <- names(header_level_ns)
  if (is.null(level_names) || !nzchar(level_names[1])) {
    level_names <- seq_along(header_level_ns)
  }
  for (i in seq_along(header_level_ns)) {
    col_name <- paste0("stat_", i)
    level_label <- level_names[i]
    header_args[[col_name]] <- format_header_count(level_label, header_level_ns[i], "n")
  }

  tbl <- rlang::exec(gtsummary::modify_header, tbl, !!!header_args)

  if (isTRUE(prepend_observation_n_row)) {
    tbl <- bmi_demos_prepend_observation_n_row(
      tbl,
      overall_n = header_overall_n,
      level_ns = header_level_ns,
      pooled_n = if (!is.null(observation_n_pooled_n)) observation_n_pooled_n else header_overall_n
    )
  }

  tbl
}

bmi_demos_build_stacked_sex_wave_table <- function(
  df_tbl,
  table_vars,
  strata_var = "sex_group",
  by_var = "study_wave",
  header_count_descriptor = "observations"
) {
  if (!strata_var %in% names(df_tbl)) {
    stop("strata_var '", strata_var, "' not found in table data.", call. = FALSE)
  }
  if (!by_var %in% names(df_tbl)) {
    stop("by_var '", by_var, "' not found in table data.", call. = FALSE)
  }

  wave_order <- intersect(bmi_demos_session_labels(), unique(as.character(df_tbl[[by_var]])))
  pooled_n <- nrow(df_tbl)

  build_stratum_table <- function(.data, ...) {
    wave_counts <- .data %>%
      dplyr::count(.data[[by_var]], name = "n") %>%
      dplyr::arrange(.data[[by_var]])
    wave_ns <- stats::setNames(wave_counts$n, as.character(wave_counts[[by_var]]))
    wave_ns <- wave_ns[intersect(wave_order, names(wave_ns))]

    stratum_tbl <- bmi_demos_build_table(
      df_tbl = .data %>% dplyr::select(dplyr::all_of(c(by_var, table_vars))),
      by_var = by_var,
      table_vars = table_vars,
      header_overall_n = nrow(.data),
      header_level_ns = wave_ns,
      header_overall_label = "Overall",
      header_level_prefix = "N",
      header_count_descriptor = header_count_descriptor,
      include_header_counts = FALSE,
      prepend_observation_n_row = TRUE,
      observation_n_pooled_n = pooled_n
    )

    stratum_tbl
  }

  strata_counts <- df_tbl %>%
    dplyr::count(.data[[strata_var]], name = "n")

  tbl <- gtsummary::tbl_strata(
    data = df_tbl %>% dplyr::select(dplyr::all_of(c(strata_var, by_var, table_vars))),
    strata = !!rlang::sym(strata_var),
    include = dplyr::all_of(table_vars),
    .combine_with = "tbl_stack",
    .combine_args = list(quiet = TRUE),
    .header = "{strata}",
    .tbl_fun = build_stratum_table
  )

  if (!is.null(tbl$df_strata) && nrow(tbl$df_strata)) {
    strata_col <- intersect(names(tbl$df_strata), c(strata_var, "strata_1"))[[1]]
    header_labels <- vapply(
      seq_len(nrow(tbl$df_strata)),
      function(i) {
        level <- as.character(tbl$df_strata[[strata_col]][[i]])
        n_obs <- strata_counts$n[match(level, as.character(strata_counts[[strata_var]]))]
        suffix <- if (!is.null(header_count_descriptor) && nzchar(header_count_descriptor)) {
          paste0(" ", header_count_descriptor)
        } else {
          ""
        }
        paste0(level, "\n(N = ", bmi_demos_format_count(n_obs), suffix, ")")
      },
      character(1)
    )
    names(header_labels) <- paste0("stacking_header_", seq_along(header_labels))
    tbl <- rlang::exec(gtsummary::modify_spanning_header, tbl, !!!as.list(header_labels))
  }

  tbl
}

bmi_demos_superscript <- function(letter) {
  sup_map <- c(
    a = "\u1d43",
    b = "\u1d47",
    c = "\u1d9c",
    d = "\u1d48",
    e = "\u1d49",
    f = "\u1da0"
  )
  unname(sup_map[[letter]])
}

bmi_demos_clear_gtsummary_footnotes <- function(tbl) {
  empty_rows <- function(x) {
    if (is.null(x) || !nrow(x)) {
      return(x)
    }
    x[0, , drop = FALSE]
  }

  tbl$table_styling$footnote_header <- empty_rows(tbl$table_styling$footnote_header)
  tbl$table_styling$footnote_body <- empty_rows(tbl$table_styling$footnote_body)
  tbl$table_styling$footnote_spanning_header <- empty_rows(tbl$table_styling$footnote_spanning_header)
  tbl$table_styling$abbreviation <- empty_rows(tbl$table_styling$abbreviation)
  tbl
}

bmi_demos_apply_footnote_marks <- function(tbl, entries) {
  for (entry in entries) {
    attach <- entry$attach
    if (is.null(attach)) {
      next
    }
    body_targets <- setdiff(attach, c("title", "abbreviations"))
    if (!length(body_targets)) {
      next
    }
    sup <- bmi_demos_superscript(entry$letter)
    tbl$table_body <- tbl$table_body %>%
      dplyr::mutate(
        label = dplyr::if_else(
          .data$variable %in% body_targets & .data$row_type == "label",
          paste0(.data$label, sup),
          .data$label
        )
      )
  }
  tbl
}

bmi_demos_format_title <- function(table_number, table_title, entries) {
  title_sup <- vapply(
    entries,
    function(entry) {
      if (!is.null(entry$attach) && "title" %in% entry$attach) {
        return(bmi_demos_superscript(entry$letter))
      }
      ""
    },
    character(1)
  )
  paste0("Table ", table_number, ". ", table_title, paste(title_sup, collapse = ""))
}

bmi_demos_footnote_config <- function(
  include_bmip = TRUE,
  extra_entries = NULL,
  stacked_sex_wave_obs_n = FALSE
) {
  footnote_a <- if (isTRUE(stacked_sex_wave_obs_n)) {
    paste(
      "Categorical values are presented as No. (%), with percentages computed within each column",
      "of the corresponding sex block; continuous variables are presented as mean (SD).",
      "For the number of observations row, percentages in the Overall column are based on all",
      "observations in the analytic sample; percentages in study wave columns are based on all",
      "observations within each sex."
    )
  } else {
    paste(
      "Categorical values are presented as No. (%);",
      "continuous variables are presented as mean (SD)."
    )
  }

  # Letters follow first appearance top-to-bottom in Table 1:
  # a = title / No. of observations; b = CDC (if present); c = race; d = ethnicity;
  # later extras (e.g. sample-size note) continue the alphabet.
  entries <- list(
    list(
      letter = "a",
      text = footnote_a,
      attach = if (isTRUE(stacked_sex_wave_obs_n)) c("title", "__obs_n__") else "title"
    )
  )

  if (isTRUE(include_bmip)) {
    entries[[length(entries) + 1]] <- list(
      letter = "b",
      text = paste(
        "CDC extended BMI percentile was derived using CDC age- and sex-specific",
        "growth reference data.",
        "CDC weight classes were derived from CDC extended BMI percentile as underweight",
        "(<5th percentile), healthy weight (5th to <85th percentile), overweight",
        "(85th to <95th percentile), and obesity (>=95th percentile)."
      ),
      attach = c("bmip_score", "cdc_weight_class")
    )
  }

  race_letter <- if (isTRUE(include_bmip)) "c" else "b"
  ethn_letter <- if (isTRUE(include_bmip)) "d" else "c"
  entries <- c(
    entries,
    list(
      list(
        letter = race_letter,
        text = paste(
          "Race (NIH 7 categories) was self-reported by the caregiver.",
          "Descriptive statistics are reported for self-reported race only for",
          "interpretation of sample diversity."
        ),
        attach = "race_nih"
      ),
      list(
        letter = ethn_letter,
        text = "Ethnicity (2 categories) was self-reported by the caregiver.",
        attach = "ethnicity"
      )
    )
  )

  if (!is.null(extra_entries)) {
    if (is.list(extra_entries) && !is.null(extra_entries$letter)) {
      extra_entries <- list(extra_entries)
    }
    used <- vapply(entries, function(e) e$letter, character(1))
    next_i <- length(used)
    alphabet <- letters
    for (i in seq_along(extra_entries)) {
      if (is.null(extra_entries[[i]]$letter) || !nzchar(extra_entries[[i]]$letter)) {
        next_i <- next_i + 1L
        extra_entries[[i]]$letter <- alphabet[[next_i]]
      }
    }
    entries <- c(entries, extra_entries)
  }

  abbreviations <- paste(
    "Abbreviations: BMI, body mass index",
    "(calculated as weight in kilograms divided by height in meters squared); CDC, Centers for",
    "Disease Control and Prevention; FU, Follow-up; GED, Generalized Education Degree; SD,",
    "standard deviation."
  )

  list(entries = entries, abbreviations = abbreviations)
}

# Sort lettered footnotes alphabetically (a, b, c, ...), notes/abbreviations last.
bmi_demos_order_footnote_entries <- function(entries, tbl = NULL) {
  if (!length(entries)) {
    return(entries)
  }
  letters_used <- vapply(entries, function(e) e$letter, character(1))
  entries[order(letters_used, seq_along(entries))]
}

bmi_demos_format_footnotes <- function(footnote_config) {
  entry_lines <- if (length(footnote_config$entries)) {
    vapply(
      footnote_config$entries,
      function(entry) {
        paste0(bmi_demos_superscript(entry$letter), " ", entry$text)
      },
      character(1)
    )
  } else {
    character()
  }

  general_lines <- if (!is.null(footnote_config$general_notes) && length(footnote_config$general_notes)) {
    vapply(
      footnote_config$general_notes,
      function(note) {
        paste0("Note. ", note)
      },
      character(1)
    )
  } else {
    character()
  }

  # Abbreviations always last, below lettered footnotes / notes.
  abbrev <- footnote_config$abbreviations
  if (is.null(abbrev)) {
    abbrev <- character()
  }
  c(entry_lines, general_lines, abbrev)
}

bmi_demos_build_gt_table <- function(tbl) {
  tbl %>%
    gtsummary::as_gt() %>%
    gt::tab_options(
      source_notes.font.size = gt::px(10),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      heading.border.bottom.style = "solid",
      column_labels.border.top.style = "solid",
      column_labels.border.bottom.style = "solid"
    )
}

bmi_table_style_config <- function() {
  list(
    font_family = Sys.getenv("BMI_IMAGING_TABLE_FONT_FAMILY", "Arial"),
    title_font_size = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_TITLE_FONT_SIZE", "12")),
    part_label_font_size = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_PART_FONT_SIZE", "11")),
    footnote_font_size = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_FOOTNOTE_FONT_SIZE", "8")),
    body_font_size = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_BODY_FONT_SIZE", "10")),
    label_font_size = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_LABEL_FONT_SIZE", "9")),
    row_padding_px = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_ROW_PADDING_PX", "3")),
    brain_region_width_px = as.numeric(Sys.getenv("BMI_IMAGING_TABLE_BRAIN_REGION_WIDTH_PX", "140")),
    part_page_break = tolower(Sys.getenv("BMI_IMAGING_TABLE_PART_PAGE_BREAK", "0")) %in%
      c("1", "true", "yes", "on")
  )
}

bmi_table_docx_add_text <- function(
  doc,
  text,
  font_family = "Arial",
  font_size = 12,
  bold = FALSE,
  spacing_after = 0
) {
  if (!nzchar(text)) {
    return(officer::body_add_par(doc, "", style = "Normal"))
  }

  prop <- officer::fp_text(
    font.family = font_family,
    font.size = font_size,
    bold = bold
  )
  fpar_obj <- if (spacing_after > 0) {
    officer::fpar(
      officer::ftext(text, prop),
      fp_p = officer::fp_par(padding.bottom = spacing_after)
    )
  } else {
    officer::fpar(officer::ftext(text, prop))
  }
  officer::body_add_fpar(doc, fpar_obj)
}

bmi_table_apply_gt_table_style <- function(gt_tbl, table_style = bmi_table_style_config()) {
  font_family <- table_style$font_family
  body_font_size <- table_style$body_font_size
  label_font_size <- table_style$label_font_size
  title_font_size <- table_style$title_font_size
  row_padding_px <- table_style$row_padding_px

  gt::tab_options(
    gt_tbl,
    table.font.size = body_font_size,
    table.font.names = font_family,
    heading.title.font.size = title_font_size,
    column_labels.font.size = label_font_size,
    data_row.padding = row_padding_px,
    table.width = gt::pct(100),
    table.layout = "auto"
  )
}

bmi_demos_export_css <- function(title_font_size_pt = 12) {
  paste(
    ".bmi-table-title {",
    paste0("  font-size: ", title_font_size_pt, "pt;"),
    "  font-weight: normal;",
    "  margin: 0 0 8px 0;",
    "}",
    "@media print {",
    "  table.gt_table thead {",
    "    display: table-row-group;",
    "  }",
    "}",
    ".gt_table_container {",
    "  overflow: visible;",
    "}",
    sep = "\n"
  )
}

bmi_demos_strip_word_header_repeat <- function(word_xml) {
  word_xml <- gsub("<w:tblHeader\\s*/>", "", word_xml, perl = TRUE)
  gsub("<w:tblHeader>\\s*</w:tblHeader>", "", word_xml, perl = TRUE)
}

bmi_table_word_escape_xml_text <- function(text) {
  text <- gsub("&", "&amp;", text, fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  text <- gsub(">", "&gt;", text, fixed = TRUE)
  text
}

bmi_table_word_escape_regex <- function(text) {
  gsub("([][\\\\()*+?{}|^$.])", "\\\\\\1", text, perl = TRUE)
}

bmi_table_gt_spanner_line_lists <- function(gt_tbl) {
  spanners <- gt_tbl[["_spanners"]]
  if (is.null(spanners) || !nrow(spanners)) {
    return(list())
  }

  line_lists <- lapply(spanners$spanner_label, function(lbl) {
    if (inherits(lbl, "md") || inherits(lbl, "html")) {
      lbl <- as.character(lbl)
    }
    lines <- strsplit(as.character(lbl), "\n", fixed = TRUE)[[1]]
    lines[nzchar(lines)]
  })
  line_lists[vapply(line_lists, length, integer(1)) >= 2L]
}

bmi_table_gt_column_label_line_lists <- function(gt_tbl) {
  boxhead <- gt_tbl[["_boxhead"]]
  if (is.null(boxhead) || !nrow(boxhead)) {
    return(list())
  }

  line_lists <- lapply(boxhead$column_label, function(lbl) {
    if (inherits(lbl, "md") || inherits(lbl, "html")) {
      lbl <- as.character(lbl)
    }
    lines <- strsplit(as.character(lbl), "\n", fixed = TRUE)[[1]]
    lines[nzchar(lines)]
  })
  line_lists[vapply(line_lists, length, integer(1)) >= 2L]
}

bmi_table_gt_multiline_label_line_lists <- function(gt_tbl) {
  c(bmi_table_gt_spanner_line_lists(gt_tbl), bmi_table_gt_column_label_line_lists(gt_tbl))
}

bmi_table_word_insert_flat_label_breaks <- function(word_xml, line_lists) {
  for (lines in line_lists) {
    if (length(lines) < 2L) {
      next
    }

    escaped <- vapply(lines, bmi_table_word_escape_xml_text, character(1))
    flat <- paste(escaped, collapse = " ")
    wt_pat <- paste0(
      "<w:t([^>]*)>\\s*",
      bmi_table_word_escape_regex(flat),
      "\\s*</w:t>"
    )

    while (TRUE) {
      wt_m <- regexpr(wt_pat, word_xml, perl = TRUE)
      if (wt_m < 0L) {
        break
      }

      pre <- substr(word_xml, 1L, wt_m - 1L)
      r_starts <- gregexpr("<w:r(?:\\s|>)", pre, perl = TRUE)[[1]]
      r_starts <- r_starts[r_starts > 0L]
      if (!length(r_starts)) {
        break
      }
      r_start <- max(r_starts)

      post_wt_end <- wt_m + attr(wt_m, "match.length") - 1L
      r_end_rel <- regexpr("</w:r>", substr(word_xml, post_wt_end, nchar(word_xml)), fixed = TRUE)
      if (r_end_rel < 0L) {
        break
      }
      r_end <- post_wt_end + r_end_rel + nchar("</w:r>") - 1L

      matched_r <- substr(word_xml, r_start, r_end)
      r_pat <- paste0(
        "<w:r([^>]*)>\\s*<w:rPr>([\\s\\S]*?)</w:rPr>\\s*<w:t([^>]*)>\\s*",
        bmi_table_word_escape_regex(flat),
        "\\s*</w:t>\\s*</w:r>"
      )
      groups <- regmatches(matched_r, regexec(r_pat, matched_r, perl = TRUE))[[1]]
      if (length(groups) < 4L) {
        break
      }

      r_attrs <- groups[[2]]
      rpr <- groups[[3]]
      t_attrs <- groups[[4]]

      rebuilt <- escaped[[1]]
      if (length(escaped) > 1L) {
        for (line in escaped[-1L]) {
          rebuilt <- paste0(
            rebuilt,
            "</w:t></w:r><w:br/><w:r",
            r_attrs,
            "><w:rPr>",
            rpr,
            "</w:rPr><w:t",
            t_attrs,
            ">",
            line
          )
        }
      }
      rebuilt <- paste0(
        "<w:r",
        r_attrs,
        "><w:rPr>",
        rpr,
        "</w:rPr><w:t",
        t_attrs,
        ">",
        rebuilt,
        "</w:t></w:r>"
      )

      word_xml <- paste0(
        substr(word_xml, 1L, r_start - 1L),
        rebuilt,
        substr(word_xml, r_end + 1L, nchar(word_xml))
      )
    }
  }

  word_xml
}

bmi_table_word_insert_multiline_label_breaks <- function(word_xml, spanner_line_lists, column_line_lists) {
  if (length(spanner_line_lists)) {
    word_xml <- bmi_table_word_insert_spanner_breaks(word_xml, spanner_line_lists)
  }
  if (length(column_line_lists)) {
    word_xml <- bmi_table_word_insert_flat_label_breaks(word_xml, column_line_lists)
  }
  word_xml
}

bmi_table_word_remove_spacer_runs <- function(word_xml) {
  spacer_run <- paste0(
    "<w:r xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\">",
    "\\s*<w:rPr>[\\s\\S]*?</w:rPr>\\s*",
    "<w:t xml:space=\"preserve\"> </w:t>\\s*</w:r>"
  )

  gsub(spacer_run, "", word_xml, perl = TRUE)
}

bmi_table_word_insert_spanner_breaks <- function(word_xml, spanner_line_lists) {
  word_xml <- bmi_table_word_remove_spacer_runs(word_xml)

  for (lines in spanner_line_lists) {
    if (length(lines) < 2L) {
      next
    }

    for (i in seq_len(length(lines) - 1L)) {
      from <- bmi_table_word_escape_xml_text(lines[[i]])
      to <- bmi_table_word_escape_xml_text(lines[[i + 1L]])
      pat <- paste0(
        "(", bmi_table_word_escape_regex(from), "</w:t>\\s*</w:r>)",
        "\\s*",
        "(<w:r[^>]*>\\s*<w:rPr>[\\s\\S]*?</w:rPr>\\s*<w:t[^>]*>\\s*",
        bmi_table_word_escape_regex(to),
        "</w:t>\\s*</w:r>)"
      )
      repl <- "\\1<w:br/>\\2"

      repeat {
        new_xml <- gsub(pat, repl, word_xml, perl = TRUE)
        if (identical(new_xml, word_xml)) {
          break
        }
        word_xml <- new_xml
      }
    }
  }

  word_xml
}

bmi_demos_gt_to_word_xml <- function(gt_tbl) {
  word_xml <- gt::as_word(gt_tbl, split = TRUE)
  if (length(word_xml) > 1L) {
    word_xml <- paste(word_xml, collapse = "")
  }
  word_xml <- bmi_demos_strip_word_header_repeat(word_xml)
  spanner_line_lists <- bmi_table_gt_spanner_line_lists(gt_tbl)
  column_line_lists <- bmi_table_gt_column_label_line_lists(gt_tbl)
  if (length(spanner_line_lists) || length(column_line_lists)) {
    word_xml <- bmi_table_word_insert_multiline_label_breaks(
      word_xml,
      spanner_line_lists = spanner_line_lists,
      column_line_lists = column_line_lists
    )
  }
  word_xml
}

bmi_demos_write_html <- function(
  gt_tbl,
  formatted_title,
  path,
  document_title = "Demographics table"
) {
  table_html <- gt::as_raw_html(gt_tbl, inline_css = TRUE)
  html_doc <- paste0(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>", document_title, "</title>",
    "<style>", bmi_demos_export_css(), "</style>",
    "</head>",
    "<body>",
    "<div class=\"bmi-table-title\">", formatted_title, "</div>",
    table_html,
    "</body>",
    "</html>"
  )
  writeLines(html_doc, path, useBytes = TRUE)
  path
}

bmi_demos_write_docx <- function(gt_tbl, formatted_title, path) {
  if (!requireNamespace("officer", quietly = TRUE)) {
    warning("Package 'officer' is required for manuscript docx export; using gt::gtsave().", call. = FALSE)
    gt::gtsave(gt_tbl, path)
    return(path)
  }

  word_xml <- bmi_demos_gt_to_word_xml(gt_tbl)
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, formatted_title, style = "Normal")
  doc <- officer::body_add_par(doc, "", style = "Normal")
  doc <- officer::body_add_xml(doc, word_xml)
  print(doc, target = path)
  path
}

bmi_table_save_gt_tables_docx <- function(table_bundles, path) {
  if (!length(table_bundles)) {
    return(invisible(NULL))
  }
  if (!requireNamespace("officer", quietly = TRUE)) {
    warning("Package 'officer' is required for multi-table docx export.", call. = FALSE)
    return(invisible(NULL))
  }

  doc <- officer::read_docx()
  for (i in seq_along(table_bundles)) {
    bundle <- table_bundles[[i]]
    gt_tbl <- bundle$gt_tbl
    if (nzchar(bundle$footnote_text)) {
      gt_tbl <- gt_tbl %>% gt::tab_source_note(source_note = bundle$footnote_text)
    }
    word_xml <- bmi_demos_gt_to_word_xml(gt_tbl)
    if (i > 1L) {
      doc <- officer::body_add_break(doc)
    }
    doc <- officer::body_add_par(doc, bundle$formatted_title, style = "Normal")
    doc <- officer::body_add_par(doc, "", style = "Normal")
    doc <- officer::body_add_xml(doc, word_xml)
  }
  print(doc, target = path)
  path
}

bmi_table_write_multipart_gt_docx <- function(
  document_title,
  table_parts,
  footnote_text,
  path,
  table_style = bmi_table_style_config()
) {
  if (!requireNamespace("officer", quietly = TRUE)) {
    warning("Package 'officer' is required for multipart docx export.", call. = FALSE)
    return(invisible(NULL))
  }

  doc <- officer::read_docx()
  doc <- bmi_table_docx_add_text(
    doc,
    document_title,
    font_family = table_style$font_family,
    font_size = table_style$title_font_size,
    spacing_after = 6
  )

  for (i in seq_along(table_parts)) {
    part <- table_parts[[i]]
    word_xml <- bmi_demos_gt_to_word_xml(part$gt_tbl)
    if (i > 1L) {
      if (isTRUE(table_style$part_page_break)) {
        doc <- officer::body_add_break(doc)
      } else {
        doc <- bmi_table_docx_add_text(
          doc,
          "",
          font_family = table_style$font_family,
          font_size = table_style$body_font_size
        )
      }
    }
    doc <- bmi_table_docx_add_text(
      doc,
      part$part_label,
      font_family = table_style$font_family,
      font_size = table_style$part_label_font_size,
      spacing_after = 4
    )
    doc <- officer::body_add_xml(doc, word_xml)
  }

  if (nzchar(footnote_text)) {
    doc <- bmi_table_docx_add_text(
      doc,
      "",
      font_family = table_style$font_family,
      font_size = table_style$footnote_font_size
    )
    for (line in strsplit(footnote_text, "\n\n", fixed = TRUE)[[1]]) {
      doc <- bmi_table_docx_add_text(
        doc,
        line,
        font_family = table_style$font_family,
        font_size = table_style$footnote_font_size,
        spacing_after = 4
      )
    }
  }

  print(doc, target = path)
  path
}

bmi_table_write_multipart_gt_html <- function(
  document_title,
  table_parts,
  footnote_text,
  path,
  document_name = "table",
  table_style = bmi_table_style_config()
) {
  part_blocks <- vapply(
    table_parts,
    function(part) {
      paste0(
        "<div class=\"bmi-table-part-label\">", part$part_label, "</div>",
        gt::as_raw_html(part$gt_tbl, inline_css = TRUE)
      )
    },
    character(1)
  )
  footnote_block <- if (nzchar(footnote_text)) {
    paste0(
      "<div class=\"bmi-table-footnotes\">",
      paste(
        vapply(
          strsplit(footnote_text, "\n\n", fixed = TRUE)[[1]],
          function(line) paste0("<p>", line, "</p>"),
          character(1)
        ),
        collapse = ""
      ),
      "</div>"
    )
  } else {
    ""
  }

  html_doc <- paste0(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>", document_name, "</title>",
    "<style>",
    bmi_demos_export_css(title_font_size_pt = table_style$title_font_size),
    "body { font-family: ", table_style$font_family, ", sans-serif; }",
    ".bmi-table-part-label {",
    "  font-size: ", table_style$part_label_font_size, "pt;",
    "  font-weight: normal;",
    "  margin: 16px 0 8px 0;",
    "}",
    ".bmi-table-footnotes {",
    "  font-size: ", table_style$footnote_font_size, "pt;",
    "  margin-top: 16px;",
    "}",
    ".bmi-table-footnotes p { margin: 0 0 8px 0; }",
    "</style>",
    "</head>",
    "<body>",
    "<div class=\"bmi-table-title\">", document_title, "</div>",
    paste(part_blocks, collapse = ""),
    footnote_block,
    "</body>",
    "</html>"
  )
  writeLines(html_doc, path, useBytes = TRUE)
  path
}

bmi_table_save_multipart_gt_table <- function(
  table_number,
  table_title,
  footnote_config,
  table_parts,
  out_dir,
  file_prefix,
  table_style = bmi_table_style_config()
) {
  if (!length(table_parts)) {
    stop("table_parts must contain at least one table part.", call. = FALSE)
  }

  formatted_title <- bmi_demos_format_title(
    table_number,
    table_title,
    footnote_config$entries
  )
  footnote_lines <- bmi_demos_format_footnotes(footnote_config)
  footnote_text <- paste(footnote_lines, collapse = "\n\n")

  paths <- list(
    docx = file.path(out_dir, paste0(file_prefix, ".docx")),
    html = file.path(out_dir, paste0(file_prefix, ".html"))
  )

  save_gt_output <- function(save_fn, path) {
    tryCatch(
      {
        save_fn(path)
        path
      },
      error = function(e) {
        warning("Could not save ", basename(path), ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  for (i in seq_along(table_parts)) {
    part <- table_parts[[i]]
    part_paths <- list(
      csv = file.path(out_dir, paste0(file_prefix, "_", part$part_id, ".csv"))
    )
    display_df <- bmi_demos_plain_export_df(part$display_df)
    label_col <- names(display_df)[1]
    title_row <- display_df[1, , drop = FALSE]
    title_row[] <- NA_character_
    title_row[[label_col]][[1]] <- paste(formatted_title, part$part_label, sep = " | ")
    readr::write_csv(dplyr::bind_rows(title_row, display_df), part_paths$csv)
    if (requireNamespace("writexl", quietly = TRUE)) {
      xlsx_path <- file.path(out_dir, paste0(file_prefix, "_", part$part_id, ".xlsx"))
      writexl::write_xlsx(display_df, xlsx_path)
      part_paths$xlsx <- xlsx_path
    }
    table_parts[[i]]$paths <- part_paths
  }

  docx_written <- save_gt_output(
    function(path) {
      bmi_table_write_multipart_gt_docx(
        document_title = formatted_title,
        table_parts = table_parts,
        footnote_text = footnote_text,
        path = path,
        table_style = table_style
      )
    },
    paths$docx
  )
  html_written <- save_gt_output(
    function(path) {
      bmi_table_write_multipart_gt_html(
        document_title = formatted_title,
        table_parts = table_parts,
        footnote_text = footnote_text,
        path = path,
        document_name = file_prefix,
        table_style = table_style
      )
    },
    paths$html
  )

  cat("\n", file_prefix, " saved to:\n", sep = "")
  for (part in table_parts) {
    cat("  ", part$paths$csv, "\n")
    if (!is.null(part$paths$xlsx) && file.exists(part$paths$xlsx)) {
      cat("  ", part$paths$xlsx, "\n")
    }
  }
  if (!is.null(html_written)) cat("  ", paths$html, "\n")
  if (!is.null(docx_written)) cat("  ", paths$docx, "\n")

  invisible(list(
    formatted_title = formatted_title,
    footnote_text = footnote_text,
    paths = paths,
    table_parts = table_parts
  ))
}

# Strip gtsummary markdown emphasis (**bold**, __italic__) from flat exports.
bmi_demos_plain_export_df <- function(df) {
  strip_md <- function(x) {
    x <- as.character(x)
    x <- gsub("\\*\\*", "", x, perl = TRUE)
    x <- gsub("__", "", x, perl = TRUE)
    x
  }
  out <- df
  names(out) <- strip_md(names(out))
  for (nm in names(out)) {
    if (is.character(out[[nm]]) || is.factor(out[[nm]])) {
      out[[nm]] <- strip_md(out[[nm]])
    }
  }
  out
}

bmi_table_save_gt_table <- function(
  gt_tbl,
  table_number,
  table_title,
  footnote_config,
  out_dir,
  file_prefix,
  csv_df = NULL,
  write_docx = TRUE
) {
  entries <- footnote_config$entries
  formatted_title <- bmi_demos_format_title(table_number, table_title, entries)
  footnote_lines <- bmi_demos_format_footnotes(footnote_config)
  footnote_text <- paste(footnote_lines, collapse = "\n\n")

  if (!inherits(gt_tbl, "gt_tbl")) {
    stop("gt_tbl must be a gt table object.", call. = FALSE)
  }
  if (nzchar(footnote_text)) {
    gt_tbl <- gt_tbl %>% gt::tab_source_note(source_note = footnote_text)
  }

  paths <- list(
    csv = file.path(out_dir, paste0(file_prefix, ".csv")),
    html = file.path(out_dir, paste0(file_prefix, ".html")),
    docx = file.path(out_dir, paste0(file_prefix, ".docx")),
    xlsx = file.path(out_dir, paste0(file_prefix, ".xlsx"))
  )

  if (is.null(csv_df)) {
    csv_df <- gt_tbl[["_data"]]
  }
  csv_df <- bmi_demos_plain_export_df(csv_df)
  label_col <- names(csv_df)[1]
  title_row <- csv_df[1, , drop = FALSE]
  title_row[] <- NA_character_
  title_row[[label_col]][[1]] <- formatted_title
  readr::write_csv(dplyr::bind_rows(title_row, csv_df), paths$csv)

  save_gt_output <- function(save_fn, path) {
    tryCatch(
      {
        save_fn(path)
        path
      },
      error = function(e) {
        warning("Could not save ", basename(path), ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  html_written <- save_gt_output(
    function(path) {
      bmi_demos_write_html(
        gt_tbl,
        formatted_title,
        path,
        document_title = file_prefix
      )
    },
    paths$html
  )
  docx_written <- NULL
  if (isTRUE(write_docx)) {
    docx_written <- save_gt_output(
      function(path) bmi_demos_write_docx(gt_tbl, formatted_title, path),
      paths$docx
    )
  }
  xlsx_written <- NULL
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(csv_df, paths$xlsx)
    xlsx_written <- paths$xlsx
  }

  cat("\n", file_prefix, " saved to:\n", sep = "")
  cat("  ", paths$csv, "\n")
  if (!is.null(html_written)) cat("  ", paths$html, "\n")
  if (!is.null(docx_written)) cat("  ", paths$docx, "\n")
  if (!is.null(xlsx_written)) cat("  ", paths$xlsx, "\n")

  invisible(list(gt = gt_tbl, paths = paths, formatted_title = formatted_title))
}

bmi_demos_save_table <- function(
  tbl,
  table_number,
  table_title,
  footnote_config,
  out_dir,
  file_prefix
) {
  entries <- footnote_config$entries
  tbl <- bmi_demos_apply_footnote_marks(tbl, entries)
  tbl <- bmi_demos_clear_gtsummary_footnotes(tbl)
  footnote_config$entries <- bmi_demos_order_footnote_entries(entries, tbl)

  gt_tbl <- bmi_demos_build_gt_table(tbl)
  tbl_df <- gtsummary::as_tibble(tbl)

  saved <- bmi_table_save_gt_table(
    gt_tbl = gt_tbl,
    table_number = table_number,
    table_title = table_title,
    footnote_config = footnote_config,
    out_dir = out_dir,
    file_prefix = file_prefix,
    csv_df = tbl_df
  )

  invisible(list(tbl = tbl, gt = saved$gt, paths = saved$paths))
}


bmi_demos_to_participant_level <- function(df, index_session = "baseline") {
  session_levels <- bmi_demos_session_levels()
  df <- df %>%
    dplyr::mutate(
      participant_id = as.character(participant_id),
      session_id = as.character(session_id),
      session_ord = match(session_id, session_levels)
    )

  if (index_session == "baseline") {
    idx <- df %>%
      dplyr::filter(session_id == "ses-00A") %>%
      dplyr::group_by(participant_id) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    missing_ids <- setdiff(unique(df$participant_id), idx$participant_id)
    if (length(missing_ids)) {
      fallback <- df %>%
        dplyr::filter(participant_id %in% missing_ids) %>%
        dplyr::group_by(participant_id) %>%
        dplyr::arrange(session_ord) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      idx <- dplyr::bind_rows(idx, fallback)
    }
  } else if (index_session == "first") {
    idx <- df %>%
      dplyr::group_by(participant_id) %>%
      dplyr::arrange(session_ord) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  } else {
    stop("index_session must be 'baseline' or 'first'", call. = FALSE)
  }

  visit_counts <- df %>%
    dplyr::distinct(participant_id, session_id) %>%
    dplyr::count(participant_id, name = "n_visits")

  idx %>%
    dplyr::left_join(visit_counts, by = "participant_id") %>%
    dplyr::select(-session_ord)
}

bmi_demos_format_count <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

bmi_demos_format_count_pct <- function(n, denominator) {
  n <- as.integer(n)
  denominator <- as.integer(denominator)
  count_text <- bmi_demos_format_count(n)
  if (!length(denominator) || is.na(denominator) || denominator < 1L) {
    return(count_text)
  }
  pct <- 100 * n / denominator
  digits <- if (denominator < 100L) 0L else 1L
  pct_text <- format(round(pct, digits), nsmall = digits, trim = TRUE)
  paste0(count_text, " (", pct_text, "%)")
}

bmi_demos_sample_size_summary <- function(df, sex_col = "sex_group", wave_col = "study_wave") {
  sex_vals <- df[[sex_col]]
  wave_counts_by_sex <- NULL
  if (wave_col %in% names(df)) {
    wave_counts_by_sex <- df %>%
      dplyr::count(.data[[sex_col]], .data[[wave_col]], name = "n") %>%
      dplyr::arrange(.data[[sex_col]], .data[[wave_col]])
  }

  list(
    n_obs = nrow(df),
    n_obs_female = sum(sex_vals == "Female", na.rm = TRUE),
    n_obs_male = sum(sex_vals == "Male", na.rm = TRUE),
    n_participants = dplyr::n_distinct(df$participant_id),
    n_participants_female = dplyr::n_distinct(df$participant_id[sex_vals == "Female"]),
    n_participants_male = dplyr::n_distinct(df$participant_id[sex_vals == "Male"]),
    wave_counts = if (wave_col %in% names(df)) {
      df %>%
        dplyr::count(.data[[wave_col]], name = "n") %>%
        dplyr::arrange(.data[[wave_col]])
    } else {
      NULL
    },
    wave_counts_by_sex = wave_counts_by_sex
  )
}

bmi_demos_format_sex_wave_counts <- function(wave_counts_by_sex, sex_col = "sex_group", wave_col = "study_wave") {
  if (is.null(wave_counts_by_sex) || !nrow(wave_counts_by_sex)) {
    return("")
  }

  sex_levels <- unique(as.character(wave_counts_by_sex[[sex_col]]))
  vapply(
    sex_levels,
    function(sex_label) {
      rows <- wave_counts_by_sex[as.character(wave_counts_by_sex[[sex_col]]) == sex_label, , drop = FALSE]
      wave_text <- paste(
        vapply(
          seq_len(nrow(rows)),
          function(i) {
            paste0(rows[[wave_col]][[i]], " ", bmi_demos_format_count(rows$n[[i]]))
          },
          character(1)
        ),
        collapse = ", "
      )
      paste0(sex_label, ": ", wave_text)
    },
    character(1)
  ) %>%
    paste(collapse = "; ")
}

bmi_demos_format_wave_counts <- function(wave_counts) {
  if (is.null(wave_counts) || !nrow(wave_counts)) {
    return("")
  }
  wave_col <- setdiff(names(wave_counts), "n")[[1]]
  paste(
    vapply(
      seq_len(nrow(wave_counts)),
      function(i) {
        paste0(wave_counts[[wave_col]][[i]], ": ", bmi_demos_format_count(wave_counts$n[[i]]))
      },
      character(1)
    ),
    collapse = "; "
  )
}
