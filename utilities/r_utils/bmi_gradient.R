# Gradient table helpers — reads gradient_analysis CSV exports

bmi_gradient_default_rois <- function() {
  c(
    "Bilateral-Accumbens-area",
    "Bilateral-ATR",
    "Bilateral-Caudate",
    "Bilateral-CgC",
    "Fmin",
    "Bilateral-Pallidum",
    "Bilateral-CgH",
    "Bilateral-Putamen",
    "Bilateral-Thalamus-Proper"
  )
}

bmi_gradient_clean_roi_name <- function(name) {
  name <- gsub("Bilateral-", "", name, fixed = TRUE)
  name <- gsub("ATR", "Anterior thalamic radiations", name, fixed = TRUE)
  name <- gsub("CgC", "Cingulate cingulum", name, fixed = TRUE)
  name <- gsub("CgH", "Parahippocampal cingulum", name, fixed = TRUE)
  name <- gsub("-Proper", "", name, fixed = TRUE)
  name <- gsub("Accumbens-area", "Accumbens", name, fixed = TRUE)
  name <- gsub("CC_Anterior", "Corpus callosum (anterior)", name, fixed = TRUE)
  name <- gsub("Fmin", "Forceps minor", name, fixed = TRUE)
  name <- gsub("Fmaj", "Forceps major", name, fixed = TRUE)
  name <- gsub("VentralDC", "Ventral diencephalon", name, fixed = TRUE)
  name
}

bmi_gradient_warn_missing_rois <- function(wide, rois, context = NULL) {
  present <- unique(as.character(wide$roi))
  missing <- setdiff(rois, present)
  if (!length(missing)) {
    return(invisible(character(0)))
  }
  prefix <- if (is.null(context) || !nzchar(context)) {
    "Requested ROIs missing from gradient CSV"
  } else {
    paste0("Requested ROIs missing from gradient CSV (", context, ")")
  }
  warning(
    prefix, ": ", paste(missing, collapse = ", "),
    "\n  Re-run tables/table2_gradient_csvs.ipynb with SAVE_CSV = TRUE after s9 ROI outputs exist.",
    call. = FALSE
  )
  invisible(missing)
}

bmi_gradient_cdc_segment_labels <- function() {
  c(
    cdc_underweight = "Underweight (<5th)",
    cdc_healthy = "Healthy (5th to <85th)",
    cdc_overweight = "Overweight (85th to <95th)",
    cdc_obesity = "Obesity (>=95th)"
  )
}

bmi_gradient_cdc_class_parts <- function() {
  list(
    cdc_underweight = list(name = "Underweight", pct = "(< 5th)"),
    cdc_healthy = list(name = "Healthy", pct = "(5th to <85th)"),
    cdc_overweight = list(name = "Overweight", pct = "(85th to <95th)"),
    cdc_obesity = list(name = "Obesity", pct = "(>= 95th)")
  )
}

bmi_gradient_format_axis_range <- function(lo, hi, digits = 1) {
  paste0(
    "[",
    format(round(lo, digits), nsmall = digits, trim = TRUE),
    ", ",
    format(round(hi, digits), nsmall = digits, trim = TRUE),
    "]"
  )
}

bmi_gradient_load_cdc_mapping <- function(grad_dir, cdc_method, correction) {
  candidates <- c(
    file.path(grad_dir, paste0("cdc_bmip_to_perc_median_mapping_", cdc_method, "_", correction, ".csv")),
    file.path(grad_dir, paste0("cdc_bmip_to_perc_median_mapping_", correction, ".csv"))
  )
  for (path in candidates) {
    if (file.exists(path)) {
      return(readr::read_csv(path, show_col_types = FALSE))
    }
  }
  NULL
}

bmi_gradient_cdc_spanner_label <- function(seg, pheno, cdc_mapping = NULL, sex = NULL) {
  parts <- bmi_gradient_cdc_class_parts()[[seg]]
  if (is.null(parts)) {
    return(gsub("_", " ", seg, fixed = TRUE))
  }

  lines <- c(parts$name, parts$pct)
  if (identical(pheno, "perc_median") && !is.null(cdc_mapping) && nrow(cdc_mapping)) {
    class_name <- sub("^cdc_", "", seg)
    rows <- cdc_mapping[cdc_mapping$cdc_class == class_name, , drop = FALSE]
    if (!is.null(sex) && nrow(rows)) {
      sex_key <- tolower(sex)
      sex_letter <- switch(sex_key, female = "F", male = "M", substr(sex_key, 1L, 1L))
      if ("designmat_sex" %in% names(rows)) {
        rows <- rows[rows$designmat_sex == sex_key, , drop = FALSE]
      }
      if (!nrow(rows) && "sex" %in% names(rows)) {
        rows <- rows[toupper(rows$sex) == toupper(sex_letter), , drop = FALSE]
      }
    }
    if (nrow(rows)) {
      if (!is.null(sex) && nrow(rows) == 1L) {
        lo <- rows$perc_median_range_low[[1]]
        hi <- rows$perc_median_range_high[[1]]
      } else {
        lo <- min(rows$perc_median_range_low, na.rm = TRUE)
        hi <- max(rows$perc_median_range_high, na.rm = TRUE)
      }
      if (is.finite(lo) && is.finite(hi)) {
        lines <- c(lines, bmi_gradient_format_axis_range(lo, hi))
      }
    }
  }
  paste(lines, collapse = "\n")
}

bmi_gradient_cdc_segment_ranges <- function() {
  list(
    cdc_underweight = c(0, 4),
    cdc_healthy = c(5, 84),
    cdc_overweight = c(85, 94),
    cdc_obesity = c(95, 100)
  )
}

bmi_gradient_results_dir <- function(ceph_root, analysis_dir, results_run) {
  file.path(ceph_root, analysis_dir, "results", results_run, "gradient_analysis")
}

bmi_gradient_load_wide <- function(path) {
  bmi_assert_file(path, "Gradient summary wide CSV")
  readr::read_csv(path, show_col_types = FALSE)
}

bmi_gradient_cdc_weight_class_comparisons <- function() {
  list(
    list(
      label = "obesity vs healthy weight",
      numerator = "cdc_obesity",
      denominator = "cdc_healthy"
    ),
    list(
      label = "obesity vs underweight",
      numerator = "cdc_obesity",
      denominator = "cdc_underweight"
    ),
    list(
      label = "overweight vs healthy weight",
      numerator = "cdc_overweight",
      denominator = "cdc_healthy"
    ),
    list(
      label = "overweight vs underweight",
      numerator = "cdc_overweight",
      denominator = "cdc_underweight"
    )
  )
}

bmi_gradient_collect_bmip_fold_changes <- function(
  wide,
  rois = bmi_gradient_default_rois(),
  pheno = "bmip",
  modality = "rni",
  value_col = "mean_signed",
  comparisons = bmi_gradient_cdc_weight_class_comparisons(),
  positive_gradients_only = TRUE
) {
  mod <- wide %>%
    dplyr::filter(
      .data$pheno == .env$pheno,
      .data$modality == .env$modality,
      .data$roi %in% .env$rois
    )

  if (!nrow(mod)) {
    return(tibble::tibble())
  }

  rows <- list()
  for (comp in comparisons) {
    num_col <- paste0(comp$numerator, "_", value_col)
    den_col <- paste0(comp$denominator, "_", value_col)
    if (!num_col %in% names(mod) || !den_col %in% names(mod)) {
      next
    }

    chunk <- mod %>%
      dplyr::transmute(
        comparison = comp$label,
        sex = .data$sex,
        roi = .data$roi,
        brain_region = vapply(.data$roi, bmi_gradient_clean_roi_name, character(1)),
        numerator = .data[[num_col]],
        denominator = .data[[den_col]],
        fold_change = .data[[num_col]] / .data[[den_col]]
      ) %>%
      dplyr::filter(
        is.finite(.data$fold_change),
        .data$denominator != 0
      )

    if (isTRUE(positive_gradients_only)) {
      chunk <- chunk %>%
        dplyr::filter(.data$numerator > 0, .data$denominator > 0)
    }

    rows[[length(rows) + 1L]] <- chunk
  }

  dplyr::bind_rows(rows)
}

bmi_gradient_format_gradient_sci <- function(x) {
  vapply(
    as.numeric(x),
    function(val) {
      if (!is.finite(val)) {
        return("NA")
      }
      format(signif(val, 3), scientific = TRUE, trim = TRUE)
    },
    character(1)
  )
}

bmi_gradient_print_fold_change_extremes <- function(fold_df) {
  if (!nrow(fold_df)) {
    cat("\nBMIP gradient fold-change summary: no valid comparisons found.\n")
    return(invisible(NULL))
  }

  fmt_fc <- function(x) {
    format(round(as.numeric(x), 2), nsmall = 2, trim = TRUE)
  }
  fmt_sex <- function(x) {
    paste0(toupper(substr(as.character(x), 1L, 1L)), substr(as.character(x), 2L, nchar(as.character(x))))
  }

  min_row <- fold_df %>%
    dplyr::slice_min(.data$fold_change, n = 1L, with_ties = FALSE)
  max_row <- fold_df %>%
    dplyr::slice_max(.data$fold_change, n = 1L, with_ties = FALSE)

  cat("\nBMIP gradient fold-change summary\n")
  cat(
    "  Fold change = ratio of mean signed BMI-RNI gradients (numerator class / denominator class).\n"
  )
  cat("  Only structure-sex-comparison cells with positive numerator and denominator gradients are included.\n")
  cat(
    "  Comparisons: ",
    paste(unique(fold_df$comparison), collapse = "; "),
    "\n",
    sep = ""
  )
  cat(
    "  Structures: ",
    paste(sort(unique(fold_df$brain_region)), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(sprintf(
    paste0(
      "  Smallest fold change: %s (%s; %s, %s)\n",
      "    numerator = %s; denominator = %s\n"
    ),
    fmt_fc(min_row$fold_change[[1]]),
    min_row$comparison[[1]],
    min_row$brain_region[[1]],
    fmt_sex(min_row$sex[[1]]),
    bmi_gradient_format_gradient_sci(min_row$numerator[[1]]),
    bmi_gradient_format_gradient_sci(min_row$denominator[[1]])
  ))
  cat(sprintf(
    paste0(
      "  Largest fold change: %s (%s; %s, %s)\n",
      "    numerator = %s; denominator = %s\n"
    ),
    fmt_fc(max_row$fold_change[[1]]),
    max_row$comparison[[1]],
    max_row$brain_region[[1]],
    fmt_sex(max_row$sex[[1]]),
    bmi_gradient_format_gradient_sci(max_row$numerator[[1]]),
    bmi_gradient_format_gradient_sci(max_row$denominator[[1]])
  ))

  invisible(list(min = min_row, max = max_row, fold_changes = fold_df))
}

bmi_gradient_load_segment_definitions <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

bmi_gradient_segment_plan_from_defs <- function(defs) {
  if (is.null(defs) || !nrow(defs)) {
    return(list())
  }
  stats::setNames(
    lapply(seq_len(nrow(defs)), function(i) {
      c(defs$range_low[[i]], defs$range_high[[i]])
    }),
    defs$segment
  )
}

bmi_gradient_segment_plan_from_wide <- function(wide, value_col = "mean_signed") {
  pat <- paste0("_(", value_col, "|n_points|mean_abs|median_abs|max_abs)$")
  segs <- unique(sub(pat, "", grep(paste0("_", value_col, "$"), names(wide), value = TRUE)))
  segs <- segs[!segs %in% c("pheno", "xname", "designmat", "sex", "laterality", "modality", "roi", "file_path", "mapping_method")]
  if (!length(segs)) {
    return(list())
  }
  stats::setNames(vector("list", length(segs)), segs)
}

bmi_gradient_segment_display_label <- function(seg, segment_plan, cdc_labels = bmi_gradient_cdc_segment_labels()) {
  if (seg %in% names(cdc_labels)) {
    return(unname(cdc_labels[[seg]]))
  }
  gsub("_", " ", seg, fixed = TRUE)
}

bmi_gradient_scale_value <- function(x, exponent) {
  x <- as.numeric(x)
  dplyr::if_else(
    is.na(x),
    NA_character_,
    format(round(x / (10^exponent), 2), nsmall = 2, trim = TRUE)
  )
}

bmi_gradient_build_display_df <- function(
  wide,
  pheno,
  segments,
  segment_plan,
  rois,
  modality = "rni",
  value_col = "mean_signed",
  exponent = -5
) {
  mod <- wide %>%
    dplyr::filter(
      .data$pheno == .env$pheno,
      .data$modality == .env$modality,
      .data$roi %in% .env$rois
    )

  female <- mod %>% dplyr::filter(.data$sex == "female")
  male <- mod %>% dplyr::filter(.data$sex == "male")

  roi_order <- intersect(rois, unique(mod$roi))
  out <- tibble::tibble(
    brain_region = vapply(roi_order, bmi_gradient_clean_roi_name, character(1))
  )

  for (seg in segments) {
    src_col <- paste0(seg, "_", value_col)
    if (!src_col %in% names(mod)) {
      next
    }
    f_vals <- female[[src_col]][match(roi_order, female$roi)]
    m_vals <- male[[src_col]][match(roi_order, male$roi)]
    out[[paste0(seg, "__female")]] <- bmi_gradient_scale_value(f_vals, exponent)
    out[[paste0(seg, "__male")]] <- bmi_gradient_scale_value(m_vals, exponent)
  }

  out
}

bmi_gradient_build_display_df_sex <- function(
  wide,
  pheno,
  segments,
  segment_plan,
  rois,
  sex,
  modality = "rni",
  value_col = "mean_signed",
  exponent = -5
) {
  sex <- tolower(sex)
  mod <- wide %>%
    dplyr::filter(
      .data$pheno == .env$pheno,
      .data$modality == .env$modality,
      .data$roi %in% .env$rois,
      tolower(.data$sex) == .env$sex
    )

  roi_order <- intersect(rois, unique(mod$roi))
  out <- tibble::tibble(
    brain_region = vapply(roi_order, bmi_gradient_clean_roi_name, character(1))
  )

  for (seg in segments) {
    src_col <- paste0(seg, "_", value_col)
    if (!src_col %in% names(mod)) {
      next
    }
    vals <- mod[[src_col]][match(roi_order, mod$roi)]
    out[[seg]] <- bmi_gradient_scale_value(vals, exponent)
  }

  out
}

bmi_gradient_build_gt_table_sex <- function(
  display_df,
  segments,
  segment_plan,
  sex,
  exponent = -5,
  stub_label = "Brain region",
  table_style = bmi_table_style_config(),
  pheno = NULL,
  cdc_mapping = NULL
) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for gradient tables.", call. = FALSE)
  }

  seg_present <- segments[segments %in% names(display_df)]
  seg_present <- seg_present[seg_present != "brain_region"]
  if (!length(seg_present)) {
    stop("No gradient segment columns found in display data.", call. = FALSE)
  }

  use_cdc_labels <- !is.null(pheno) && all(grepl("^cdc_", seg_present))
  col_map <- c(brain_region = stub_label)
  for (seg in seg_present) {
    col_map[[seg]] <- if (isTRUE(use_cdc_labels)) {
      gt::md(bmi_gradient_cdc_spanner_label(seg, pheno, cdc_mapping = cdc_mapping, sex = sex))
    } else {
      lo <- segment_plan[[seg]][1]
      hi <- segment_plan[[seg]][2]
      if (!is.null(lo) && !is.null(hi) && !is.na(lo) && !is.na(hi)) {
        gt::md(paste0(
          bmi_gradient_segment_display_label(seg, segment_plan),
          "\n",
          bmi_gradient_format_axis_range(lo, hi)
        ))
      } else {
        bmi_gradient_segment_display_label(seg, segment_plan)
      }
    }
  }

  gt_tbl <- gt::gt(display_df, rowname_col = NULL)
  gt_tbl <- gt::cols_label(gt_tbl, .list = col_map)
  gt_tbl <- bmi_table_apply_gt_table_style(gt_tbl, table_style = table_style)
  gt_tbl <- gt::cols_align(gt_tbl, align = "left", columns = brain_region)
  gt_tbl <- gt::cols_align(gt_tbl, align = "center", columns = tidyselect::everything())
  gt_tbl <- gt::tab_style(
    gt_tbl,
    style = gt::cell_text(whitespace = "pre-line", align = "center"),
    locations = gt::cells_column_labels()
  )
  region_w <- as.numeric(table_style$brain_region_width_px)
  gt_tbl <- gt::cols_width(
    gt_tbl,
    stats::as.formula(paste0("brain_region ~ ", region_w))
  )

  gt_tbl
}

bmi_gradient_build_gt_table <- function(
  display_df,
  segments,
  segment_plan,
  exponent = -5,
  stub_label = "Brain region",
  table_style = bmi_table_style_config(),
  pheno = NULL,
  cdc_mapping = NULL
) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for gradient tables.", call. = FALSE)
  }

  seg_present <- segments[segments %in% sub("__female$", "", grep("__female$", names(display_df), value = TRUE))]
  if (!length(seg_present)) {
    stop("No gradient segment columns found in display data.", call. = FALSE)
  }

  col_map <- c(brain_region = stub_label)
  for (seg in seg_present) {
    col_map[[paste0(seg, "__female")]] <- "Female"
    col_map[[paste0(seg, "__male")]] <- "Male"
  }

  gt_tbl <- gt::gt(display_df, rowname_col = NULL)
  gt_tbl <- gt::cols_label(gt_tbl, .list = col_map)
  gt_tbl <- bmi_table_apply_gt_table_style(gt_tbl, table_style = table_style)
  gt_tbl <- gt::cols_align(gt_tbl, align = "left", columns = brain_region)
  gt_tbl <- gt::cols_align(gt_tbl, align = "center", columns = tidyselect::everything())
  gt_tbl <- gt::tab_style(
    gt_tbl,
    style = gt::cell_text(whitespace = "pre-line", align = "center"),
    locations = gt::cells_column_labels()
  )
  gt_tbl <- gt::tab_style(
    gt_tbl,
    style = gt::cell_text(whitespace = "pre-line", align = "center"),
    locations = gt::cells_column_spanners()
  )
  region_w <- as.numeric(table_style$brain_region_width_px)
  gt_tbl <- gt::cols_width(
    gt_tbl,
    stats::as.formula(paste0("brain_region ~ ", region_w))
  )

  use_cdc_spanners <- !is.null(pheno) && all(grepl("^cdc_", seg_present))

  for (seg in seg_present) {
    f_col <- paste0(seg, "__female")
    m_col <- paste0(seg, "__male")
    span_label <- if (isTRUE(use_cdc_spanners)) {
      gt::md(bmi_gradient_cdc_spanner_label(seg, pheno, cdc_mapping = cdc_mapping))
    } else {
      lo <- segment_plan[[seg]][1]
      hi <- segment_plan[[seg]][2]
      if (!is.null(lo) && !is.null(hi) && !is.na(lo) && !is.na(hi)) {
        gt::md(paste0(
          bmi_gradient_segment_display_label(seg, segment_plan),
          "\n",
          bmi_gradient_format_axis_range(lo, hi)
        ))
      } else {
        bmi_gradient_segment_display_label(seg, segment_plan)
      }
    }
    gt_tbl <- gt::tab_spanner(
      gt_tbl,
      label = span_label,
      columns = tidyselect::all_of(c(f_col, m_col))
    )
  }

  gt_tbl
}

bmi_gradient_footnote_config <- function(
  pheno,
  correction = "fdr",
  include_bmip_note = TRUE,
  cdc_mapping_method = NULL
) {
  pheno_note <- if (identical(pheno, "bmip")) {
    "Gradients were summarized on the CDC extended BMI percentile axis (bmip)."
  } else {
    "Gradients were summarized on the percent-from-median BMI axis (perc_median)."
  }

  cdc_note <- if (!is.null(cdc_mapping_method)) {
    paste0(
      "CDC weight-class segments for perc_median used the ",
      cdc_mapping_method,
      " bmip-to-perc_median mapping from the imaging cohort."
    )
  } else {
    NULL
  }

  entries <- list(
    list(
      letter = "a",
      text = paste(
        "Values are the mean first derivative of the BMI-RNI smooth",
        "(change in the partial BMI effect on RNI per BMI unit) within each BMI segment",
        "and are presented in scientific notation (\u00d710\u207b\u2075)."
      ),
      attach = "title"
    ),
    list(
      letter = "b",
      text = paste(pheno_note, "Multiple comparisons were corrected using", toupper(correction), "."),
      attach = NULL
    )
  )

  if (!is.null(cdc_note)) {
    entries[[length(entries) + 1]] <- list(
      letter = "c",
      text = cdc_note,
      attach = NULL
    )
  }

  abbreviations <- paste(
    "Abbreviations: BMI, body mass index; CDC, Centers for Disease Control and Prevention;",
    "RNI, restricted normalized isotropic signal fraction."
  )

  list(entries = entries, abbreviations = abbreviations)
}

bmi_gradient_cdc_part_labels <- function() {
  c(
    bmip = "A. CDC Extended BMI Percentile",
    perc_median = "B. Percent From Median BMI"
  )
}

bmi_gradient_cdc_part_labels_by_sex <- function() {
  c(
    female_bmip = "A. Female: CDC Extended BMI Percentile",
    male_bmip = "B. Male: CDC Extended BMI Percentile",
    female_perc_median = "C. Female: Percent From Median BMI",
    male_perc_median = "D. Male: Percent From Median BMI"
  )
}

bmi_gradient_cdc_table2_title <- function() {
  "Mean Rate of Change in the BMI-RNI Association by CDC Weight Class and Sex"
}

bmi_gradient_cdc_table2_by_sex_title <- function() {
  "Mean Rate of Change in the BMI-RNI Association by CDC Weight Class, Stratified by Sex"
}

bmi_gradient_cdc_table2_footnote_config <- function(correction = "fdr", cdc_method = "minmax") {
  correction_abbr <- switch(
    tolower(correction),
    fdr = "FDR, false discovery rate",
    bonferroni = "Bonferroni correction",
    all_voxels = "all-voxels threshold",
    toupper(correction)
  )

  list(
    entries = list(),
    general_notes = c(
      paste(
        "Values are the mean first derivative of the BMI-RNI smooth",
        "(change in the partial BMI effect on RNI per BMI unit) within each CDC weight class",
        "and are presented in scientific notation (\u00d710\u207b\u2075). CDC weight classes were",
        "defined as underweight (<5th percentile), healthy weight (5th to <85th percentile),",
        "overweight (85th to <95th percentile), and obesity (\u226595th percentile)."
      ),
      paste(
        "Positive values indicate that the BMI-associated contribution to RNI increases with",
        "BMI within that range; negative values indicate a decrease. Derivatives were averaged",
        "across", toupper(correction), "-significant voxels within each region of interest.",
        "Results are shown separately for females and males, derived from sex-stratified",
        "voxelwise linear mixed-effects models."
      ),
      paste0(
        "Part A summarizes results on the CDC extended BMI percentile axis. Part B summarizes",
        "results on the percent-from-median BMI axis, with CDC weight-class windows mapped",
        "from CDC extended BMI percentiles. Bracketed ranges in",
        "part B are sex-pooled perc_median windows (minimum lower bound to maximum upper",
        "bound across females and males)."
      )
    ),
    abbreviations = paste(
      "Abbreviations: BMI, body mass index; CDC, Centers for Disease Control and Prevention;",
      correction_abbr, "; RNI, restricted normalized isotropic signal fraction."
    )
  )
}

bmi_gradient_cdc_table2_by_sex_footnote_config <- function(correction = "fdr", cdc_method = "minmax") {
  correction_abbr <- switch(
    tolower(correction),
    fdr = "FDR, false discovery rate",
    bonferroni = "Bonferroni correction",
    all_voxels = "all-voxels threshold",
    toupper(correction)
  )

  list(
    entries = list(),
    general_notes = c(
      paste(
        "Values are the mean first derivative of the BMI-RNI smooth",
        "(change in the partial BMI effect on RNI per BMI unit) within each CDC weight class",
        "and are presented in scientific notation (\u00d710\u207b\u2075). CDC weight classes were",
        "defined as underweight (<5th percentile), healthy weight (5th to <85th percentile),",
        "overweight (85th to <95th percentile), and obesity (\u226595th percentile)."
      ),
      paste(
        "Positive values indicate that the BMI-associated contribution to RNI increases with",
        "BMI within that range; negative values indicate a decrease. Derivatives were averaged",
        "across", toupper(correction), "-significant voxels within each region of interest.",
        "Each part shows results for a single sex, derived from sex-stratified",
        "voxelwise linear mixed-effects models."
      ),
      paste0(
        "Parts A and C summarize results on the CDC extended BMI percentile axis for females",
        "and males, respectively. Parts B and D summarize results on the percent-from-median",
        "BMI axis, with CDC weight-class windows mapped from CDC extended BMI percentiles.",
        "Bracketed ranges in parts B and D are sex-specific perc_median windows."
      )
    ),
    abbreviations = paste(
      "Abbreviations: BMI, body mass index; CDC, Centers for Disease Control and Prevention;",
      correction_abbr, "; RNI, restricted normalized isotropic signal fraction."
    )
  )
}

bmi_gradient_default_title <- function(pheno, table_kind = c("main", "cdc")) {
  table_kind <- match.arg(table_kind)
  pheno_label <- if (identical(pheno, "bmip")) {
    "CDC Extended BMI Percentile"
  } else {
    "Percent From Median BMI"
  }
  if (table_kind == "cdc") {
    return(paste0(
      "Mean Rate of Change in the BMI-RNI Association by Brain Region and CDC Weight Class (",
      pheno_label, "), Stratified by Sex"
    ))
  }
  paste0(
    "Mean Rate of Change in the BMI-RNI Association by Brain Region and BMI Segment (",
    pheno_label, "), Stratified by Sex"
  )
}
