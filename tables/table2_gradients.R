################################
# tables/table2_gradients.R
# Manuscript Table 2 — CDC weight-class gradient tables.
#
# Optional env:
  # BMI_IMAGING_GRADIENT_CORRECTION=fdr|bonferroni|all_voxels
  # BMI_IMAGING_GRADIENT_PHENOS=bmip,perc_median
  # BMI_IMAGING_GRADIENT_CDC_METHOD=minmax
  # BMI_IMAGING_GRADIENT_TABLE_NUMBER=2
  # BMI_IMAGING_TABLE_FONT_FAMILY=Arial
  # BMI_IMAGING_TABLE_BODY_FONT_SIZE=10
  # BMI_IMAGING_TABLE_LABEL_FONT_SIZE=9
  # BMI_IMAGING_TABLE_TITLE_FONT_SIZE=12
  # BMI_IMAGING_TABLE_FOOTNOTE_FONT_SIZE=8
  # BMI_IMAGING_TABLE_PART_PAGE_BREAK=0   set 1 to force part B onto a new page
  # BMI_IMAGING_GRADIENT_BY_SEX=1         set 0 to skip sex-stratified Table 2 outputs
#
################################

required_pkgs <- c("tidyverse", "gt")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(
      "Package '", p, "' is required. Example:\n",
      "  conda install -c conda-forge r-tidyverse r-gt r-officer",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_demos_data.R"))
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "utilities", "r_utils", "bmi_gradient.R"))

ceph_root <- Sys.getenv("BMI_IMAGING_CEPH_ROOT", bmi_default_ceph_root())
results_run <- Sys.getenv("BMI_IMAGING_RESULTS_RUN", unset = "")
if (!nzchar(results_run)) {
  results_run <- paste0("results_", bmi_run_tag())
}

correction <- tolower(Sys.getenv("BMI_IMAGING_GRADIENT_CORRECTION", unset = "fdr"))
phenos <- strsplit(Sys.getenv("BMI_IMAGING_GRADIENT_PHENOS", "bmip,perc_median"), ",", fixed = TRUE)[[1]]
phenos <- tolower(trimws(phenos))
phenos <- intersect(c("bmip", "perc_median"), phenos)
cdc_method <- Sys.getenv("BMI_IMAGING_GRADIENT_CDC_METHOD", unset = "minmax")
table_number <- Sys.getenv("BMI_IMAGING_GRADIENT_TABLE_NUMBER", unset = "2")
exponent <- as.integer(Sys.getenv("BMI_IMAGING_GRADIENT_EXPONENT", "-5"))
table_style <- bmi_table_style_config()
gradient_table_style <- table_style
selected_rois <- bmi_gradient_default_rois()
part_labels <- bmi_gradient_cdc_part_labels()
part_labels_by_sex <- bmi_gradient_cdc_part_labels_by_sex()
build_by_sex <- tolower(Sys.getenv("BMI_IMAGING_GRADIENT_BY_SEX", "1")) %in%
  c("1", "true", "yes", "on")

out_dir <- file.path(ceph_root, "shared", "gradients")
bmi_ensure_dir(out_dir)

pheno_meta <- list(
  bmip = list(analysis_dir = "bmip"),
  perc_median = list(analysis_dir = "perc_median")
)

build_cdc_table_part <- function(wide, pheno, part_id, part_label, grad_dir, table_style = gradient_table_style) {
  cdc_segments <- names(bmi_gradient_cdc_segment_labels())
  cdc_plan <- bmi_gradient_cdc_segment_ranges()
  cdc_mapping <- if (identical(pheno, "perc_median")) {
    bmi_gradient_load_cdc_mapping(grad_dir, cdc_method, correction)
  } else {
    NULL
  }

  display_df <- bmi_gradient_build_display_df(
    wide = wide,
    pheno = pheno,
    segments = cdc_segments,
    segment_plan = cdc_plan,
    rois = selected_rois,
    exponent = exponent
  )
  if (!nrow(display_df)) {
    return(NULL)
  }

  list(
    part_id = part_id,
    part_label = part_label,
    pheno = pheno,
    gt_tbl = bmi_gradient_build_gt_table(
      display_df = display_df,
      segments = cdc_segments,
      segment_plan = cdc_plan,
      exponent = exponent,
      table_style = table_style,
      pheno = pheno,
      cdc_mapping = cdc_mapping
    ),
    display_df = display_df
  )
}

build_cdc_table_part_sex <- function(
  wide,
  pheno,
  sex,
  part_id,
  part_label,
  grad_dir,
  table_style = gradient_table_style
) {
  cdc_segments <- names(bmi_gradient_cdc_segment_labels())
  cdc_plan <- bmi_gradient_cdc_segment_ranges()
  cdc_mapping <- if (identical(pheno, "perc_median")) {
    bmi_gradient_load_cdc_mapping(grad_dir, cdc_method, correction)
  } else {
    NULL
  }

  display_df <- bmi_gradient_build_display_df_sex(
    wide = wide,
    pheno = pheno,
    segments = cdc_segments,
    segment_plan = cdc_plan,
    rois = selected_rois,
    sex = sex,
    exponent = exponent
  )
  if (!nrow(display_df)) {
    return(NULL)
  }

  list(
    part_id = part_id,
    part_label = part_label,
    pheno = pheno,
    sex = sex,
    gt_tbl = bmi_gradient_build_gt_table_sex(
      display_df = display_df,
      segments = cdc_segments,
      segment_plan = cdc_plan,
      sex = sex,
      exponent = exponent,
      table_style = table_style,
      pheno = pheno,
      cdc_mapping = cdc_mapping
    ),
    display_df = display_df
  )
}

cat("\nGradient Table 2 (CDC weight classes)\n")
cat("  results_run:", results_run, "\n")
cat("  correction:  ", correction, "\n")
cat("  cdc_method:  ", cdc_method, "\n")
cat("  by_sex:      ", if (build_by_sex) "yes" else "no", "\n")
cat("  output:      ", out_dir, "\n\n")

table_parts <- list()
table_parts_by_sex <- list()
wide_by_pheno <- list()

for (pheno in phenos) {
  if (!pheno %in% names(pheno_meta)) {
    warning("Unknown pheno profile: ", pheno, call. = FALSE)
    next
  }
  if (!pheno %in% names(part_labels)) {
    warning("No part label configured for pheno: ", pheno, call. = FALSE)
    next
  }

  grad_dir <- bmi_gradient_results_dir(ceph_root, pheno_meta[[pheno]]$analysis_dir, results_run)
  tag <- paste0(pheno, "_", correction)
  wide_cdc_path <- file.path(grad_dir, paste0("gradient_summary_wide_cdc_", cdc_method, "_", tag, ".csv"))

  if (!file.exists(wide_cdc_path)) {
    warning(
      "Missing CDC wide gradient CSV: ", wide_cdc_path,
      "\n  Run s11_compute_gradient_differences.ipynb export cell (SAVE_CSV = TRUE) for this results run.",
      call. = FALSE
    )
    next
  }

  wide_cdc <- bmi_gradient_load_wide(wide_cdc_path)
  wide_by_pheno[[pheno]] <- wide_cdc
  bmi_gradient_warn_missing_rois(wide_cdc, selected_rois, context = pheno)

  if (identical(pheno, "bmip")) {
    fold_changes_bmip <- bmi_gradient_collect_bmip_fold_changes(
      wide = wide_cdc,
      rois = selected_rois
    )
    bmi_gradient_print_fold_change_extremes(fold_changes_bmip)
  }

  part <- build_cdc_table_part(
    wide = wide_cdc,
    pheno = pheno,
    part_id = pheno,
    part_label = part_labels[[pheno]],
    grad_dir = grad_dir
  )
  if (is.null(part)) {
    warning("No rows for CDC gradient table part (", pheno, "); skipping.", call. = FALSE)
    next
  }

  cat("Part:", part$part_label, "| pheno:", pheno, "\n")
  table_parts[[length(table_parts) + 1L]] <- part
}

if (isTRUE(build_by_sex)) {
  sex_pheno_plan <- list(
    list(sex = "female", pheno = "bmip"),
    list(sex = "male", pheno = "bmip"),
    list(sex = "female", pheno = "perc_median"),
    list(sex = "male", pheno = "perc_median")
  )

  for (plan in sex_pheno_plan) {
    sex <- plan$sex
    pheno <- plan$pheno
    if (!pheno %in% phenos) {
      next
    }
    if (!pheno %in% names(wide_by_pheno)) {
      next
    }

    part_id <- paste0(sex, "_", pheno)
    part_label <- part_labels_by_sex[[part_id]]
    if (is.null(part_label)) {
      warning("No sex-stratified part label configured for: ", part_id, call. = FALSE)
      next
    }

    grad_dir <- bmi_gradient_results_dir(ceph_root, pheno_meta[[pheno]]$analysis_dir, results_run)
    part <- build_cdc_table_part_sex(
      wide = wide_by_pheno[[pheno]],
      pheno = pheno,
      sex = sex,
      part_id = part_id,
      part_label = part_label,
      grad_dir = grad_dir
    )
    if (is.null(part)) {
      warning("No rows for sex-stratified CDC gradient table part (", part_id, "); skipping.", call. = FALSE)
      next
    }

    cat("Sex-stratified part:", part$part_label, "| sex:", sex, "| pheno:", pheno, "\n")
    table_parts_by_sex[[length(table_parts_by_sex) + 1L]] <- part
  }
}

if (!length(table_parts)) {
  stop(
    "No CDC gradient table parts were built. Expected s11 CSV exports such as ",
    "gradient_summary_wide_cdc_", cdc_method, "_bmip_", correction, ".csv under each ",
    "pheno's results/", results_run, "/gradient_analysis/. ",
    "Re-run s11_compute_gradient_differences.ipynb with SAVE_CSV = TRUE.",
    call. = FALSE
  )
}

file_prefix <- paste0(
  "table", table_number, "_gradient_cdc_",
  paste(vapply(table_parts, function(x) x$pheno, character(1)), collapse = "_"),
  "_", cdc_method, "_", correction
)

bmi_table_save_multipart_gt_table(
  table_number = table_number,
  table_title = bmi_gradient_cdc_table2_title(),
  footnote_config = bmi_gradient_cdc_table2_footnote_config(correction, cdc_method),
  table_parts = table_parts,
  out_dir = out_dir,
  file_prefix = file_prefix,
  table_style = table_style
)

if (length(table_parts_by_sex)) {
  file_prefix_by_sex <- paste0(
    "table", table_number, "_gradient_cdc_by_sex_",
    paste(phenos, collapse = "_"),
    "_", cdc_method, "_", correction
  )

  bmi_table_save_multipart_gt_table(
    table_number = table_number,
    table_title = bmi_gradient_cdc_table2_by_sex_title(),
    footnote_config = bmi_gradient_cdc_table2_by_sex_footnote_config(correction, cdc_method),
    table_parts = table_parts_by_sex,
    out_dir = out_dir,
    file_prefix = file_prefix_by_sex,
    table_style = table_style
  )
}

cat(
  "\nDone. Wrote Table ", table_number, " (", length(table_parts), " combined-sex part(s)",
  if (length(table_parts_by_sex)) {
    paste0(" + ", length(table_parts_by_sex), " sex-stratified part(s)")
  } else {
    ""
  },
  ") to ", out_dir, ".\n",
  sep = ""
)
