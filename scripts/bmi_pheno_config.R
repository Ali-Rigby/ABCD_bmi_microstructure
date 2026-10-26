# Phenotype profiles for unified BMI imaging pipeline.
# Set BMI_PHENO_VARIABLE before pheno-specific scripts: "bmip" or "perc_median".

bmi_pheno_profile <- function(pheno = Sys.getenv("BMI_PHENO_VARIABLE", unset = "")) {
  if (!nzchar(pheno)) {
    stop("Set BMI_PHENO_VARIABLE to 'bmip' or 'perc_median' before running pheno scripts.", call. = FALSE)
  }
  pheno <- trimws(pheno)

  profiles <- list(
    bmip = list(
      pheno_var = "bmip",
      vector_prefix = "bmi_percentile_vector",
      vector_file = "bmi_percentile_vector.txt",
      bf_col_prefix = "bf_bmi_",
      bf_table = "table_bf_bmi_perc.txt",
      bf_table_interp = "table_bf_bmi_perc_interp.txt",
      design_tag = "Sbmi",
      xname = "bmi",
      linear_term = "bmip",
      bf_grep = "bf_bmi",
      analysis_dir = "bmip",
      merged_datafile = "datafile6.0_with_bmi_age_pds_bfs.txt",
      label = "CDC extended BMI percentile"
    ),
    perc_median = list(
      pheno_var = "perc_median",
      vector_prefix = "perc_median_vector",
      vector_file = "perc_median_vector.txt",
      vector_file = "perc_median_vector.txt",
      bf_col_prefix = "bf_perc_median_",
      bf_table = "table_bf_perc_median.txt",
      bf_table_interp = "table_bf_perc_median_interp.txt",
      design_tag = "Sperc",
      xname = "perc_median",
      linear_term = "perc_median",
      bf_grep = "bf_perc_median",
      analysis_dir = "perc_median",
      merged_datafile = "datafile6.0_with_perc_median_age_pds_bfs.txt",
      label = "% from median (age/sex)"
    )
  )

  if (!pheno %in% names(profiles)) {
    stop("BMI_PHENO_VARIABLE must be 'bmip' or 'perc_median'.", call. = FALSE)
  }
  profiles[[pheno]]
}
