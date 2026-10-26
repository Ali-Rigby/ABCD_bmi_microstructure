# Detect >100-day discrepancies between anthropometry and admin ages per observation.

suppressMessages({
  library(arrow)
  library(dplyr)
  library(readr)
})

imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
PATHS <- bmi_shared_paths()
tab <- PATHS$abcd_tabulated

path_anthr <- file.path(tab, "ph_y_anthr.parquet")
path_admin <- file.path(tab, "mr_y_adm__qtn.parquet")
path_dyn <- file.path(tab, "ab_g_dyn.parquet")

id_cols <- c("participant_id", "session_id")
col_age_anthr <- "ph_y_anthr_age"
col_age_admin <- "mr_y_adm__qtn__pre__r01_age"
col_age_dyn   <- "ab_g_dyn__visit_age"

read_needed <- function(path, cols) {
  ds <- tryCatch(arrow::open_dataset(path), error = function(e) NULL)
  if (!is.null(ds)) {
    return(ds %>% dplyr::select(dplyr::all_of(cols)) %>% dplyr::collect())
  }
  df <- arrow::read_parquet(path, as_data_frame = TRUE)
  dplyr::select(df, dplyr::all_of(intersect(cols, names(df))))
}

message("Reading anthropometry parquet…")
anthr <- read_needed(path_anthr, c(id_cols, col_age_anthr))
message("Reading admin parquet…")
admin <- read_needed(path_admin, c(id_cols, col_age_admin))
message("Reading AB dyn parquet…")
dyn <- read_needed(path_dyn, c(id_cols, col_age_dyn))

stopifnot(col_age_anthr %in% names(anthr))
stopifnot(col_age_admin %in% names(admin))
stopifnot(all(id_cols %in% names(anthr)))
stopifnot(all(id_cols %in% names(admin)))
stopifnot(col_age_dyn %in% names(dyn))
stopifnot(all(id_cols %in% names(dyn)))

anthr <- anthr %>% mutate(!!col_age_anthr := suppressWarnings(as.numeric(.data[[col_age_anthr]])))
admin <- admin %>% mutate(!!col_age_admin := suppressWarnings(as.numeric(.data[[col_age_admin]])))
dyn <- dyn %>% mutate(!!col_age_dyn := suppressWarnings(as.numeric(.data[[col_age_dyn]])))

dat <- anthr %>% inner_join(admin, by = id_cols) %>% left_join(dyn, by = id_cols)

days_per_year <- 365.25
threshold_days <- 100
threshold_years <- threshold_days / days_per_year

dat <- dat %>% mutate(
  age_diff_years = abs(.data[[col_age_anthr]] - .data[[col_age_admin]]),
  age_diff_days  = age_diff_years * days_per_year,
  discrepant_100d = !is.na(age_diff_days) & is.finite(age_diff_days) & age_diff_days > threshold_days,
  anthr_vs_dyn_days = ifelse(!is.na(.data[[col_age_dyn]]),
                             abs(.data[[col_age_anthr]] - .data[[col_age_dyn]]) * days_per_year, NA_real_),
  admin_vs_dyn_days = ifelse(!is.na(.data[[col_age_dyn]]),
                             abs(.data[[col_age_admin]] - .data[[col_age_dyn]]) * days_per_year, NA_real_)
)

total <- nrow(dat)
n_discrepant <- sum(dat$discrepant_100d, na.rm = TRUE)
pct <- if (total > 0) round(100 * n_discrepant / total, 2) else NA_real_
message(sprintf("Observations joined: %d", total))
message(sprintf(">100-day discrepancies: %d (%.2f%%)", n_discrepant, pct))

out_dir <- PATHS$shared_datafiles
bmi_ensure_dir(out_dir)
out_path <- file.path(out_dir, "age_scan_anthro_discrepancy_report.csv")
dat_out <- dat %>%
  mutate(threshold_days = threshold_days, threshold_years = threshold_years) %>%
  arrange(desc(discrepant_100d), participant_id, session_id)
readr::write_csv(dat_out, out_path)
message(sprintf("Report written: %s", out_path))

discrepant_out <- dat_out %>%
  filter(discrepant_100d) %>%
  select(all_of(id_cols), !!col_age_anthr, !!col_age_admin, !!col_age_dyn,
         age_diff_days, anthr_vs_dyn_days, admin_vs_dyn_days)

out_path_disc <- file.path(out_dir, "age_scan_anthro_discrepancy_discrepant_only.csv")
readr::write_csv(discrepant_out, out_path_disc)
message(sprintf("Discrepant-only report written: %s", out_path_disc))
print(discrepant_out)
