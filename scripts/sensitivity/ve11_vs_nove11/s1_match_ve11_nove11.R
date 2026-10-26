################################
# 1:1 matching: smaller arm -> larger arm on perc_median, then age.
# VE11 vs noVE11 on full-sample 6.0 imaging base (all waves).
#
# Prerequisite: standard s1 on $CEPH/shared/datafiles.
# Source env/run_ve11_vs_noVE11_6_0_shared.env before running.
################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

imaging_repo_root <- Sys.getenv(
  "BMI_IMAGING_REPO_ROOT",
  "/home/arigby/projects/ABCD_bmi_microstructure"
)
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
core <- bmi_project_paths_core()

caliper_perc <- as.numeric(Sys.getenv("BMI_MATCH_PERC_MEDIAN_CALIPER", unset = "10"))
caliper_age <- as.numeric(Sys.getenv("BMI_MATCH_AGE_CALIPER", unset = "0.25"))
if (!is.finite(caliper_perc) || caliper_perc <= 0) {
  stop("Set BMI_MATCH_PERC_MEDIAN_CALIPER to a positive number.", call. = FALSE)
}
if (!is.finite(caliper_age) || caliper_age <= 0) {
  stop("Set BMI_MATCH_AGE_CALIPER to a positive number.", call. = FALSE)
}

ve11_pat <- trimws(Sys.getenv("BMI_VE11_SOFTWARE_PATTERN", unset = "syngo[ _]MR[ _]E11"))
if (!nzchar(ve11_pat)) {
  stop("Set BMI_VE11_SOFTWARE_PATTERN to a non-empty regex.", call. = FALSE)
}

input_base <- trimws(Sys.getenv(
  "BMI_VE11_NOVE11_MATCH_INPUT_BASE",
  unset = "datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_dmri_incl.txt"
))
infile <- file.path(core$ceph_project, "shared", "datafiles", input_base)
bmi_assert_file(infile, "Full-sample imaging base (run standard s1 first)")

analysis_root <- trimws(Sys.getenv("BMI_VE11_NOVE11_ANALYSIS_ROOT", unset = "ve11_vs_noVE11_6_0"))
analysis_base <- file.path(core$ceph_project, analysis_root)
shared_out <- file.path(analysis_base, "shared", "datafiles")
matching_dir <- file.path(analysis_base, "matching")
qc_dir <- file.path(matching_dir, "qc")
bmi_ensure_dir(shared_out)
bmi_ensure_dir(qc_dir)

infile_df <- read.delim(infile, stringsAsFactors = FALSE)
bmi_assert_cols(
  infile_df,
  c(
    "participant_id", "session_id", "perc_median", "bmip", "ab_g_dyn__visit_age",
    "ab_g_dyn__design_mr__serial_software", "ab_g_dyn__design_mr__serial",
    "ab_g_dyn__design_mr__software", "ab_g_stc__cohort_sex"
  ),
  "Imaging base datafile"
)

pool <- infile_df %>%
  mutate(
    bmi_ve11_arm = dplyr::case_when(
      is.na(ab_g_dyn__design_mr__serial_software) ~ NA_character_,
      grepl(ve11_pat, as.character(ab_g_dyn__design_mr__serial_software), perl = TRUE) ~ "VE11",
      TRUE ~ "noVE11"
    ),
    visit_key = paste(participant_id, session_id, sep = "|")
  ) %>%
  filter(!is.na(bmi_ve11_arm), !is.na(perc_median), !is.na(ab_g_dyn__visit_age))

arm_n <- table(pool$bmi_ve11_arm)
n_ve11 <- as.integer(if ("VE11" %in% names(arm_n)) arm_n[["VE11"]] else 0L)
n_nove11 <- as.integer(if ("noVE11" %in% names(arm_n)) arm_n[["noVE11"]] else 0L)
cat(sprintf(
  "Matching pool: VE11=%d, noVE11=%d | calipers: perc_median±%.2f, age±%.2f yr\n",
  n_ve11, n_nove11, caliper_perc, caliper_age
))
cat("VE11 software pattern:", ve11_pat, "\n")

if (n_ve11 == 0L || n_nove11 == 0L) {
  stop("Need both VE11 and noVE11 visits in the matching pool.", call. = FALSE)
}

# Match smaller arm into larger (same spirit as XA30 -> VE11)
if (n_nove11 <= n_ve11) {
  donor_arm <- "noVE11"
  pool_arm <- "VE11"
} else {
  donor_arm <- "VE11"
  pool_arm <- "noVE11"
}
cat(sprintf("Match direction: %s (donor) -> %s (pool)\n", donor_arm, pool_arm))

donors <- pool %>% filter(bmi_ve11_arm == donor_arm) %>% arrange(participant_id, session_id)
targets <- pool %>% filter(bmi_ve11_arm == pool_arm) %>% arrange(participant_id, session_id)

used_keys <- character(0)
pairs <- vector("list", nrow(donors))
n_unmatched <- 0L

for (i in seq_len(nrow(donors))) {
  row_d <- donors[i, ]
  candidates <- targets %>%
    filter(!visit_key %in% used_keys) %>%
    filter(
      abs(perc_median - row_d$perc_median) <= caliper_perc,
      abs(ab_g_dyn__visit_age - row_d$ab_g_dyn__visit_age) <= caliper_age
    ) %>%
    mutate(
      match_dist = abs(perc_median - row_d$perc_median) +
        0.1 * abs(ab_g_dyn__visit_age - row_d$ab_g_dyn__visit_age)
    ) %>%
    arrange(match_dist, participant_id, session_id)

  if (!nrow(candidates)) {
    n_unmatched <- n_unmatched + 1L
    next
  }

  best <- candidates[1, ]
  used_keys <- c(used_keys, best$visit_key)
  pairs[[i]] <- bind_rows(
    row_d %>% mutate(match_role = donor_arm),
    best %>% mutate(match_role = pool_arm)
  )
}

pairs <- pairs[!vapply(pairs, is.null, logical(1))]
if (!length(pairs)) {
  stop("No matched pairs found. Widen calipers or inspect input base.", call. = FALSE)
}

matched_pairs <- bind_rows(pairs) %>%
  group_by(participant_id, session_id, match_role) %>%
  slice(1) %>%
  ungroup()

n_pairs <- length(pairs)
cat(sprintf(
  "Matched pairs: %d / %d %s visits (unmatched %s=%d)\n",
  n_pairs, nrow(donors), donor_arm, donor_arm, n_unmatched
))

pair_ids <- matched_pairs %>%
  mutate(pair_id = rep(seq_len(n_pairs), each = 2)) %>%
  select(
    pair_id, match_role, participant_id, session_id, bmi_ve11_arm,
    perc_median, bmip, ab_g_dyn__visit_age, ab_g_stc__cohort_sex,
    ab_g_dyn__design_mr__serial, ab_g_dyn__design_mr__software,
    ab_g_dyn__design_mr__serial_software
  )

write.table(
  pair_ids, file.path(matching_dir, "matched_pairs.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

matched_visits <- matched_pairs %>%
  distinct(participant_id, session_id, bmi_ve11_arm)

write.table(
  matched_visits %>% filter(bmi_ve11_arm == "VE11") %>% select(participant_id, session_id),
  file.path(matching_dir, "matched_ve11_visits.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  matched_visits %>% filter(bmi_ve11_arm == "noVE11") %>% select(participant_id, session_id),
  file.path(matching_dir, "matched_nove11_visits.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Join arm label onto full imaging rows for matched visits
matched_combined <- infile_df %>%
  inner_join(
    matched_visits,
    by = c("participant_id", "session_id")
  )

out_combined <- file.path(shared_out, "datafile_matched_combined.txt")
write.table(matched_combined, file = out_combined, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote matched combined table:", out_combined, "(", nrow(matched_combined), "rows )\n")

balance <- matched_combined %>%
  group_by(bmi_ve11_arm) %>%
  summarise(
    n = n(),
    age_mean = mean(ab_g_dyn__visit_age, na.rm = TRUE),
    age_sd = sd(ab_g_dyn__visit_age, na.rm = TRUE),
    perc_mean = mean(perc_median, na.rm = TRUE),
    perc_sd = sd(perc_median, na.rm = TRUE),
    bmip_mean = mean(bmip, na.rm = TRUE),
    bmip_sd = sd(bmip, na.rm = TRUE),
    pct_female = mean(ab_g_stc__cohort_sex == "F", na.rm = TRUE) * 100,
    n_sessions = dplyr::n_distinct(session_id),
    .groups = "drop"
  )
write.csv(balance, file.path(qc_dir, "post_match_balance_summary.csv"), row.names = FALSE)
print(balance)

session_counts <- matched_combined %>%
  count(bmi_ve11_arm, session_id, name = "n") %>%
  arrange(bmi_ve11_arm, session_id)
write.csv(session_counts, file.path(qc_dir, "post_match_session_counts.csv"), row.names = FALSE)

match_meta <- data.frame(
  donor_arm = donor_arm,
  pool_arm = pool_arm,
  n_pairs = n_pairs,
  n_unmatched_donor = n_unmatched,
  caliper_perc_median = caliper_perc,
  caliper_age = caliper_age,
  ve11_pattern = ve11_pat,
  input_file = infile,
  stringsAsFactors = FALSE
)
write.csv(match_meta, file.path(qc_dir, "match_meta.csv"), row.names = FALSE)

cat("Done.\n")
