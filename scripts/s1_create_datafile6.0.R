################################

# This example script describes how to combine variables from ABCD instruments 
# and save as a .txt file. It can be used as input to the script 
# run_makeDesign_example.R which describes how to create design
# matrices for the FEMA package. Users should save a local copy and edit paths 
# appropriately. 

################################
# The following R packages need to be loaded

for (p in c("tidyverse", "psych", "plyr", "dplyr", "pracma", "arrow", "splines", "data.table")) {
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

# Central paths + validation helpers (see scripts/config_project.R)
# Unified pipeline: source env/run_shared.csh before running s1.
imaging_repo_root <- Sys.getenv("BMI_IMAGING_REPO_ROOT", "/home/arigby/projects/ABCD_bmi_microstructure")
source(file.path(normalizePath(imaging_repo_root, mustWork = FALSE), "scripts", "config_project.R"))
PATHS <- bmi_shared_paths()

# Load CDC reference data and cdcanthro function from local files
bmi_assert_file(PATHS$cdc_ref_rda, "CDC reference .rda")
load(PATHS$cdc_ref_rda)  # Loads cdc__ref__data
bmi_assert_file(PATHS$cdcanthro_r, "cdcanthro.R")
source(PATHS$cdcanthro_r)

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- PATHS$abcd_tabulated
bmi_assert_dir(inpath, "ABCD tabulated directory")

# Define the full path to the output RDS file 
outpath <- PATHS$shared_datafiles
bmi_ensure_dir(outpath)

fname_pheno_qc <- "datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_removed"
fname_imaging <- "datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_dmri_incl"
outmatfile_pheno_qc <- paste0(outpath, "/", fname_pheno_qc)
outmatfile_imaging <- paste0(outpath, "/", fname_imaging)

# Define the path to tge cmig_utils/r directory, R needs to be able to 
# parse functions from this directory
funcpath <- PATHS$cmig_tools_r
bmi_assert_dir(funcpath, "cmig_tools R utilities directory")

################################
# Data from DEAP
agefile <- paste0(inpath, '/ab_g_dyn.parquet')
generalfile <- paste0(inpath, '/ab_g_stc.parquet')
demosfile <- paste0(inpath, '/ab_p_demo.parquet')
mriinfofile <- paste0(inpath, '/mr_y_adm__info.parquet')
mriqcfile <- paste0(inpath, '/mr_y_qc__incl.parquet')
mrithkfile <- paste0(inpath, '/mr_y_smri__thk__dsk.parquet')
mriareafile <- paste0(inpath, '/mr_y_smri__area__dsk.parquet')
mrivolfile <- paste0(inpath, '/mr_y_smri__vol__aseg.parquet')
physfile <- paste0(inpath, '/ph_y_anthr.parquet')
pds_parent_file <- paste0(inpath, '/ph_p_pds.parquet')
pds_youth_file <- paste0(inpath, '/ph_y_pds.parquet')

################################
# Load DYN file and extract age + dyn-sourced covariates (scanner, edu, income6)
dyn <- arrow::read_parquet(agefile)
agevars <- c('participant_id', 'session_id', 'ab_g_dyn__visit_age')
age <- dyn[, agevars]

dynvars <- c('participant_id','session_id',
             'ab_g_dyn__design_mr__serial',
             'ab_g_dyn__design_mr__software',
             'ab_g_dyn__cohort_edu__cgs',
             'ab_g_dyn__cohort_income__hhold__6lvl')
dyn_keep <- dyn[, dynvars]

outmat <- join(age, dyn_keep, by = c('participant_id','session_id'))

# Combine scanner serial and software into one variable
outmat <- outmat %>%
  mutate(ab_g_dyn__design_mr__serial_software = ifelse(
    is.na(ab_g_dyn__design_mr__serial) | is.na(ab_g_dyn__design_mr__software),
    NA_character_,
    paste0(ab_g_dyn__design_mr__serial, '_', ab_g_dyn__design_mr__software)
  ))

################################
# Load GENERAL file and extract fam id, sex, genetic PCs 
general <- arrow::read_parquet(generalfile)
genpcvars <- grep('ab_g_stc__gen_pc__', names(general), value = TRUE)
genpcvars <- sort(genpcvars)[1:10]
generalvars <- c('participant_id', 'ab_g_stc__design_id__fam', 'ab_g_stc__cohort_sex', genpcvars)
general <- general[, generalvars]

# Convert 1 → 'M' and 2 → 'F' in ab_g_stc__cohort_sex
general <- general %>%
  mutate(ab_g_stc__cohort_sex = recode(ab_g_stc__cohort_sex,
                                       `1` = "M",
                                       `2` = "F"))

outmat <- join(outmat, general, by = c('participant_id'))

################################
# Load MRI QC file and extract include vars
mriqc <- arrow::read_parquet(mriqcfile)
inclvar <- grep('mr_y_qc__incl', names(mriqc), value=TRUE)
mriqcvar <- c('participant_id', 'session_id', inclvar)
mriqc <- mriqc[, mriqcvar]

outmat <- join(outmat, mriqc, by = c('participant_id', 'session_id'))

# ################################
# # Load MRI thickness file and extract 
# thk <- arrow::read_parquet(mrithkfile)
# thkvars <- c('participant_id', 'session_id', 'mr_y_smri__thk__dsk_mean')
# thk <- thk[, thkvars]

# outmat <- join(outmat, thk, by = c('participant_id', 'session_id'))

# ################################
# # Load MRI area file and extract
# area <- arrow::read_parquet(mriareafile)
# areavars <- c('participant_id', 'session_id', 'mr_y_smri__area__dsk_sum')
# area <- area[, areavars]

# outmat <- join(outmat, area, by = c('participant_id', 'session_id'))

# ################################
# # Load MRI volume file and extract
# vol <- arrow::read_parquet(mrivolfile)
# volvars <- c('participant_id', 'session_id', 'mr_y_smri__vol__aseg__icv_sum')
# vol <- vol[, volvars]

# outmat <- join(outmat, vol, by = c('participant_id', 'session_id'))

################################
# Load pds file and compute average of youth and parent scores
pds_parent <- arrow::read_parquet(pds_parent_file)
pds_youth <- arrow::read_parquet(pds_youth_file)
pds_parent <- pds_parent[, c('participant_id', 'session_id', 'ph_p_pds__f_categ', 'ph_p_pds__m_categ')]
pds_youth <- pds_youth[, c('participant_id', 'session_id', 'ph_y_pds__f_categ', 'ph_y_pds__m_categ')]

# Merge PDS data to check for regression in scores across timepoints
pds_all <- full_join(pds_youth, pds_parent, by = c("participant_id", "session_id"))

# Add sex information for stratified analysis
pds_all <- pds_all %>%
  dplyr::left_join(general[, c('participant_id', 'ab_g_stc__cohort_sex')], by = 'participant_id')

# Check for PDS regression before computing averages
cat("\n--- PDS Regression Check ---\n")

# Initialize data frames to store results
regression_summary <- data.frame()
regression_details <- data.frame()
all_trajectories <- data.frame()

# Define instruments to analyze
instruments <- list(
  list(name = "ph_p_pds__f_categ + ph_p_pds__m_categ", 
       female_var = "ph_p_pds__f_categ", 
       male_var = "ph_p_pds__m_categ"),
  list(name = "ph_y_pds__f_categ + ph_y_pds__m_categ", 
       female_var = "ph_y_pds__f_categ", 
       male_var = "ph_y_pds__m_categ")
)

for (instrument in instruments) {
  female_var <- instrument$female_var
  male_var <- instrument$male_var
  
  # Get female data
  pds_temp_f <- pds_all %>%
    dplyr::filter(ab_g_stc__cohort_sex == "F") %>%
    dplyr::select(participant_id, session_id, ab_g_stc__cohort_sex, all_of(female_var)) %>%
    dplyr::filter(!is.na(.data[[female_var]])) %>%
    dplyr::mutate(pds_numeric = as.numeric(.data[[female_var]])) %>%
    dplyr::arrange(participant_id, session_id)
  
  # Get male data
  pds_temp_m <- pds_all %>%
    dplyr::filter(ab_g_stc__cohort_sex == "M") %>%
    dplyr::select(participant_id, session_id, ab_g_stc__cohort_sex, all_of(male_var)) %>%
    dplyr::filter(!is.na(.data[[male_var]])) %>%
    dplyr::mutate(pds_numeric = as.numeric(.data[[male_var]])) %>%
    dplyr::arrange(participant_id, session_id)
  
  # Combine for full dataset
  pds_full_dataset <- bind_rows(pds_temp_f, pds_temp_m)
  
  # Print combined instrument total
  cat(sprintf("\n%s:\n", instrument$name))
  cat("  Combined instrument total (M+F):\n")
  
  if (nrow(pds_full_dataset) == 0) {
    cat("    No data available\n")
  } else {
    n_total_obs <- nrow(pds_full_dataset)
    
    # Get total observations by session (for denominator in session breakdown)
    session_totals_combined <- pds_full_dataset %>%
      dplyr::group_by(session_id) %>%
      dplyr::summarise(total_obs = dplyr::n(), .groups = 'drop')
    
    # Group by participant and identify regressions
    pds_check <- pds_full_dataset %>%
      dplyr::group_by(participant_id) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::mutate(
        is_regression = pds_numeric < dplyr::lag(pds_numeric)
      ) %>%
      dplyr::ungroup()
    
    # Get participant-level summary
    pds_participant_summary <- pds_check %>%
      dplyr::group_by(participant_id) %>%
      dplyr::summarise(
        has_any_regression = any(is_regression, na.rm = TRUE),
        .groups = 'drop'
      )
    
    n_participants_with_regression <- sum(pds_participant_summary$has_any_regression, na.rm = TRUE)
    n_total_participants <- nrow(pds_participant_summary)
    
    # Get observation-level summary
    regression_obs <- pds_check %>%
      dplyr::filter(is_regression == TRUE)
    
    n_observations_with_regression <- nrow(regression_obs)
    
    # Get session breakdown - regressions at each session relative to total obs at that session
    session_breakdown_combined <- regression_obs %>%
      dplyr::group_by(session_id) %>%
      dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
      dplyr::left_join(session_totals_combined, by = "session_id") %>%
      dplyr::mutate(percentage = 100 * count / total_obs)
    
    # Print summary
    cat(sprintf("    Participants with regression: %d / %d (%.1f%%)\n", 
                n_participants_with_regression,
                n_total_participants,
                100 * n_participants_with_regression / n_total_participants))
    cat(sprintf("    Observations with regression: %d / %d (%.1f%%)\n",
                n_observations_with_regression,
                n_total_obs,
                100 * n_observations_with_regression / n_total_obs))
    
    # Save summary statistics
    regression_summary <- rbind(regression_summary, data.frame(
      instrument = instrument$name,
      analysis_type = "Combined",
      sex = "Both",
      n_participants_regression = n_participants_with_regression,
      n_total_participants = n_total_participants,
      pct_participants_regression = 100 * n_participants_with_regression / n_total_participants,
      n_observations_regression = n_observations_with_regression,
      n_total_observations = n_total_obs,
      pct_observations_regression = 100 * n_observations_with_regression / n_total_obs
    ))
    
    # Save all trajectories (including non-regression)
    # Select only common columns to avoid column name mismatches when rbinding
    trajectories_temp <- pds_full_dataset %>%
      select(participant_id, session_id, ab_g_stc__cohort_sex, pds_numeric) %>%
      mutate(instrument = instrument$name, sex = ab_g_stc__cohort_sex)
    all_trajectories <- rbind(all_trajectories, trajectories_temp)
    
    # Save detailed regression observations
    if (n_observations_with_regression > 0) {
      # Select only common columns to avoid column name mismatches when rbinding
      regression_details_temp <- regression_obs %>%
        select(participant_id, session_id, ab_g_stc__cohort_sex, pds_numeric, is_regression) %>%
        mutate(instrument = instrument$name,
               analysis_type = "Combined",
               sex = ab_g_stc__cohort_sex)
      regression_details <- rbind(regression_details, regression_details_temp)
    }
    
    if (n_observations_with_regression > 0) {
      cat("    Session breakdown of regression observations:\n")
      for (i in 1:nrow(session_breakdown_combined)) {
        cat(sprintf("      %s: %d observations (%.1f%%)\n", 
                    session_breakdown_combined$session_id[i],
                    session_breakdown_combined$count[i],
                    session_breakdown_combined$percentage[i]))
      }
    }
  }
  
  # Now analyze each sex separately
  for (sex_label in c("Females", "Males")) {
    sex_code <- if (sex_label == "Females") "F" else "M"
    pds_var <- if (sex_label == "Females") female_var else male_var
    pds_sex_specific <- if (sex_label == "Females") pds_temp_f else pds_temp_m
    
    cat(sprintf("\n%s:\n", pds_var))
    
    # Sex-specific stats with percentages relative to combined total
    cat(sprintf("  %s (relative to instrument total M+F):\n", sex_label))
    
    if (nrow(pds_full_dataset) == 0) {
      cat("    No combined data available\n")
    } else if (nrow(pds_sex_specific) == 0) {
      cat(sprintf("    No %s data available\n", tolower(sex_label)))
    } else {
      # Calculate total numbers from combined dataset for denominators
      n_total_obs_combined <- nrow(pds_full_dataset)
      
      pds_check_combined <- pds_full_dataset %>%
        dplyr::group_by(participant_id) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup()
      
      pds_participant_summary_combined <- pds_check_combined %>%
        dplyr::group_by(participant_id) %>%
        dplyr::summarise(has_any_regression = any(is_regression <- pds_numeric < dplyr::lag(pds_numeric), na.rm = TRUE),
                        .groups = 'drop')
      
      n_total_participants_combined <- nrow(pds_participant_summary_combined)
      
      # Calculate sex-specific regression stats
      pds_check_sex <- pds_sex_specific %>%
        dplyr::group_by(participant_id) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::mutate(
          is_regression = pds_numeric < dplyr::lag(pds_numeric)
        ) %>%
        dplyr::ungroup()
      
      pds_participant_summary_sex <- pds_check_sex %>%
        dplyr::group_by(participant_id) %>%
        dplyr::summarise(
          has_any_regression = any(is_regression, na.rm = TRUE),
          .groups = 'drop'
        )
      
      n_participants_with_regression_sex <- sum(pds_participant_summary_sex$has_any_regression, na.rm = TRUE)
      
      regression_obs_sex <- pds_check_sex %>%
        dplyr::filter(is_regression == TRUE)
      
      n_observations_with_regression_sex <- nrow(regression_obs_sex)
      
      # Print summary with percentages relative to combined total
      cat(sprintf("    Participants with regression: %d / %d (%.1f%%)\n", 
                  n_participants_with_regression_sex,
                  n_total_participants_combined,
                  100 * n_participants_with_regression_sex / n_total_participants_combined))
      cat(sprintf("    Observations with regression: %d / %d (%.1f%%)\n",
                  n_observations_with_regression_sex,
                  n_total_obs_combined,
                  100 * n_observations_with_regression_sex / n_total_obs_combined))
      
      # Save sex-specific summary (relative to combined)
      regression_summary <- rbind(regression_summary, data.frame(
        instrument = pds_var,
        analysis_type = "Sex-specific (vs combined)",
        sex = sex_label,
        n_participants_regression = n_participants_with_regression_sex,
        n_total_participants = n_total_participants_combined,
        pct_participants_regression = 100 * n_participants_with_regression_sex / n_total_participants_combined,
        n_observations_regression = n_observations_with_regression_sex,
        n_total_observations = n_total_obs_combined,
        pct_observations_regression = 100 * n_observations_with_regression_sex / n_total_obs_combined
      ))
      
      if (n_observations_with_regression_sex > 0) {
        session_breakdown_sex <- regression_obs_sex %>%
          dplyr::group_by(session_id) %>%
          dplyr::summarise(count = dplyr::n(), .groups = 'drop')
        
        # Join with combined total observations to get percentages relative to M+F total at each session
        session_breakdown_relative <- session_breakdown_sex %>%
          dplyr::left_join(session_totals_combined, by = "session_id") %>%
          dplyr::mutate(percentage = 100 * count / total_obs)
        
        cat("    Session breakdown of regression observations:\n")
        for (i in 1:nrow(session_breakdown_relative)) {
          cat(sprintf("      %s: %d observations (%.1f%%)\n", 
                      session_breakdown_relative$session_id[i],
                      session_breakdown_relative$count[i],
                      session_breakdown_relative$percentage[i]))
        }
      }
    }
    
    # Sex-specific stats with percentages relative to sex-specific total
    cat(sprintf("  %s only:\n", sex_label))
    
    if (nrow(pds_sex_specific) == 0) {
      cat("    No data available\n")
    } else {
      n_total_obs_sex <- nrow(pds_sex_specific)
      
      # Get total observations by session for sex-specific data
      session_totals_sex <- pds_sex_specific %>%
        dplyr::group_by(session_id) %>%
        dplyr::summarise(total_obs = dplyr::n(), .groups = 'drop')
      
      # Group by participant and identify regressions
      pds_check <- pds_sex_specific %>%
        dplyr::group_by(participant_id) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::mutate(
          is_regression = pds_numeric < dplyr::lag(pds_numeric)
        ) %>%
        dplyr::ungroup()
      
      # Get participant-level summary
      pds_participant_summary <- pds_check %>%
        dplyr::group_by(participant_id) %>%
        dplyr::summarise(
          has_any_regression = any(is_regression, na.rm = TRUE),
          .groups = 'drop'
        )
      
      n_participants_with_regression <- sum(pds_participant_summary$has_any_regression, na.rm = TRUE)
      n_total_participants <- nrow(pds_participant_summary)
      
      # Get observation-level summary
      regression_obs <- pds_check %>%
        dplyr::filter(is_regression == TRUE)
      
      n_observations_with_regression <- nrow(regression_obs)
      
      # Get session breakdown - regressions at each session relative to total sex-specific obs at that session
      session_breakdown <- regression_obs %>%
        dplyr::group_by(session_id) %>%
        dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
        dplyr::left_join(session_totals_sex, by = "session_id") %>%
        dplyr::mutate(percentage = 100 * count / total_obs)
      
      # Print summary
      cat(sprintf("    Participants with regression: %d / %d (%.1f%%)\n", 
                  n_participants_with_regression,
                  n_total_participants,
                  100 * n_participants_with_regression / n_total_participants))
      cat(sprintf("    Observations with regression: %d / %d (%.1f%%)\n",
                  n_observations_with_regression,
                  n_total_obs_sex,
                  100 * n_observations_with_regression / n_total_obs_sex))
      
      # Save sex-specific only summary
      regression_summary <- rbind(regression_summary, data.frame(
        instrument = pds_var,
        analysis_type = "Sex-specific only",
        sex = sex_label,
        n_participants_regression = n_participants_with_regression,
        n_total_participants = n_total_participants,
        pct_participants_regression = 100 * n_participants_with_regression / n_total_participants,
        n_observations_regression = n_observations_with_regression,
        n_total_observations = n_total_obs_sex,
        pct_observations_regression = 100 * n_observations_with_regression / n_total_obs_sex
      ))
      
      if (n_observations_with_regression > 0) {
        cat("    Session breakdown of regression observations:\n")
        for (i in 1:nrow(session_breakdown)) {
          cat(sprintf("      %s: %d observations (%.1f%%)\n", 
                      session_breakdown$session_id[i],
                      session_breakdown$count[i],
                      session_breakdown$percentage[i]))
        }
      }
    }
  }
}
cat("\n")

# Save PDS regression data to CSV files
cat("--- Saving PDS Regression Data ---\n")
pds_output_dir <- file.path(outpath, "pds_regression_analysis")
if (!dir.exists(pds_output_dir)) {
  dir.create(pds_output_dir, recursive = TRUE)
}

# Save summary statistics
write.csv(regression_summary, 
          file.path(pds_output_dir, "pds_regression_summary.csv"), 
          row.names = FALSE)
cat(sprintf("Saved summary to: %s/pds_regression_summary.csv\n", pds_output_dir))

# Save detailed regression observations
if (nrow(regression_details) > 0) {
  write.csv(regression_details, 
            file.path(pds_output_dir, "pds_regression_details.csv"), 
            row.names = FALSE)
  cat(sprintf("Saved regression details to: %s/pds_regression_details.csv\n", pds_output_dir))
}

# Save all trajectories
write.csv(all_trajectories, 
          file.path(pds_output_dir, "pds_all_trajectories.csv"), 
          row.names = FALSE)
cat(sprintf("Saved all trajectories to: %s/pds_all_trajectories.csv\n", pds_output_dir))

cat("\n")

# Now compute averaged PDS scores
pds_youth <- pds_youth %>%
  mutate(pds_categ_y_all = as.numeric(coalesce(ph_y_pds__m_categ, ph_y_pds__f_categ)))

pds_parent <- pds_parent %>%
  mutate(pds_categ_p_all = as.numeric(coalesce(ph_p_pds__m_categ, ph_p_pds__f_categ)))

pds_combined <- full_join(pds_youth, pds_parent, by = c("participant_id", "session_id"))

pds_combined <- pds_combined %>%
  mutate(pds_categ_y_p_all = rowMeans(select(., pds_categ_y_all, pds_categ_p_all), na.rm = TRUE))

outmat <- join(outmat, pds_combined[, c('participant_id', 'session_id', 'pds_categ_y_all', 'pds_categ_p_all', 'pds_categ_y_p_all')], 
            by = c('participant_id', 'session_id'))

################################
# Load phys file and extract height and weight for BMI computation
phys = arrow::read_parquet(physfile)
# Select relevant columns and compute BMI and TMI
physvars <- c('participant_id', 'session_id', 'ph_y_anthr__height_mean', 'ph_y_anthr__weight_mean')
phys <- phys[, physvars]

phys <- phys %>%
  mutate(
        weight_lb = as.numeric(ph_y_anthr__weight_mean),
        height_in = as.numeric(ph_y_anthr__height_mean),
        weightkg = weight_lb * 0.453592,
        heightm = height_in * 0.0254,
        anthro_bmi_corr = weightkg / (heightm^2),
        anthro_tmi_calc = weightkg / (heightm^3)
  ) %>%
  select(participant_id, session_id, anthro_bmi_corr, weightkg, heightm)

outmat <- join(outmat, phys, by = c('participant_id', 'session_id'))

################################
################################
# ANTHROPOMETRIC QC (decreasing height, then percentile extremes, then duplicate report)
################################

cat("\n========================================\n")
cat("STARTING ANTHROPOMETRIC QUALITY CONTROL\n")
cat("========================================\n\n")

# Remove rows where height, weight, or BMI is missing (initial cleanup)
cat("--- Initial Missing Data Removal ---\n")
n_before_missing <- nrow(outmat)
outmat <- outmat %>%
  filter(!is.na(heightm) & !is.na(weightkg) & !is.na(anthro_bmi_corr))
n_after_missing <- nrow(outmat)
cat(sprintf("Removed %d observations with missing height/weight/BMI\n", n_before_missing - n_after_missing))
cat(sprintf("Remaining observations: %d\n\n", n_after_missing))

# Longitudinal height plausibility before cross-sectional percentile trimming
outmat <- outmat %>% dplyr::arrange(participant_id, ab_g_dyn__visit_age)

cat("--- Decreasing Height QC ---\n")

tmpAllHeight <- outmat %>%
  dplyr::group_by(participant_id) %>%
  dplyr::arrange(ab_g_dyn__visit_age) %>%
  dplyr::summarise(
    heights = list(heightm),
    is_sorted = list(heightm == cummax(heightm)),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(has_decrease = sapply(is_sorted, function(x) !all(x)))

participants_with_dec_height <- tmpAllHeight %>%
  dplyr::filter(has_decrease) %>%
  dplyr::pull(participant_id)

cat(sprintf("Participants with decreasing height: %d\n", length(participants_with_dec_height)))

rows_to_remove <- c()
if (length(participants_with_dec_height) > 0) {
  for (pid in participants_with_dec_height) {
    participant_data <- outmat %>%
      dplyr::filter(participant_id == pid) %>%
      dplyr::arrange(ab_g_dyn__visit_age)

    tmpHt <- participant_data$heightm
    tmpAge <- participant_data$ab_g_dyn__visit_age

    if (!identical(tmpAge, sort(tmpAge))) {
      stop("Age is not sorted for participant ", pid)
    }

    if (identical(tmpHt, sort(tmpHt))) {
      stop("Height is already sorted for participant ", pid, "; possible bug!")
    }

    sorted_order <- order(tmpHt)
    position_diffs <- seq_along(tmpHt) - sorted_order
    rows_out_of_order <- which(position_diffs != 0)

    if (length(rows_out_of_order) > 0) {
      all_indices <- which(outmat$participant_id == pid)
      indices_to_remove <- all_indices[rows_out_of_order]
      rows_to_remove <- c(rows_to_remove, indices_to_remove)
    }
  }
}

n_removed_dec_height <- length(rows_to_remove)
if (n_removed_dec_height > 0) {
  outmat <- outmat[-rows_to_remove, ]
}

cat(sprintf("Removed %d observations due to decreasing height\n", n_removed_dec_height))
cat(sprintf("Remaining observations: %d\n\n", nrow(outmat)))

# Percentile-based extreme value removal by age bins and sex
cat("--- Percentile-Based Extreme Value Removal ---\n")

# Define age bins (matching MATLAB: 8-10, 10-11, 11-12, ..., 16-18)
age_bins_lower <- c(8, 10:16)
age_bins_upper <- c(10:16, 18)

# Split by sex
outmat_males <- outmat %>% filter(ab_g_stc__cohort_sex == "M")
outmat_females <- outmat %>% filter(ab_g_stc__cohort_sex == "F")

cat(sprintf("Processing %d age bins for males and females\n", length(age_bins_lower)))

# Save original values before percentile QC (to determine direction later)
outmat_males$orig_heightm <- outmat_males$heightm
outmat_males$orig_weightkg <- outmat_males$weightkg
outmat_males$orig_bmi <- outmat_males$anthro_bmi_corr
outmat_females$orig_heightm <- outmat_females$heightm
outmat_females$orig_weightkg <- outmat_females$weightkg
outmat_females$orig_bmi <- outmat_females$anthro_bmi_corr

# Process each age bin
for (i in seq_along(age_bins_lower)) {
  age_lower <- age_bins_lower[i]
  age_upper <- age_bins_upper[i]
  
  # Filter data for this age bin
  males_age_bin <- outmat_males %>%
    filter(ab_g_dyn__visit_age >= age_lower & ab_g_dyn__visit_age < age_upper)
  
  females_age_bin <- outmat_females %>%
    filter(ab_g_dyn__visit_age >= age_lower & ab_g_dyn__visit_age < age_upper)
  
  # Calculate and apply percentiles for males
  if (nrow(males_age_bin) > 0) {
    # Height
    ht_lower_m <- quantile(males_age_bin$heightm, 0.005, na.rm = TRUE)
    ht_upper_m <- quantile(males_age_bin$heightm, 0.995, na.rm = TRUE)
    
    # Weight
    wt_lower_m <- quantile(males_age_bin$weightkg, 0.005, na.rm = TRUE)
    wt_upper_m <- quantile(males_age_bin$weightkg, 0.995, na.rm = TRUE)
    
    # BMI
    bmi_lower_m <- quantile(males_age_bin$anthro_bmi_corr, 0.005, na.rm = TRUE)
    bmi_upper_m <- quantile(males_age_bin$anthro_bmi_corr, 0.995, na.rm = TRUE)
    
    # Apply to males
    loc_age_males <- outmat_males$ab_g_dyn__visit_age >= age_lower & 
                     outmat_males$ab_g_dyn__visit_age < age_upper
    
    outmat_males$heightm[loc_age_males & (outmat_males$heightm < ht_lower_m | outmat_males$heightm > ht_upper_m)] <- NA
    outmat_males$weightkg[loc_age_males & (outmat_males$weightkg < wt_lower_m | outmat_males$weightkg > wt_upper_m)] <- NA
    outmat_males$anthro_bmi_corr[loc_age_males & (outmat_males$anthro_bmi_corr < bmi_lower_m | outmat_males$anthro_bmi_corr > bmi_upper_m)] <- NA
  }
  
  # Calculate and apply percentiles for females
  if (nrow(females_age_bin) > 0) {
    # Height
    ht_lower_f <- quantile(females_age_bin$heightm, 0.005, na.rm = TRUE)
    ht_upper_f <- quantile(females_age_bin$heightm, 0.995, na.rm = TRUE)
    
    # Weight
    wt_lower_f <- quantile(females_age_bin$weightkg, 0.005, na.rm = TRUE)
    wt_upper_f <- quantile(females_age_bin$weightkg, 0.995, na.rm = TRUE)
    
    # BMI
    bmi_lower_f <- quantile(females_age_bin$anthro_bmi_corr, 0.005, na.rm = TRUE)
    bmi_upper_f <- quantile(females_age_bin$anthro_bmi_corr, 0.995, na.rm = TRUE)
    
    # Apply to females
    loc_age_females <- outmat_females$ab_g_dyn__visit_age >= age_lower & 
                       outmat_females$ab_g_dyn__visit_age < age_upper
    
    outmat_females$heightm[loc_age_females & (outmat_females$heightm < ht_lower_f | outmat_females$heightm > ht_upper_f)] <- NA
    outmat_females$weightkg[loc_age_females & (outmat_females$weightkg < wt_lower_f | outmat_females$weightkg > wt_upper_f)] <- NA
    outmat_females$anthro_bmi_corr[loc_age_females & (outmat_females$anthro_bmi_corr < bmi_lower_f | outmat_females$anthro_bmi_corr > bmi_upper_f)] <- NA
  }
}

# Merge males and females back together
outmat <- bind_rows(outmat_males, outmat_females)

# Classify each flagged observation by measure and direction
outmat <- outmat %>%
  dplyr::mutate(
    ht_flagged  = is.na(heightm),
    wt_flagged  = is.na(weightkg),
    bmi_flagged = is.na(anthro_bmi_corr),
    ht_below  = ht_flagged  & !is.na(orig_heightm) & (orig_heightm < orig_heightm), # placeholder
    ht_above  = ht_flagged  & !is.na(orig_heightm),
    wt_below  = wt_flagged  & !is.na(orig_weightkg),
    wt_above  = wt_flagged  & !is.na(orig_weightkg),
    bmi_below = bmi_flagged & !is.na(orig_bmi),
    bmi_above = bmi_flagged & !is.na(orig_bmi)
  )

# Values below median that were flagged → below 0.5th; above median → above 99.5th.
outmat <- outmat %>%
  dplyr::group_by(ab_g_stc__cohort_sex,
                  age_bin = cut(ab_g_dyn__visit_age,
                                breaks = c(8, 10:16, 18),
                                right = FALSE,
                                include.lowest = TRUE)) %>%
  dplyr::mutate(
    ht_median  = median(orig_heightm, na.rm = TRUE),
    wt_median  = median(orig_weightkg, na.rm = TRUE),
    bmi_median = median(orig_bmi, na.rm = TRUE),
    ht_below   = ht_flagged  & (orig_heightm < ht_median),
    ht_above   = ht_flagged  & (orig_heightm >= ht_median),
    wt_below   = wt_flagged  & (orig_weightkg < wt_median),
    wt_above   = wt_flagged  & (orig_weightkg >= wt_median),
    bmi_below  = bmi_flagged & (orig_bmi < bmi_median),
    bmi_above  = bmi_flagged & (orig_bmi >= bmi_median)
  ) %>%
  dplyr::ungroup()

# --- Per-measure summary by session and sex ---
for (measure in c("Height", "Weight", "BMI")) {
  flag_col  <- switch(measure, Height = "ht_flagged",  Weight = "wt_flagged",  BMI = "bmi_flagged")
  below_col <- switch(measure, Height = "ht_below",    Weight = "wt_below",    BMI = "bmi_below")
  above_col <- switch(measure, Height = "ht_above",    Weight = "wt_above",    BMI = "bmi_above")

  cat(sprintf("\n  %s extremes by session and sex:\n", measure))
  summ <- outmat %>%
    dplyr::group_by(session_id, ab_g_stc__cohort_sex) %>%
    dplyr::summarise(
      n_total    = dplyr::n(),
      n_flagged  = sum(.data[[flag_col]], na.rm = TRUE),
      n_below    = sum(.data[[below_col]], na.rm = TRUE),
      n_above    = sum(.data[[above_col]], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(session_id, ab_g_stc__cohort_sex)

  for (i in 1:nrow(summ)) {
    cat(sprintf("    %s | %s: %d / %d flagged (%.1f%%) — %d below 0.5th, %d above 99.5th\n",
                summ$session_id[i], summ$ab_g_stc__cohort_sex[i],
                summ$n_flagged[i], summ$n_total[i],
                100 * summ$n_flagged[i] / summ$n_total[i],
                summ$n_below[i], summ$n_above[i]))
  }
}

# --- Overlap summary ---
cat("\n  Overlap of flagged observations:\n")
overlap <- outmat %>%
  dplyr::filter(ht_flagged | wt_flagged | bmi_flagged) %>%
  dplyr::summarise(
    total_flagged     = dplyr::n(),
    ht_only           = sum(ht_flagged & !wt_flagged & !bmi_flagged),
    wt_only           = sum(!ht_flagged & wt_flagged & !bmi_flagged),
    bmi_only          = sum(!ht_flagged & !wt_flagged & bmi_flagged),
    ht_and_wt_only    = sum(ht_flagged & wt_flagged & !bmi_flagged),
    ht_and_bmi_only   = sum(ht_flagged & !wt_flagged & bmi_flagged),
    wt_and_bmi_only   = sum(!ht_flagged & wt_flagged & bmi_flagged),
    all_three         = sum(ht_flagged & wt_flagged & bmi_flagged)
  )
cat(sprintf("    Total flagged observations: %d\n", overlap$total_flagged))
cat(sprintf("    Height only:        %d\n", overlap$ht_only))
cat(sprintf("    Weight only:        %d\n", overlap$wt_only))
cat(sprintf("    BMI only:           %d\n", overlap$bmi_only))
cat(sprintf("    Height + Weight:    %d\n", overlap$ht_and_wt_only))
cat(sprintf("    Height + BMI:       %d\n", overlap$ht_and_bmi_only))
cat(sprintf("    Weight + BMI:       %d\n", overlap$wt_and_bmi_only))
cat(sprintf("    All three:          %d\n", overlap$all_three))
cat("\n")

# Clean up temp columns
outmat <- outmat %>%
  dplyr::select(-orig_heightm, -orig_weightkg, -orig_bmi,
                -ht_flagged, -wt_flagged, -bmi_flagged,
                -ht_below, -ht_above, -wt_below, -wt_above,
                -bmi_below, -bmi_above, -ht_median, -wt_median,
                -bmi_median, -age_bin)

# Remove rows where height, weight, or BMI is now missing after percentile QC
n_before_removal <- nrow(outmat)
outmat <- outmat %>%
  filter(!is.na(heightm) & !is.na(weightkg) & !is.na(anthro_bmi_corr))
n_after_removal <- nrow(outmat)

cat(sprintf("Removed %d observations with extreme values (0.5/99.5 percentiles)\n", n_before_removal - n_after_removal))
cat(sprintf("Remaining observations: %d\n\n", n_after_removal))

# Sort by participant_id and age
outmat <- outmat %>% arrange(participant_id, ab_g_dyn__visit_age)

# Check for duplicate anthropometric values when age changes
cat("--- Duplicate Anthropometric Value Check ---\n")

dup_check <- outmat %>%
  dplyr::group_by(participant_id) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_unique_age = dplyr::n_distinct(ab_g_dyn__visit_age),
    n_unique_bmi = dplyr::n_distinct(anthro_bmi_corr, na.rm = TRUE),
    n_unique_height = dplyr::n_distinct(heightm, na.rm = TRUE),
    n_unique_weight = dplyr::n_distinct(weightkg, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    dup_bmi = n_obs > n_unique_bmi & n_unique_age == n_obs,
    dup_height = n_obs > n_unique_height & n_unique_age == n_obs,
    dup_weight = n_obs > n_unique_weight & n_unique_age == n_obs
  )

n_dup_bmi <- sum(dup_check$dup_bmi, na.rm = TRUE)
n_dup_height <- sum(dup_check$dup_height, na.rm = TRUE)
n_dup_weight <- sum(dup_check$dup_weight, na.rm = TRUE)

cat(sprintf("Participants with duplicate BMI (age changes): %d\n", n_dup_bmi))
cat(sprintf("Participants with duplicate height (age changes): %d\n", n_dup_height))
cat(sprintf("Participants with duplicate weight (age changes): %d\n", n_dup_weight))
cat("(These are flagged but not removed)\n\n")

cat("========================================\n")
cat("ANTHROPOMETRIC QC COMPLETE\n")
cat("========================================\n\n")

################################
# Compute CDC BMI Percentiles
################################

cat("--- Computing CDC BMI Percentiles ---\n")

# Prepare data for cdcanthro
outmat <- outmat %>%
  mutate(
    age_in_months = ab_g_dyn__visit_age * 12,
    sex = ab_g_stc__cohort_sex,
    height_cm = heightm * 100,
    weight_kg = weightkg,
    bmi = anthro_bmi_corr
  )

# Compute percentiles (cdcanthro requires data.frame, not tibble)
# all = TRUE needed for adj_perc_median (age-20-adjusted % from median)
outmat <- cdcanthro(
  data = as.data.frame(outmat),
  age = age_in_months,
  wt = weight_kg,
  ht = height_cm,
  bmi = bmi,
  all = TRUE
)

# Keep CDC BMI metrics (bmiz = extended z-score paired with bmip); drop other cdcanthro outputs
outmat <- outmat %>%
  select(-age_in_months, -sex, -height_cm, -weight_kg, -bmi,
         -waz, -wap, -haz, -hap, -p50, -p95,
         -original_bmiz,
         -mod_bmiz, -mod_waz, -mod_haz,
         -bmi_l, -bmi_m, -bmi_s, -sigma,
         -dist_median, -adj_dist_median, -log_perc_median, -adj_log_perc_median)

cat(sprintf(
  "Added bmip, bmiz, original_bmip, perc_median, bmip95 for %d observations\n\n",
  nrow(outmat)
))

################################
# Apply session filter, discrepant CSV exclusion, and covariate completeness filter
################################

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# Keep only selected sessions
cat("--- Applying Session Filter ---\n")
keep_sessions <- c("ses-00A", "ses-02A", "ses-04A", "ses-06A")
n_before_session <- nrow(outmat)
outmat <- outmat[outmat$session_id %in% keep_sessions, ]
n_after_session <- nrow(outmat)
cat(sprintf("Kept only %s sessions\n", paste(keep_sessions, collapse = ", ")))
cat(sprintf("Removed %d observations\n", n_before_session - n_after_session))
cat(sprintf("Remaining observations: %d\n\n", n_after_session))

# Remove observations involved in regression of averaged PDS across keep_sessions
cat("--- Removing PDS Average Regression Pairs (keep_sessions only) ---\n")

pds_reg_pairs <- outmat %>%
  dplyr::filter(session_id %in% keep_sessions) %>%
  dplyr::mutate(session_id = factor(session_id, levels = keep_sessions, ordered = TRUE)) %>%
  dplyr::arrange(participant_id, session_id) %>%
  dplyr::group_by(participant_id) %>%
  dplyr::mutate(
    prev_pds = dplyr::lag(pds_categ_y_p_all),
    prev_session = dplyr::lag(as.character(session_id)),
    is_regression = !is.na(prev_pds) & !is.na(pds_categ_y_p_all) &
      (pds_categ_y_p_all < prev_pds)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is_regression) %>%
  dplyr::transmute(
    participant_id,
    session_before = prev_session,
    session_after = as.character(session_id),
    pds_before = prev_pds,
    pds_after = pds_categ_y_p_all
  )

if (nrow(pds_reg_pairs) > 0) {
  pds_regression_obs <- dplyr::bind_rows(
    pds_reg_pairs %>% dplyr::transmute(participant_id, session_id = session_before),
    pds_reg_pairs %>% dplyr::transmute(participant_id, session_id = session_after)
  ) %>% dplyr::distinct()

  n_before_pds_reg <- nrow(outmat)
  outmat <- dplyr::anti_join(outmat, pds_regression_obs, by = c("participant_id", "session_id"))
  n_removed_pds_reg <- n_before_pds_reg - nrow(outmat)

  cat(sprintf("Removed %d observations involved in PDS average regression pairs\n", n_removed_pds_reg))
  cat(sprintf("Regression pairs identified: %d\n", nrow(pds_reg_pairs)))
  cat(sprintf("Remaining observations: %d\n\n", nrow(outmat)))
} else {
  cat("No PDS average regression pairs detected in keep_sessions\n\n")
}

# Exclude observations listed in discrepant ages CSV (participant_id + session_id)
cat("--- Applying Discrepant Age Exclusion ---\n")
disc_path <- bmi_age_discrepancy_csv_path(PATHS)
if (file.exists(disc_path)) {
  n_before_disc <- nrow(outmat)
  suppressWarnings({
    disc <- readr::read_csv(disc_path, show_col_types = FALSE)
  })
  keys <- intersect(c('participant_id','session_id'), names(disc))
  if (all(c('participant_id','session_id') %in% keys)) {
    outmat <- dplyr::anti_join(outmat, disc[, c('participant_id','session_id')], by = c('participant_id','session_id'))
  } else if ('participant_id' %in% names(disc)) {
    outmat <- dplyr::anti_join(outmat, dplyr::distinct(disc['participant_id']), by = 'participant_id')
  }
  n_after_disc <- nrow(outmat)
  cat(sprintf("Removed %d observations from discrepant age CSV\n", n_before_disc - n_after_disc))
  cat(sprintf("Remaining observations: %d\n\n", n_after_disc))
} else {
  cat("Discrepant age CSV not found, skipping\n\n")
}

# Covariate completeness filter
cat("--- Applying Covariate Completeness Filter ---\n")
# Use dyn-sourced scanner/software and edu/income6
is_missing_scanner <- is.na(outmat$ab_g_dyn__design_mr__serial_software)
is_missing_education <- is.na(outmat$ab_g_dyn__cohort_edu__cgs)
is_missing_income6 <- is.na(outmat$ab_g_dyn__cohort_income__hhold__6lvl)

# First 10 genetic PCs
pc_cols <- grep('^ab_g_stc__gen_pc__', names(outmat), value = TRUE)
if (length(pc_cols) > 10) pc_cols <- pc_cols[1:10]
is_missing_genetics <- if (length(pc_cols) > 0) {
  apply(as.data.frame(outmat)[, pc_cols, drop = FALSE], 1, function(r) any(is.na(r)))
} else { rep(TRUE, nrow(outmat)) }

is_missing_bmi <- is.na(outmat$anthro_bmi_corr)
is_missing_pds <- is.na(outmat$pds_categ_y_p_all)

to_delete <- is_missing_scanner | is_missing_education | is_missing_income6 | is_missing_genetics | is_missing_bmi | is_missing_pds
cat(sprintf('Removed %d observations for missing covariates\n', sum(to_delete)))
cat(sprintf('  - Missing scanner: %d\n', sum(is_missing_scanner)))
cat(sprintf('  - Missing education: %d\n', sum(is_missing_education)))
cat(sprintf('  - Missing income: %d\n', sum(is_missing_income6)))
cat(sprintf('  - Missing genetics PCs: %d\n', sum(is_missing_genetics)))
cat(sprintf('  - Missing BMI: %d\n', sum(is_missing_bmi)))
cat(sprintf('  - Missing PDS: %d\n', sum(is_missing_pds)))

outmat_filtered <- outmat[!to_delete, ]
cat(sprintf("Remaining observations: %d\n", nrow(outmat_filtered)))
cat(sprintf("Unique participants: %d\n\n", n_distinct(outmat_filtered$participant_id)))

# Write broad pheno-QC table (mapping / exploration)
bmi_assert_positive_nrow(outmat_filtered, "Filtered output (outmat_filtered)")
write.table(outmat_filtered, file = paste0(outmatfile_pheno_qc, ".txt"), sep = "\t", row.names = FALSE)
cat(paste0("File written to ", outmatfile_pheno_qc, ".txt\n"))
cat(sprintf("  Observations: %d | Unique participants: %d\n\n",
            nrow(outmat_filtered), n_distinct(outmat_filtered$participant_id)))

# Imaging-intended cohort: dMRI included + family ID (feeds s2 onward)
cat("--- Applying dMRI + family ID filter (imaging cohort) ---\n")
bmi_assert_cols(
  outmat_filtered,
  c("mr_y_qc__incl__dmri_indicator", "ab_g_stc__design_id__fam"),
  "Filtered output"
)

is_not_dmri <- is.na(outmat_filtered$mr_y_qc__incl__dmri_indicator) |
  outmat_filtered$mr_y_qc__incl__dmri_indicator != 1
is_missing_fam <- is.na(outmat_filtered$ab_g_stc__design_id__fam)
imaging_keep <- !is_not_dmri & !is_missing_fam

n_before_imaging <- nrow(outmat_filtered)
outmat_imaging <- outmat_filtered[imaging_keep, ]
n_after_imaging <- nrow(outmat_imaging)

cat(sprintf("Removed %d observations for dMRI exclusion (indicator != 1 or NA)\n", sum(is_not_dmri)))
cat(sprintf("Removed %d observations for missing family ID\n", sum(is_missing_fam)))
cat(sprintf("Removed %d observations total for imaging filter\n", n_before_imaging - n_after_imaging))
cat(sprintf("Remaining observations: %d\n", n_after_imaging))
cat(sprintf("Unique participants: %d\n\n", n_distinct(outmat_imaging$participant_id)))

bmi_assert_positive_nrow(outmat_imaging, "Imaging output (outmat_imaging)")
write.table(outmat_imaging, file = paste0(outmatfile_imaging, ".txt"), sep = "\t", row.names = FALSE)
cat(paste0("File written to ", outmatfile_imaging, ".txt\n\n"))