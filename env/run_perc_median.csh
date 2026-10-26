# tcsh/csh: source env/run_perc_median.csh  
if ($?BMI_IMAGING_REPO_ROOT) then
  set _bmi_env_dir = "$BMI_IMAGING_REPO_ROOT/env"
else
  set _bmi_env_dir = /home/arigby/projects/ABCD_bmi_microstructure/env
endif
source "${_bmi_env_dir}/run_shared.csh"
source "${_bmi_env_dir}/run_date.csh"
setenv BMI_PHENO_VARIABLE perc_median
setenv BMI_IMAGING_ANALYSIS_DIR perc_median
setenv BMI_IMAGING_PHENO_GRID_STEP 1
