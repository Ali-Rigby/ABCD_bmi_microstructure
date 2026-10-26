# tcsh/csh: source env/run_bmip.csh 
if ($?BMI_IMAGING_REPO_ROOT) then
  set _bmi_env_dir = "$BMI_IMAGING_REPO_ROOT/env"
else
  set _bmi_env_dir = /home/arigby/projects/ABCD_bmi_microstructure/env
endif
source "${_bmi_env_dir}/run_shared.csh"
source "${_bmi_env_dir}/run_date.csh"
setenv BMI_PHENO_VARIABLE bmip
setenv BMI_IMAGING_ANALYSIS_DIR bmip
# bmip uses fixed 0:100 percentile grid in s3 (101 rows); BMI_IMAGING_PHENO_GRID_STEP is not used
