# tcsh/csh: source env/run_ve11_vs_nove11.csh
# Loads shared VE11 vs noVE11 path overrides into the current shell (for interactive work).
# Bash estimation runners source the companion *.env files instead.
if ($?BMI_IMAGING_REPO_ROOT) then
  set _bmi_env_dir = "$BMI_IMAGING_REPO_ROOT/env"
else
  set _bmi_env_dir = /home/arigby/projects/ABCD_bmi_microstructure/env
endif
source "${_bmi_env_dir}/run_shared.csh"
source "${_bmi_env_dir}/run_date.csh"
setenv BMI_VE11_NOVE11_ANALYSIS_ROOT ve11_vs_noVE11_6_0
setenv BMI_IMAGING_SHARED_DATAFILES_DIR ve11_vs_noVE11_6_0/shared/datafiles
setenv BMI_IMAGING_ABCD_SYNC /space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0
unsetenv BMI_IMAGING_VOXELWISE_RSI_DIR
echo "VE11 vs noVE11 (tcsh): BMI_VE11_NOVE11_ANALYSIS_ROOT=$BMI_VE11_NOVE11_ANALYSIS_ROOT"
