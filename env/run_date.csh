# tcsh/csh: source env/run_date.csh
# Sets dated designmat/results folder names when not already defined.
#
# Override the date tag:
#   setenv BMI_IMAGING_RUN_TAG 2026_06_15
# Disable auto-dating:
#   setenv BMI_IMAGING_AUTO_RUN_DATE 0

if (! $?BMI_IMAGING_AUTO_RUN_DATE) setenv BMI_IMAGING_AUTO_RUN_DATE 1

if ("$BMI_IMAGING_AUTO_RUN_DATE" == "1" || "$BMI_IMAGING_AUTO_RUN_DATE" == "true") then
  if (! $?BMI_IMAGING_RUN_TAG) setenv BMI_IMAGING_RUN_TAG `date +%Y_%m_%d`
  if (! $?BMI_IMAGING_DESIGNMATS_DIR) then
    if (! $?BMI_IMAGING_DESIGNMATS_SUBDIR) setenv BMI_IMAGING_DESIGNMATS_SUBDIR designmats_$BMI_IMAGING_RUN_TAG
  endif
  if (! $?BMI_IMAGING_RESULTS_RUN) setenv BMI_IMAGING_RESULTS_RUN results_$BMI_IMAGING_RUN_TAG
endif
