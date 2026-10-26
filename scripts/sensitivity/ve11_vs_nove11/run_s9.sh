#!/usr/bin/env bash
# s9 ROI-averaged splines for VE11 vs noVE11 (both phenos × both arms).
set -euo pipefail

REPO="${BMI_IMAGING_REPO_ROOT:-/home/arigby/projects/ABCD_bmi_microstructure}"
cd "$REPO"

MATLAB="${MATLAB_CMD:-matlab -batch}"
PHENOS="${BMI_VE11_NOVE11_PHENOS:-bmip perc_median}"
RUN_BMIP="${BMI_IMAGING_RESULTS_RUN_BMIP:-results_2026_09_02}"
RUN_PERC="${BMI_IMAGING_RESULTS_RUN_PERC_MEDIAN:-results_2026_10_02}"

export BMI_IMAGING_AUTO_RUN_DATE=0
export BMI_IMAGING_BASE_DATAFILE_IMAGING=datafile_matched_combined.txt

for pheno in $PHENOS; do
  case "$pheno" in
    bmip) export BMI_IMAGING_RESULTS_RUN="$RUN_BMIP" ;;
    perc_median) export BMI_IMAGING_RESULTS_RUN="$RUN_PERC" ;;
    *)
      echo "Unknown pheno: $pheno" >&2
      exit 1
      ;;
  esac
  export BMI_PHENO_VARIABLE="$pheno"
  for arm in VE11 noVE11; do
    export BMI_SCANNER_ARM="$arm"
    # shellcheck source=/dev/null
    source env/run_ve11_vs_nove11_arm.env
    echo "=== s9 $pheno $arm ($BMI_IMAGING_RESULTS_RUN) ==="
    $MATLAB "run('scripts/matlab_add_paths.m'); run('scripts/s9_batch_compute_roi_averaged_splines.m')"
  done
done

echo "All requested s9 jobs finished."
