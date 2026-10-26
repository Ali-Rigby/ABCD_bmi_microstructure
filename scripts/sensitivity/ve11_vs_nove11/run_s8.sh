#!/usr/bin/env bash
# showVol extended (s8) for a single VE11/noVE11 arm.
#
# Usage:
#   bash scripts/sensitivity/ve11_vs_nove11/run_s8.sh VE11
#   BMI_VE11_NOVE11_PHENOS=perc_median bash scripts/sensitivity/ve11_vs_nove11/run_s8.sh noVE11
set -euo pipefail

REPO="${BMI_IMAGING_REPO_ROOT:-/home/arigby/projects/ABCD_bmi_microstructure}"
cd "$REPO"

ARM="${1:-VE11}"
PHENO="${BMI_VE11_NOVE11_PHENOS:-bmip}"
PHENO="${PHENO%% *}"
RESULTS_RUN="${BMI_IMAGING_RESULTS_RUN:-}"
if [[ -z "$RESULTS_RUN" ]]; then
  case "$PHENO" in
    bmip) RESULTS_RUN=results_2026_09_02 ;;
    perc_median) RESULTS_RUN=results_2026_10_02 ;;
    *) RESULTS_RUN=results_2026_09_02 ;;
  esac
fi
MATLAB="${MATLAB_CMD:-matlab -batch}"

export BMI_IMAGING_RESULTS_RUN="$RESULTS_RUN"
export BMI_S8_ANALYSIS_PROFILE="ve11_nove11_${ARM}_${PHENO}"

echo "s8: profile=${BMI_S8_ANALYSIS_PROFILE} run=${RESULTS_RUN}"
$MATLAB "run('scripts/matlab_add_paths.m'); run('scripts/sensitivity/ve11_vs_nove11/s8_run_showVol.m')"
