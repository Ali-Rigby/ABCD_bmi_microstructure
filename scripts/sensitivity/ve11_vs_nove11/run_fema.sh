#!/usr/bin/env bash
# s5 designmats + s6 FEMA + s7 Wald for VE11 and noVE11 × both phenotypes.
set -euo pipefail

REPO="${BMI_IMAGING_REPO_ROOT:-/home/arigby/projects/ABCD_bmi_microstructure}"
cd "$REPO"

export BMI_IMAGING_BASE_DATAFILE_IMAGING=datafile_matched_combined.txt
MATLAB="${MATLAB_CMD:-matlab -batch}"
Rscript="${RSCRIPT:-/home/arigby/programs/anaconda3/envs/r/bin/Rscript}"
PHENOS="${BMI_VE11_NOVE11_PHENOS:-bmip perc_median}"
SENS=scripts/sensitivity/ve11_vs_nove11

for pheno in $PHENOS; do
  export BMI_PHENO_VARIABLE="$pheno"
  for arm in VE11 noVE11; do
    export BMI_SCANNER_ARM="$arm"
    # shellcheck source=/dev/null
    source env/run_ve11_vs_nove11_arm.env
    echo "=== s5 $pheno $arm ==="
    "$Rscript" "${SENS}/s5_make_design_ve11_nove11.R"
    echo "=== s6 $pheno $arm ==="
    $MATLAB "run('scripts/matlab_add_paths.m'); run('scripts/s6_run_FEMA.m')"
    echo "=== s7 $pheno $arm ==="
    $MATLAB "run('scripts/matlab_add_paths.m'); run('scripts/s7_batch_compute_wald_stats.m')"
  done
done

echo "Done."
