#!/usr/bin/env bash
# s2–s4 for bmip and perc_median on the matched combined sample.
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
  # shellcheck source=/dev/null
  source env/run_ve11_vs_nove11_pheno.env
  echo "=== s2 $pheno ==="
  "$Rscript" scripts/s2_create_var_vecs.R
  echo "=== s3 $pheno (matched-sample age knots) ==="
  $MATLAB "run('scripts/matlab_add_paths.m'); run('${SENS}/s3_create_basis_vars.m')"
  echo "=== s4 $pheno ==="
  "$Rscript" scripts/s4_join_bfs_to_datafile.R
done

echo "BFs complete. Next: bash scripts/sensitivity/ve11_vs_nove11/run_fema.sh"
