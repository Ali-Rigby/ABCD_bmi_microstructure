#!/usr/bin/env bash
# VE11 vs noVE11: layout + matching.
# Prerequisite: main-pipeline s1 imaging base on $CEPH/shared/datafiles.
set -euo pipefail

REPO="${BMI_IMAGING_REPO_ROOT:-/home/arigby/projects/ABCD_bmi_microstructure}"
cd "$REPO"
# shellcheck source=/dev/null
source env/run_ve11_vs_nove11_shared.env

Rscript="${RSCRIPT:-/home/arigby/programs/anaconda3/envs/r/bin/Rscript}"
if [[ ! -x "$Rscript" ]]; then
  echo "Rscript not found: $Rscript (set RSCRIPT)" >&2
  exit 1
fi

bash scripts/sensitivity/ve11_vs_nove11/ensure_ceph_layout.sh

echo "=== s1 match: perc_median + age (VE11 vs noVE11) ==="
"$Rscript" scripts/sensitivity/ve11_vs_nove11/s1_match_ve11_nove11.R

export BMI_IMAGING_BASE_DATAFILE_IMAGING=datafile_matched_combined.txt

echo ""
echo "Matching complete. Review:"
echo "  ${BMI_IMAGING_CEPH_ROOT}/${BMI_VE11_NOVE11_ANALYSIS_ROOT}/matching/qc/"
echo "  Notebook: figures/efigures_ve11_vs_nove11_matched_distributions.ipynb"
echo ""
echo "Next:"
echo "  bash scripts/sensitivity/ve11_vs_nove11/run_bfs.sh"
echo "  bash scripts/sensitivity/ve11_vs_nove11/run_fema.sh"
