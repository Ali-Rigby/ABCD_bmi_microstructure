#!/usr/bin/env bash
# Create Ceph folders for ve11_vs_noVE11_6_0 matched sensitivity analysis.
set -euo pipefail

REPO="${BMI_IMAGING_REPO_ROOT:-/home/arigby/projects/ABCD_bmi_microstructure}"
# shellcheck source=/dev/null
source "${REPO}/env/run_ve11_vs_nove11_shared.env"

CEPH="${BMI_IMAGING_CEPH_ROOT}"
ROOT="${BMI_VE11_NOVE11_ANALYSIS_ROOT}"
DM_SUB="${BMI_IMAGING_DESIGNMATS_SUBDIR:-}"
RES_RUN="${BMI_IMAGING_RESULTS_RUN:-}"

mkdir -p "${CEPH}/${ROOT}/shared/datafiles"
mkdir -p "${CEPH}/${ROOT}/matching/qc"
mkdir -p "${CEPH}/${ROOT}/demographics"
echo "OK ${CEPH}/${ROOT}/shared/datafiles"
echo "OK ${CEPH}/${ROOT}/matching/qc"
echo "OK ${CEPH}/${ROOT}/demographics"

for pheno in bmip perc_median; do
  mkdir -p "${CEPH}/${ROOT}/${pheno}/datafiles"
  echo "OK ${CEPH}/${ROOT}/${pheno}/datafiles"
  for arm in VE11 noVE11; do
    for sub in designmats results; do
      d="${CEPH}/${ROOT}/${pheno}/${arm}/${sub}"
      mkdir -p "$d"
      echo "OK $d"
    done
    if [[ -n "$DM_SUB" ]]; then
      mkdir -p "${CEPH}/${ROOT}/${pheno}/${arm}/designmats/${DM_SUB}"
    fi
    if [[ -n "$RES_RUN" ]]; then
      mkdir -p "${CEPH}/${ROOT}/${pheno}/${arm}/results/${RES_RUN}"
    fi
  done
done

echo "VE11 vs noVE11 layout ready under ${CEPH}/${ROOT}/"
