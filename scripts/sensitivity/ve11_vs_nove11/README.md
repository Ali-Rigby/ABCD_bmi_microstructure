# Sensitivity analysis: VE11 vs noVE11

BMI- and age-matched comparison of **VE11** vs **noVE11** visits (combined-sex models; both `bmip` and `perc_median`).

## Design

| Item | Description |
|------|--------|
| Population | Full-sample imaging base after main-pipeline s1 |
| Matching | 1:1 on `perc_median` (±10) then age (±0.25 y) |
| Sex | Covariate only (`designmat1_*`; no s5a) |
| Imaging | `BMI_IMAGING_ABCD_SYNC` |
| noVE11 model | `s(pheno)+s(age)+s(pds)+sex+edu+income+PCs+ScanSoft` |
| VE11 model | `s(pheno)+s(age)+s(pds)+sex+edu+income+PCs+ScannerSerial` |

Ceph root: `$BMI_IMAGING_CEPH_ROOT/ve11_vs_noVE11_6_0/`.

## Run order

Prerequisite: main `s0`/`s1` has written  
`$CEPH/shared/datafiles/datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_dmri_incl.txt`.

```bash
cd /home/arigby/projects/ABCD_bmi_microstructure

# Both phenotypes by default (override with BMI_VE11_NOVE11_PHENOS=bmip)
bash scripts/sensitivity/ve11_vs_nove11/run_matching.sh
# Review: figures/efigures_ve11_vs_nove11_matched_distributions.ipynb

bash scripts/sensitivity/ve11_vs_nove11/run_bfs.sh
bash scripts/sensitivity/ve11_vs_nove11/run_fema.sh
bash scripts/sensitivity/ve11_vs_nove11/run_s9.sh

# Optional showVol GUI 
bash scripts/sensitivity/ve11_vs_nove11/run_s8.sh VE11
```

## Manuscript outputs

| Item | Entry point |
|------|-------------|
| Matched distributions | `figures/efigures_ve11_vs_nove11_matched_distributions.ipynb` |
| ROI spline overlays | `figures/efigures_ve11_vs_nove11_roi_splines.ipynb` |
| Beta scatterplots | `figures/efigures_ve11_vs_nove11_beta_scatterplots.ipynb` |
| showVol screenshot bundles | `figures/matlab/showvol_results_for_screenshots_ve11_vs_nove11.m` |
| showVol PNGs | `figures/matlab/showvol_screenshots_ve11_vs_nove11.m` |
| Demographics (analytic sample) | `tables/etable2_ve11_vs_nove11_demographics.R` |

Analytic sample = FEMA `iid`/`eid` visits in **bmip ∩ perc_median** within each arm
(Overall = pooled). Extract visit lists once if missing:

```bash
matlab -batch "run('scripts/matlab_add_paths.m'); run('scripts/sensitivity/ve11_vs_nove11/extract_iid_eid_ve11_vs_nove11.m')"

/home/arigby/programs/anaconda3/envs/sci/bin/Rscript \
  tables/etable2_ve11_vs_nove11_demographics.R
```

