# BMI and restriction spectrum imaging in ABCD

Analysis code for voxelwise associations between body mass and the restricted normalized signal fraction (RNI) from restriction spectrum imaging (RSI) in the Adolescent Brain Cognitive Development (ABCD) Study.

Two age- and sex-specific BMI phenotypes are analyzed in parallel:

- **BMI percentile** (`bmip`) — CDC extended BMI percentile  
- **Percent from median** (`perc_median`) — $100 \times (\mathrm{BMI}/M - 1)$, where $M$ is the CDC LMS median BMI for age and sex  

Models are sex-stratified mixed-effects models of RNI with natural cubic splines of the BMI phenotype, age, and pubertal development, plus household income, parental education, genetic principal components, and MRI serial_software.

## Data

Analyses use ABCD Release 6.1 tabulated data and concatenated voxelwise RNI measures. Participant-level data are restricted and **not** redistributed with this repository. Obtain data through the NIH Brain Development Cohorts (NBDC) Data Hub (DOI: https://doi.org/10.82525/ab4q-qr87). Point local paths at your tabulated and imaging directories by editing `env/run_shared.csh`.

Design matrices, FEMA results, and derived figure/table inputs live on local analysis storage (configured via `BMI_IMAGING_CEPH_ROOT`), not in this repo.

## Requirements

- **R** — tidyverse / dplyr, arrow, gtsummary, gt 
- **MATLAB** with [cmig_tools / FEMA](https://github.com/cmig-research-group/cmig_tools)  
- **Python 3** — NumPy, pandas, matplotlib, seaborn, SciPy, h5py  

A conda environment with the R and Python packages used here is in [`environment.yml`](environment.yml):

```bash
conda env create -f environment.yml
```

MATLAB and cmig_tools are installed separately.

## Processing pipeline

Before running phenotype-specific steps, source a tcsh profile (edit paths in `env/run_shared.csh` first):

```tcsh
source env/run_shared.csh
source env/run_bmip.csh          # or env/run_perc_median.csh
```

Shared visit table (run once):

| Step | Script | Description |
|------|--------|-------------|
| s0 | `scripts/s0_age_scan_anthro_discrepancy.R` | Flag age / scan / anthropometry mismatches |
| s1 | `scripts/s1_create_datafile6.0.R` | Build the visit-level analysis table |

Then, **per phenotype**:

| Step | Script | Description |
|------|--------|-------------|
| s2 | `scripts/s2_create_var_vecs.R` | Phenotype, age, and PDS evaluation grids |
| s3 | `scripts/s3_create_basis_vars.m` | Natural cubic spline basis functions |
| s4 | `scripts/s4_join_bfs_to_datafile.R` | Merge bases onto the visit table |
| s5 | `scripts/s5_make_design.R` | Combined-sex design matrix |
| s5a | `scripts/s5a_stratify_design_matrix.R` | Male and female design matrices |
| s6 | `scripts/s6_run_FEMA.m` | Voxelwise mixed-effects models |
| s7 | `scripts/s7_batch_compute_wald_stats.m` | Omnibus Wald tests (FDR and Bonferroni) |
| s8 | `scripts/s8_run_showVol.m` | Volume map exports for visualization |
| s9 | `scripts/s9_batch_compute_roi_averaged_splines.m` | ROI-averaged spline curves |
| s9b | `scripts/s9b_rank_roi_wald_stats.m` | Rank ROIs by mean Wald statistic |

After s6 for both phenotypes, build matched analytic-sample tables:

```tcsh
Rscript scripts/helpers/build_model_sample_tables.R
```

Supporting libraries live under `utilities/` (`matlab_utils/`, `r_utils/`) and are imported by the scripts above; they are not run directly.

## Figures and tables

Manuscript outputs are under `figures/` and `tables/`. These require local pipeline results (and ABCD inputs for some panels).

| Item | Entry point |
|------|-------------|
| Figure 1 — BMI metric mapping | `figures/figure1_bmi_mapping.ipynb` |
| Figure 2 — voxelwise maps (`bmip`) | showVol screenshot workflow below |
| Figure 3 — ROI spline panels | `figures/figure3_roi_splines.ipynb` |
| Table 1 — demographics | `tables/table1_demographics.R` |
| Table 2 — gradients | `tables/table2_gradient_csvs.ipynb`, then `tables/table2_gradients.R` |
| eFigures 1–4 — sample distributions | `figures/efigures_sample_distributions.ipynb` |
| eFigure 5 — voxelwise maps (`perc_median`) | same showVol workflow as Figure 2 |
| eTable 1 — sex spatial overlap | `tables/etable1_sex_spatial_overlap.py` |

### Figure 2 / eFigure 5 (showVol screenshots)

After s7/s8 and a fixed results run (e.g. `results_2026_06_29`), in MATLAB:

1. `figures/matlab/showvol_results_for_screenshots.m` — build Wald + spline-derivative volume bundles  
2. `figures/matlab/showvol_screenshots.m` — capture PNGs (uses `utilities/matlab_utils/save_showvol_images.m`)  
3. Optional: `figures/matlab/make_movies.m` — AVI/MP4 from derivative PNG sequences  
4. Compose panels in Adobe Illustrator (`bmip` → Figure 2; `perc_median` → eFigure 5)

## Sensitivity analysis: VE11 vs noVE11

Matched VE11 vs noVE11 comparison (combined-sex models; both `bmip` and `perc_median`). Full runbook:

[`scripts/sensitivity/ve11_vs_nove11/README.md`](scripts/sensitivity/ve11_vs_nove11/README.md)

```bash
bash scripts/sensitivity/ve11_vs_nove11/run_matching.sh
bash scripts/sensitivity/ve11_vs_nove11/run_bfs.sh      # both phenos by default
bash scripts/sensitivity/ve11_vs_nove11/run_fema.sh
bash scripts/sensitivity/ve11_vs_nove11/run_s9.sh
```

| Item | Entry point |
|------|-------------|
| Matched distributions | `figures/efigures_ve11_vs_nove11_matched_distributions.ipynb` |
| ROI spline overlays | `figures/efigures_ve11_vs_nove11_roi_splines.ipynb` |
| Beta scatterplots | `figures/efigures_ve11_vs_nove11_beta_scatterplots.ipynb` |
| showVol screenshots | `figures/matlab/showvol_*_ve11_vs_nove11.m` |
| eTable 2 — demographics | `tables/etable2_ve11_vs_nove11_demographics.R` |

## Repository layout

```
env/              # tcsh / bash environment profiles
scripts/          # main estimation pipeline (s0–s9b)
scripts/helpers/  # analytic-sample glue scripts
scripts/sensitivity/ve11_vs_nove11/   # VE11 vs noVE11 sensitivity track
figures/          # manuscript figure notebooks
figures/matlab/   # showVol screenshot / movie scripts
tables/           # manuscript table scripts
utilities/        # shared R and MATLAB libraries
```

## Citation

Rigby, A., Pecheva, D., Parekh, P., Smith, D. M., Becker, A., Linkersdörfer, J., Watts, R., Loughnan, R., Hagler, D. J., Rhee, K. E., Makowski, C., Dale, A. M., & Jernigan, T. L. (2026). Assessing nonlinear associations between body mass index and brain microstructure across adolescence (p. 2026.04.02.716201). bioRxiv. https://doi.org/10.64898/2026.04.02.716201

