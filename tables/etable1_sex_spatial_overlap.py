
"""eTable 1. Female–male spatial overlap of voxelwise BMI associations.

Full-sample only. Writes the published table:
  Phenotype | Pearson correlation | Dice coefficient (FDR)

Pearson: unthresholded omnibus Wald maps (vol_W_s) on the male∩female in-mask voxels.
Dice: FDR-significant voxels (vol_pW_s_FDR == 0).

Env:
  BMI_IMAGING_CEPH_ROOT
  BMI_IMAGING_RESULTS_RUN   e.g. results_2026_06_29
  BMI_PHENO_VARIABLE        bmip | perc_median | both (default both)
  BMI_IMAGING_S7_FSTEM      default rni
"""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PHENO_LABELS = {
    "bmip": "BMI percentile",
    "perc_median": "Percent from median",
}


def env(name: str, default: str = "") -> str:
    return os.environ.get(name, default).strip()


def xname_for_pheno(pheno: str) -> str:
    return "bmi" if pheno == "bmip" else "perc_median"


def find_wald_mat(results_dir: Path, sex: str, modality: str, xname: str) -> Path:
    if sex == "female":
        matches = [
            d
            for d in sorted(results_dir.glob("designmat*f_*"))
            if d.is_dir() and "female" in d.name
        ]
    else:
        matches = [
            d
            for d in sorted(results_dir.glob("designmat*m_*"))
            if d.is_dir() and "male" in d.name and "female" not in d.name
        ]
    if not matches:
        raise FileNotFoundError(
            f"No sex-stratified designmat dir for {sex} under {results_dir}"
        )
    wald = matches[0] / f"{modality}_Wald_splines_CI_output_{xname}.mat"
    if not wald.is_file():
        raise FileNotFoundError(f"Missing Wald mat: {wald}")
    return wald


def load_h5_array(path: Path, key: str) -> np.ndarray:
    with h5py.File(path, "r") as h:
        if key not in h:
            raise KeyError(f"{key} not in {path}; keys={list(h.keys())}")
        return np.asarray(h[key])


def inmask_from_pval_vol(vol_p: np.ndarray) -> np.ndarray:
    p = np.asarray(vol_p, dtype=float)
    return np.isfinite(p) & (p > 0.0) & (p <= 1.0)


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def dice(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    n_a = int(np.sum(a))
    n_b = int(np.sum(b))
    n_and = int(np.sum(a & b))
    return (2.0 * n_and / (n_a + n_b)) if (n_a + n_b) else float("nan")


def overlap_metrics(mat_m: Path, mat_f: Path) -> tuple[float, float]:
    vol_w_m = load_h5_array(mat_m, "vol_W_s").reshape(-1)
    vol_w_f = load_h5_array(mat_f, "vol_W_s").reshape(-1)
    vol_p_m = load_h5_array(mat_m, "vol_pW_s").reshape(-1)
    vol_p_f = load_h5_array(mat_f, "vol_pW_s").reshape(-1)
    vol_fdr_m = load_h5_array(mat_m, "vol_pW_s_FDR").reshape(-1)
    vol_fdr_f = load_h5_array(mat_f, "vol_pW_s_FDR").reshape(-1)

    if vol_w_m.size != vol_w_f.size:
        raise ValueError(
            f"Volume size mismatch: male={vol_w_m.size} female={vol_w_f.size}"
        )

    inter = inmask_from_pval_vol(vol_p_m) & inmask_from_pval_vol(vol_p_f)
    if int(np.sum(inter)) == 0:
        raise ValueError("Empty intersection mask")

    r = pearson(vol_w_m[inter], vol_w_f[inter])
    d = dice(vol_fdr_m[inter] == 0, vol_fdr_f[inter] == 0)
    return r, d


def main() -> None:
    ceph = Path(
        env(
            "BMI_IMAGING_CEPH_ROOT",
            "/space/ceph/1/ABCD/users/arigby/2026_bmi_imaging_analysis",
        )
    )
    results_run = env("BMI_IMAGING_RESULTS_RUN") or f"results_{date.today().strftime('%Y_%m_%d')}"
    pheno_arg = env("BMI_PHENO_VARIABLE", "both")
    phenos = ["bmip", "perc_median"] if pheno_arg.lower() == "both" else [pheno_arg]
    modality = env("BMI_IMAGING_S7_FSTEM", "rni")

    rows = []
    for pheno in phenos:
        xname = xname_for_pheno(pheno)
        results_dir = ceph / pheno / "results" / results_run
        try:
            mat_m = find_wald_mat(results_dir, "male", modality, xname)
            mat_f = find_wald_mat(results_dir, "female", modality, xname)
            r, d = overlap_metrics(mat_m, mat_f)
        except (FileNotFoundError, KeyError, ValueError) as e:
            print(f"SKIP {pheno}: {e}")
            continue

        rows.append(
            {
                "Phenotype": PHENO_LABELS.get(pheno, pheno),
                "Pearson correlation": round(r, 3),
                "Dice coefficient (FDR)": round(d, 3),
            }
        )
        print(f"{pheno}: Pearson={r:.3f}  Dice FDR={d:.3f}")

    if not rows:
        raise SystemExit(
            "No overlaps completed — check BMI_IMAGING_RESULTS_RUN and that s7 finished."
        )

    out = pd.DataFrame(rows)
    out_dir = ceph / "sex_spatial_overlap" / results_run
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "sex_spatial_overlap_manuscript.csv"
    out.to_csv(out_csv, index=False)
    print(f"\nManuscript table: {out_csv}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
