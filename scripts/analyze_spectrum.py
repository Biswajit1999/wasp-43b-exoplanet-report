from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spectra"
FIGURES = ROOT / "figures"
STATS_FILE = FIGURES / "spectrum_statistics.csv"

def flat_test(values, errors):
    values, errors = np.asarray(values, float), np.asarray(errors, float)
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    values, errors = values[good], errors[good]
    weights = 1 / errors**2
    mean = np.sum(weights * values) / np.sum(weights)
    statistic = np.sum(((values - mean) / errors)**2)
    dof = len(values) - 1
    return {"n": len(values), "mean": mean, "chi2": statistic, "dof": dof,
            "p": chi2.sf(statistic, dof)}

def offset_model_test(wavelength, values, errors, model_wavelength, model_values):
    model = np.interp(wavelength, model_wavelength, model_values)
    weights = 1 / errors**2
    offset = np.sum(weights * (values - model)) / np.sum(weights)
    statistic = np.sum(((values - model - offset) / errors)**2)
    dof = len(values) - 1
    return {"n": len(values), "offset": offset, "chi2": statistic,
            "dof": dof, "p": chi2.sf(statistic, dof), "model": model + offset}

def write_rows(rows):
    fields = sorted({key for row in rows for key in row})
    with STATS_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)

import h5py

FIGURE_FILE = FIGURES / "wasp43b_published_spectrum.png"

def linear_test(wavelength, values, errors):
    x = wavelength - np.mean(wavelength); design = np.column_stack([np.ones(len(x)), x])
    weights = 1 / errors**2; covariance = np.linalg.inv(design.T @ (weights[:, None] * design))
    coefficients = covariance @ (design.T @ (weights * values)); model = design @ coefficients
    statistic = np.sum(((values - model) / errors)**2); dof = len(values) - 2
    return {"slope_ppm_per_micron": coefficients[1], "slope_error": np.sqrt(covariance[1, 1]),
            "linear_chi2": statistic, "linear_dof": dof, "linear_p": chi2.sf(statistic, dof)}

def main():
    FIGURES.mkdir(exist_ok=True)
    with h5py.File(DATA / "fiducial_combined.h5") as handle:
        phases = handle["phase"][...]; wavelength = handle["wavelength"][...]
        values = handle["fp_fs"][...]; errors = handle["fp_fs_error"][...]
    rows = []; fig, ax = plt.subplots(figsize=(9.2, 5.3))
    for phase, spectrum, uncertainty in zip(phases, values, errors):
        flat = flat_test(spectrum, uncertainty); linear = linear_test(wavelength, spectrum, uncertainty)
        rows.append({"comparison": f"phase {phase:.2f}", **flat, **linear})
        ax.errorbar(wavelength, spectrum, yerr=uncertainty, fmt="o-", ms=3.5, lw=1.2, label=f"orbital phase {phase:.2f}")
    write_rows(rows)
    ax.set(xlabel="Wavelength [micron]", ylabel="Planet/star flux ratio [ppm]",
           title="WASP-43 b: published JWST MIRI/LRS phase-resolved emission spectra")
    ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=8); fig.tight_layout()
    fig.savefig(FIGURE_FILE, dpi=190); plt.close(fig)
    return {"rows": rows, "n": len(wavelength)}

if __name__ == "__main__":
    result = main(); print(f"WASP-43 b: four orbital phases, {result['n']} wavelength bins each")
