"""Sensitivity and identifiability audit for the four archived phase spectra.

This module deliberately limits its claim.  It asks whether a pre-eclipse
maximum is stable under reasonable analysis choices.  It does not reproduce
the source paper's full time-series, two-harmonic phase-curve inference.
"""

from __future__ import annotations

import csv
from itertools import product
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spectra"
FIGURES = ROOT / "figures"
MULTIVERSE_FILE = FIGURES / "phase_offset_multiverse.csv"
JACKKNIFE_FILE = FIGURES / "phase_offset_wavelength_jackknife.csv"
SUMMARY_FILE = FIGURES / "phase_identifiability_summary.csv"
FIGURE_FILE = FIGURES / "phase_identifiability_audit.png"
TRUSTED_MAXIMUM_MICRON = 10.5


def load_reduction(filename: str, error_rule: str = "mean") -> dict[str, np.ndarray]:
    with h5py.File(DATA / filename) as handle:
        wavelength = np.asarray(handle["wavelength"][...], dtype=float)
        phase = np.asarray(handle["phase"][...], dtype=float)
        flux = np.asarray(handle["fp_fs"][...], dtype=float)
        if "fp_fs_error" in handle:
            negative = positive = np.asarray(handle["fp_fs_error"][...], dtype=float)
        else:
            negative = np.asarray(handle["fp_fs_errorNeg"][...], dtype=float)
            positive = np.asarray(handle["fp_fs_errorPos"][...], dtype=float)
    if error_rule == "mean":
        error = 0.5 * (negative + positive)
    elif error_rule == "maximum":
        error = np.maximum(negative, positive)
    else:
        raise ValueError(f"unknown error rule: {error_rule}")
    return {"wavelength": wavelength, "phase": phase, "flux": flux, "error": error}


def aggregation_weights(wavelength: np.ndarray, error: np.ndarray, rule: str) -> np.ndarray:
    if rule == "uniform":
        raw = np.ones(len(wavelength))
    elif rule == "inverse_variance":
        raw = 1.0 / error**2
    elif rule == "trapezoid":
        raw = np.zeros(len(wavelength))
        spacing = np.diff(wavelength)
        raw[:-1] += 0.5 * spacing
        raw[1:] += 0.5 * spacing
    else:
        raise ValueError(f"unknown aggregation rule: {rule}")
    return raw / raw.sum()


def aggregate(reduction: dict[str, np.ndarray], selected: np.ndarray, rule: str) -> tuple[np.ndarray, np.ndarray]:
    wave = reduction["wavelength"][selected]
    values, uncertainties = [], []
    for flux, error in zip(reduction["flux"][:, selected], reduction["error"][:, selected]):
        weights = aggregation_weights(wave, error, rule)
        values.append(float(weights @ flux))
        uncertainties.append(float(np.sqrt(np.sum((weights * error) ** 2))))
    return np.asarray(values), np.asarray(uncertainties)


def harmonic_design(phase: np.ndarray, harmonics: int) -> np.ndarray:
    columns = [np.ones(len(phase))]
    for harmonic in range(1, harmonics + 1):
        columns.extend([
            np.cos(2 * np.pi * harmonic * phase),
            np.sin(2 * np.pi * harmonic * phase),
        ])
    return np.column_stack(columns)


def first_harmonic_fit(phase: np.ndarray, values: np.ndarray, errors: np.ndarray) -> dict[str, float]:
    design = harmonic_design(phase, 1)
    precision = 1.0 / errors**2
    normal = design.T @ (precision[:, None] * design)
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ (design.T @ (precision * values))
    model = design @ coefficients
    phase_max = float(np.arctan2(coefficients[2], coefficients[1]) / (2 * np.pi) % 1.0)
    offset = float(((phase_max - 0.5 + 0.5) % 1.0 - 0.5) * 360.0)
    rng = np.random.default_rng(4300)
    draws = rng.multivariate_normal(coefficients, covariance, size=20_000)
    draw_phase = np.arctan2(draws[:, 2], draws[:, 1]) / (2 * np.pi) % 1.0
    draw_offset = ((draw_phase - 0.5 + 0.5) % 1.0 - 0.5) * 360.0
    return {
        "offset_deg": offset,
        "formal_error_deg": float(0.5 * (np.percentile(draw_offset, 84) - np.percentile(draw_offset, 16))),
        "semi_amplitude_ppm": float(np.hypot(coefficients[1], coefficients[2])),
        "chi_square": float(np.sum(((values - model) / errors) ** 2)),
        "dof": int(len(values) - design.shape[1]),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> dict[str, object]:
    FIGURES.mkdir(exist_ok=True)
    files = {"fiducial": "fiducial_combined.h5", "eureka_v1": "eureka_v1.h5"}
    lower_bounds = (5.25, 5.75, 6.25)
    upper_bounds = (8.25, 8.75, 9.25, 9.75, 10.25)
    aggregation_rules = ("trapezoid", "uniform", "inverse_variance")
    error_rules = ("mean", "maximum")
    multiverse: list[dict[str, object]] = []
    for reduction_name, lower, upper, aggregation_rule, error_rule in product(
        files, lower_bounds, upper_bounds, aggregation_rules, error_rules
    ):
        reduction = load_reduction(files[reduction_name], error_rule)
        selected = (reduction["wavelength"] >= lower) & (reduction["wavelength"] <= upper)
        if selected.sum() < 5:
            continue
        values, errors = aggregate(reduction, selected, aggregation_rule)
        fit = first_harmonic_fit(reduction["phase"], values, errors)
        multiverse.append({
            "reduction": reduction_name,
            "wavelength_min_micron": lower,
            "wavelength_max_micron": upper,
            "n_wavelength_bins": int(selected.sum()),
            "aggregation_rule": aggregation_rule,
            "error_rule": error_rule,
            **fit,
        })
    write_csv(MULTIVERSE_FILE, multiverse)

    jackknife: list[dict[str, object]] = []
    for reduction_name, filename in files.items():
        reduction = load_reduction(filename, "mean")
        trusted = reduction["wavelength"] <= TRUSTED_MAXIMUM_MICRON
        for omitted in reduction["wavelength"][trusted]:
            selected = trusted & (reduction["wavelength"] != omitted)
            values, errors = aggregate(reduction, selected, "trapezoid")
            jackknife.append({
                "reduction": reduction_name,
                "omitted_wavelength_micron": float(omitted),
                **first_harmonic_fit(reduction["phase"], values, errors),
            })
    write_csv(JACKKNIFE_FILE, jackknife)

    offsets = np.asarray([row["offset_deg"] for row in multiverse], dtype=float)
    fiducial_offsets = np.asarray([
        row["offset_deg"] for row in multiverse if row["reduction"] == "fiducial"
    ], dtype=float)
    eureka_offsets = np.asarray([
        row["offset_deg"] for row in multiverse if row["reduction"] == "eureka_v1"
    ], dtype=float)
    jackknife_offsets = np.asarray([row["offset_deg"] for row in jackknife], dtype=float)
    phase = load_reduction(files["fiducial"])["phase"]
    one_harmonic = harmonic_design(phase, 1)
    two_harmonic = harmonic_design(phase, 2)
    summary = [
        {"quantity": "multiverse_designs", "value": len(multiverse), "interpretation": "predeclared analysis combinations"},
        {"quantity": "offset_min_deg", "value": float(offsets.min()), "interpretation": "negative means before eclipse"},
        {"quantity": "offset_median_deg", "value": float(np.median(offsets)), "interpretation": "across all designs"},
        {"quantity": "offset_max_deg", "value": float(offsets.max()), "interpretation": "negative means before eclipse"},
        {"quantity": "all_designs_pre_eclipse", "value": bool(np.all(offsets < 0)), "interpretation": "directional robustness only"},
        {"quantity": "fiducial_offset_range_deg", "value": f"{fiducial_offsets.min():.6f} to {fiducial_offsets.max():.6f}", "interpretation": "analysis-choice envelope"},
        {"quantity": "eureka_offset_range_deg", "value": f"{eureka_offsets.min():.6f} to {eureka_offsets.max():.6f}", "interpretation": "analysis-choice envelope"},
        {"quantity": "jackknife_offset_range_deg", "value": f"{jackknife_offsets.min():.6f} to {jackknife_offsets.max():.6f}", "interpretation": "leave-one-wavelength-out, trusted bins"},
        {"quantity": "phase_observations", "value": len(phase), "interpretation": "archive phase bins"},
        {"quantity": "first_harmonic_parameters", "value": one_harmonic.shape[1], "interpretation": "one residual degree of freedom"},
        {"quantity": "first_harmonic_rank", "value": int(np.linalg.matrix_rank(one_harmonic)), "interpretation": "identified under the assumed model"},
        {"quantity": "two_harmonic_parameters", "value": two_harmonic.shape[1], "interpretation": "more parameters than observations"},
        {"quantity": "two_harmonic_rank", "value": int(np.linalg.matrix_rank(two_harmonic)), "interpretation": "rank deficient; source model cannot be reproduced"},
    ]
    write_csv(SUMMARY_FILE, summary)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    ax = axes[0]
    for name, color in (("fiducial", "#0369a1"), ("eureka_v1", "#b45309")):
        subset = np.asarray([row["offset_deg"] for row in multiverse if row["reduction"] == name])
        ax.hist(subset, bins=np.linspace(offsets.min() - 0.2, offsets.max() + 0.2, 20), alpha=0.62, label=name, color=color)
    ax.axvline(0, color="#111827", lw=1.3)
    ax.set(xlabel="First-harmonic maximum relative to eclipse [deg]", ylabel="Analysis designs", title="180-choice sensitivity multiverse")
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.2)

    ax = axes[1]
    for index, name in enumerate(files):
        rows = [row for row in jackknife if row["reduction"] == name]
        ax.plot([row["omitted_wavelength_micron"] for row in rows], [row["offset_deg"] for row in rows], "o-", label=name, color=("#0369a1", "#b45309")[index])
    ax.axhline(0, color="#111827", lw=1.3)
    ax.set(xlabel="Omitted wavelength [micron]", ylabel="Maximum relative to eclipse [deg]", title="Leave-one-wavelength-out audit")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.suptitle("WASP-43 b: a pre-eclipse maximum is robust; its precision is not", fontsize=15, weight="bold")
    fig.savefig(FIGURE_FILE, dpi=190)
    plt.close(fig)
    return {"multiverse": multiverse, "jackknife": jackknife, "summary": summary}


if __name__ == "__main__":
    result = main()
    offsets = np.asarray([row["offset_deg"] for row in result["multiverse"]])
    print(f"{len(offsets)} designs: offset {offsets.min():.2f} to {offsets.max():.2f} deg; all pre-eclipse={bool(np.all(offsets < 0))}")
