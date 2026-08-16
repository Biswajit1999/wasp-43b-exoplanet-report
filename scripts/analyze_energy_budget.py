"""Band-limited thermal and circulation diagnostics for WASP-43 b.

The public four-phase MIRI/LRS spectra are converted to blackbody-star colour
temperatures, integrated into a 5.25-10.25 micron phase curve, and compared
between the fiducial and Eureka reductions.  A Cowan & Agol (2011) albedo-
recirculation grid is shown only as an illustrative mapping because band
brightness temperatures are not bolometric effective temperatures.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, h, k


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spectra"
FIGURES = ROOT / "figures"
SUMMARY_FILE = FIGURES / "energy_budget_statistics.csv"
TEMPERATURE_FILE = FIGURES / "phase_brightness_temperatures.csv"
FIGURE_FILE = FIGURES / "wasp43b_energy_budget.png"

R_EARTH_M = 6.371e6
R_SUN_M = 6.957e8
R_SUN_AU = 0.00465047


def system_parameters() -> dict[str, float]:
    with (ROOT / "data" / "system_parameters.csv").open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    return {
        "planet_radius_earth": float(row["pl_rade"]),
        "star_radius_sun": float(row["st_rad"]),
        "star_temperature_k": float(row["st_teff"]),
        "semimajor_axis_au": float(row["pl_orbsmax"]),
    }


def load_reduction(filename: str) -> dict[str, np.ndarray]:
    with h5py.File(DATA / filename) as handle:
        wavelength = np.asarray(handle["wavelength"][...], dtype=float)
        phase = np.asarray(handle["phase"][...], dtype=float)
        flux = np.asarray(handle["fp_fs"][...], dtype=float)
        if "fp_fs_error" in handle:
            error = np.asarray(handle["fp_fs_error"][...], dtype=float)
        else:
            error = 0.5 * (
                np.asarray(handle["fp_fs_errorNeg"][...], dtype=float)
                + np.asarray(handle["fp_fs_errorPos"][...], dtype=float)
            )
    return {"wavelength": wavelength, "phase": phase, "flux": flux, "error": error}


def brightness_temperature(
    wavelength_micron: np.ndarray,
    flux_ratio_ppm: np.ndarray,
    star_temperature_k: float,
    area_ratio: float,
) -> np.ndarray:
    """Invert a blackbody planet/star monochromatic flux ratio."""
    wavelength_m = np.asarray(wavelength_micron, dtype=float) * 1e-6
    flux_ratio = np.maximum(np.asarray(flux_ratio_ppm, dtype=float) * 1e-6, 1e-14)
    stellar_exponent = h * c / (wavelength_m * k * star_temperature_k)
    planet_exponent = np.log1p(np.expm1(stellar_exponent) / (flux_ratio / area_ratio))
    return h * c / (wavelength_m * k * planet_exponent)


def temperature_spectrum(
    reduction: dict[str, np.ndarray],
    parameters: dict[str, float],
) -> tuple[np.ndarray, np.ndarray]:
    area_ratio = (
        parameters["planet_radius_earth"] * R_EARTH_M
        / (parameters["star_radius_sun"] * R_SUN_M)
    ) ** 2
    temperatures = []
    errors = []
    for values, uncertainty in zip(reduction["flux"], reduction["error"]):
        central = brightness_temperature(
            reduction["wavelength"], values,
            parameters["star_temperature_k"], area_ratio,
        )
        upper = brightness_temperature(
            reduction["wavelength"], values + uncertainty,
            parameters["star_temperature_k"], area_ratio,
        )
        lower = brightness_temperature(
            reduction["wavelength"], np.maximum(values - uncertainty, 1e-6),
            parameters["star_temperature_k"], area_ratio,
        )
        temperatures.append(central)
        errors.append(0.5 * (upper - lower))
    return np.asarray(temperatures), np.asarray(errors)


def weighted_temperature(values: np.ndarray, errors: np.ndarray) -> dict[str, float]:
    weights = 1.0 / errors**2
    mean = float(np.sum(weights * values) / np.sum(weights))
    formal_error = float(np.sqrt(1.0 / np.sum(weights)))
    chi_square = float(np.sum(((values - mean) / errors) ** 2))
    dof = len(values) - 1
    # Spectral structure makes the constant-temperature scatter much larger
    # than the propagated per-bin error.  Preserve it in the summary error.
    scatter_scale = np.sqrt(max(chi_square / dof, 1.0))
    return {
        "mean": mean,
        "formal_error": formal_error,
        "error": formal_error * scatter_scale,
        "chi_square": chi_square,
        "dof": dof,
        "scatter_scale": scatter_scale,
    }


def integration_weights(wavelength: np.ndarray) -> np.ndarray:
    x = np.asarray(wavelength, dtype=float)
    weights = np.zeros(len(x))
    spacing = np.diff(x)
    weights[:-1] += 0.5 * spacing
    weights[1:] += 0.5 * spacing
    return weights / (x[-1] - x[0])


def band_average(reduction: dict[str, np.ndarray], maximum_wavelength: float = 10.5) -> dict[str, np.ndarray]:
    selected = reduction["wavelength"] <= maximum_wavelength
    wavelength = reduction["wavelength"][selected]
    weights = integration_weights(wavelength)
    values = reduction["flux"][:, selected] @ weights
    errors = np.sqrt(np.sum((reduction["error"][:, selected] * weights) ** 2, axis=1))
    return {"phase": reduction["phase"], "value": values, "error": errors, "weights": weights}


def sinusoid_fit(phase: np.ndarray, values: np.ndarray, errors: np.ndarray) -> dict[str, object]:
    design = np.column_stack([
        np.ones(len(phase)), np.cos(2 * np.pi * phase), np.sin(2 * np.pi * phase)
    ])
    precision = 1.0 / errors**2
    covariance = np.linalg.pinv(design.T @ (precision[:, None] * design))
    coefficients = covariance @ (design.T @ (precision * values))
    model = design @ coefficients
    chi_square = float(np.sum(((values - model) / errors) ** 2))
    phase_max = float(np.arctan2(coefficients[2], coefficients[1]) / (2 * np.pi) % 1.0)
    orbital_offset = float(((phase_max - 0.5 + 0.5) % 1.0 - 0.5) * 360.0)

    rng = np.random.default_rng(43)
    draws = rng.multivariate_normal(coefficients, covariance, size=50_000)
    draw_phase = np.arctan2(draws[:, 2], draws[:, 1]) / (2 * np.pi) % 1.0
    draw_offset = ((draw_phase - 0.5 + 0.5) % 1.0 - 0.5) * 360.0
    draw_amplitude = np.hypot(draws[:, 1], draws[:, 2])
    return {
        "coefficients": coefficients,
        "covariance": covariance,
        "model": model,
        "chi_square": chi_square,
        "dof": len(values) - len(coefficients),
        "phase_max": phase_max,
        "orbital_offset_deg": orbital_offset,
        "orbital_offset_error_deg": 0.5 * (np.percentile(draw_offset, 84) - np.percentile(draw_offset, 16)),
        "eastward_hotspot_proxy_deg": -orbital_offset,
        "semi_amplitude_ppm": float(np.hypot(coefficients[1], coefficients[2])),
        "semi_amplitude_error_ppm": 0.5 * (np.percentile(draw_amplitude, 84) - np.percentile(draw_amplitude, 16)),
    }


def energy_grid(
    day_temperature: float,
    day_error: float,
    night_temperature: float,
    night_error: float,
    substellar_temperature: float,
) -> dict[str, object]:
    albedo = np.linspace(0.0, 0.8, 401)
    redistribution = np.linspace(0.001, 1.0, 400)
    aa, ee = np.meshgrid(albedo, redistribution)
    common = substellar_temperature * np.power(1.0 - aa, 0.25)
    model_day = common * np.power(2.0 / 3.0 - 5.0 * ee / 12.0, 0.25)
    model_night = common * np.power(ee / 4.0, 0.25)
    chi_square = (
        ((day_temperature - model_day) / day_error) ** 2
        + ((night_temperature - model_night) / night_error) ** 2
    )
    best_index = np.unravel_index(np.argmin(chi_square), chi_square.shape)
    minimum = float(chi_square[best_index])
    region = chi_square <= minimum + 2.30
    return {
        "albedo": albedo,
        "redistribution": redistribution,
        "chi_square": chi_square,
        "minimum_chi_square": minimum,
        "best_albedo": float(aa[best_index]),
        "best_redistribution": float(ee[best_index]),
        "albedo_68_min": float(aa[region].min()),
        "albedo_68_max": float(aa[region].max()),
        "redistribution_68_min": float(ee[region].min()),
        "redistribution_68_max": float(ee[region].max()),
    }


def main() -> dict[str, object]:
    FIGURES.mkdir(exist_ok=True)
    parameters = system_parameters()
    reductions = {
        "fiducial": load_reduction("fiducial_combined.h5"),
        "eureka_v1": load_reduction("eureka_v1.h5"),
    }
    analysis: dict[str, dict[str, object]] = {}
    temperature_rows = []
    for name, reduction in reductions.items():
        temperatures, temperature_errors = temperature_spectrum(reduction, parameters)
        valid = reduction["wavelength"] <= 10.5
        summaries = [
            weighted_temperature(values[valid], errors[valid])
            for values, errors in zip(temperatures, temperature_errors)
        ]
        band = band_average(reduction)
        sinusoid = sinusoid_fit(band["phase"], band["value"], band["error"])
        analysis[name] = {
            "temperatures": temperatures,
            "temperature_errors": temperature_errors,
            "summaries": summaries,
            "band": band,
            "sinusoid": sinusoid,
        }
        for phase_index, phase in enumerate(reduction["phase"]):
            for wavelength_index, wavelength in enumerate(reduction["wavelength"]):
                temperature_rows.append({
                    "reduction": name,
                    "phase": phase,
                    "wavelength_micron": wavelength,
                    "planet_star_flux_ppm": reduction["flux"][phase_index, wavelength_index],
                    "flux_error_ppm": reduction["error"][phase_index, wavelength_index],
                    "blackbody_star_brightness_temperature_k": temperatures[phase_index, wavelength_index],
                    "temperature_error_k": temperature_errors[phase_index, wavelength_index],
                    "used_in_5p25_10p25_summary": bool(wavelength <= 10.5),
                })

    star_radius_au = parameters["star_radius_sun"] * R_SUN_AU
    substellar_temperature = parameters["star_temperature_k"] * np.sqrt(
        star_radius_au / parameters["semimajor_axis_au"]
    )
    fiducial_summaries = analysis["fiducial"]["summaries"]
    night = fiducial_summaries[0]
    day = fiducial_summaries[2]
    grid = energy_grid(
        day["mean"], day["error"], night["mean"], night["error"], substellar_temperature
    )
    fiducial_band = analysis["fiducial"]["band"]
    day_night_flux_ratio = float(fiducial_band["value"][2] / fiducial_band["value"][0])
    temperature_contrast = float(day["mean"] - night["mean"])
    temperature_contrast_error = float(np.hypot(day["error"], night["error"]))

    summary_rows: list[tuple[str, object, str]] = [
        ("analysis_band_min", 5.25, "micron"),
        ("analysis_band_max", 10.25, "micron; bins above 10.5 excluded"),
        ("substellar_irradiation_temperature", substellar_temperature, "K; saved stellar parameters"),
        ("fiducial_nightside_brightness_temperature", night["mean"], "K; blackbody-star approximation"),
        ("fiducial_nightside_temperature_error", night["error"], "K; inflated by spectral scatter"),
        ("fiducial_dayside_brightness_temperature", day["mean"], "K; blackbody-star approximation"),
        ("fiducial_dayside_temperature_error", day["error"], "K; inflated by spectral scatter"),
        ("fiducial_day_night_temperature_contrast", temperature_contrast, "K"),
        ("fiducial_day_night_temperature_contrast_error", temperature_contrast_error, "K"),
        ("fiducial_day_night_band_flux_ratio", day_night_flux_ratio, "dimensionless"),
        ("fiducial_phase_max", analysis["fiducial"]["sinusoid"]["phase_max"], "orbital phase"),
        ("fiducial_peak_offset_from_eclipse", analysis["fiducial"]["sinusoid"]["orbital_offset_deg"], "degrees; negative is before eclipse"),
        ("fiducial_peak_offset_error", analysis["fiducial"]["sinusoid"]["orbital_offset_error_deg"], "degrees; formal four-point sinusoid"),
        ("fiducial_eastward_hotspot_proxy", analysis["fiducial"]["sinusoid"]["eastward_hotspot_proxy_deg"], "degrees; interpretation assumes longitudinal thermal dominance"),
        ("illustrative_best_bond_albedo", grid["best_albedo"], "band-temperature Cowan-Agol mapping"),
        ("illustrative_bond_albedo_68_min", grid["albedo_68_min"], "conditional"),
        ("illustrative_bond_albedo_68_max", grid["albedo_68_max"], "conditional"),
        ("illustrative_best_redistribution", grid["best_redistribution"], "epsilon; 0 none, 1 complete"),
        ("illustrative_redistribution_68_min", grid["redistribution_68_min"], "conditional"),
        ("illustrative_redistribution_68_max", grid["redistribution_68_max"], "conditional"),
    ]
    for name in ("fiducial", "eureka_v1"):
        summary = analysis[name]["summaries"]
        sinusoid = analysis[name]["sinusoid"]
        summary_rows.extend([
            (f"{name}_night_temperature", summary[0]["mean"], "K"),
            (f"{name}_day_temperature", summary[2]["mean"], "K"),
            (f"{name}_peak_offset_from_eclipse", sinusoid["orbital_offset_deg"], "degrees"),
            (f"{name}_phase_curve_semi_amplitude", sinusoid["semi_amplitude_ppm"], "ppm"),
        ])
    with SUMMARY_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        for name, value, unit in summary_rows:
            writer.writerow([name, f"{value:.12g}" if isinstance(value, float) else value, unit])
    with TEMPERATURE_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(temperature_rows[0]))
        writer.writeheader()
        writer.writerows(temperature_rows)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    ax = axes[0, 0]
    fiducial = reductions["fiducial"]
    for index, phase in enumerate(fiducial["phase"]):
        ax.errorbar(
            fiducial["wavelength"], analysis["fiducial"]["temperatures"][index],
            yerr=analysis["fiducial"]["temperature_errors"][index],
            fmt="o-", ms=3.2, lw=1.2, label=f"phase {phase:.2f}",
        )
    ax.axvspan(10.5, fiducial["wavelength"].max() + 0.1, color="#94a3b8", alpha=0.16, label=">10.5 micron excluded from summaries")
    ax.set(xlabel="Wavelength [micron]", ylabel="Brightness temperature [K]", title="Blackbody-star colour temperatures")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    phase_grid = np.linspace(0, 1, 400)
    colors = {"fiducial": "#0f766e", "eureka_v1": "#7c3aed"}
    for name in ("fiducial", "eureka_v1"):
        band = analysis[name]["band"]
        sinusoid = analysis[name]["sinusoid"]
        coefficients = sinusoid["coefficients"]
        curve = coefficients[0] + coefficients[1] * np.cos(2 * np.pi * phase_grid) + coefficients[2] * np.sin(2 * np.pi * phase_grid)
        ax.errorbar(band["phase"], band["value"], yerr=band["error"], fmt="o", color=colors[name], label=f"{name} bins")
        ax.plot(phase_grid, curve, color=colors[name], lw=2, label=f"{name} sinusoid")
    ax.axvline(0.5, color="#475569", ls="--", lw=1.2, label="secondary eclipse")
    ax.set(xlabel="Orbital phase", ylabel="5.25-10.25 micron planet/star flux [ppm]", title="Four-phase broadband reconstruction")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = axes[1, 0]
    delta = grid["chi_square"] - grid["minimum_chi_square"]
    image = ax.contourf(grid["albedo"], grid["redistribution"], np.minimum(delta, 20), levels=np.linspace(0, 20, 21), cmap="viridis_r")
    ax.contour(grid["albedo"], grid["redistribution"], delta, levels=[2.30, 6.18], colors=["white", "#fbbf24"], linewidths=1.6)
    ax.plot(grid["best_albedo"], grid["best_redistribution"], "r*", ms=12, label="illustrative best fit")
    ax.set(xlabel="Bond albedo", ylabel="Heat redistribution epsilon", title="Conditional Cowan-Agol mapping")
    ax.legend(frameon=False, fontsize=8)
    fig.colorbar(image, ax=ax, label="Delta chi-square (clipped at 20)")

    ax = axes[1, 1]
    wavelength = fiducial["wavelength"]
    contrast = analysis["fiducial"]["temperatures"][2] - analysis["fiducial"]["temperatures"][0]
    contrast_error = np.hypot(analysis["fiducial"]["temperature_errors"][2], analysis["fiducial"]["temperature_errors"][0])
    ax.errorbar(wavelength, contrast, yerr=contrast_error, fmt="o-", color="#be123c", ms=4)
    ax.axhline(temperature_contrast, color="#475569", ls="--", label=f"band summary {temperature_contrast:.0f} K")
    ax.axvspan(10.5, wavelength.max() + 0.1, color="#94a3b8", alpha=0.16)
    ax.set(xlabel="Wavelength [micron]", ylabel="Dayside - nightside brightness temperature [K]", title="Wavelength-dependent thermal contrast")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("WASP-43 b: phase-resolved thermal contrast and energy-budget diagnostics", fontsize=15, weight="bold")
    fig.savefig(FIGURE_FILE, dpi=190)
    plt.close(fig)

    return {
        "parameters": parameters,
        "reductions": reductions,
        "analysis": analysis,
        "grid": grid,
        "substellar_temperature": substellar_temperature,
        "day_temperature": day["mean"],
        "day_temperature_error": day["error"],
        "night_temperature": night["mean"],
        "night_temperature_error": night["error"],
        "temperature_contrast": temperature_contrast,
        "temperature_contrast_error": temperature_contrast_error,
        "day_night_flux_ratio": day_night_flux_ratio,
    }


if __name__ == "__main__":
    result = main()
    sinusoid = result["analysis"]["fiducial"]["sinusoid"]
    print(
        f"WASP-43 b: Tday={result['day_temperature']:.0f} +/- {result['day_temperature_error']:.0f} K; "
        f"Tnight={result['night_temperature']:.0f} +/- {result['night_temperature_error']:.0f} K; "
        f"peak={sinusoid['orbital_offset_deg']:+.2f} deg from eclipse"
    )
