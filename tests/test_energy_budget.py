from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_energy_budget as energy


def test_blackbody_inversion_is_monotonic():
    wavelength = np.asarray([5.0, 5.0, 5.0])
    temperature = energy.brightness_temperature(wavelength, np.asarray([500, 1000, 2000]), 4400, 0.025)
    assert np.all(np.diff(temperature) > 0)


def test_sinusoid_recovers_peak_phase():
    phase = np.asarray([0.0, 0.25, 0.5, 0.75])
    truth = 4000 - 2000 * np.cos(2 * np.pi * phase) + 300 * np.sin(2 * np.pi * phase)
    result = energy.sinusoid_fit(phase, truth, np.full(4, 20.0))
    assert -15 < result["orbital_offset_deg"] < 0
    assert result["semi_amplitude_ppm"] > 1900


def test_real_energy_budget_matches_expected_physical_regime():
    result = energy.main()
    assert 1500 < result["day_temperature"] < 1700
    assert 800 < result["night_temperature"] < 1000
    assert 600 < result["temperature_contrast"] < 850
    assert 2.5 < result["day_night_flux_ratio"] < 4.0
    assert 0 <= result["grid"]["best_albedo"] <= 0.8
    assert 0 < result["grid"]["best_redistribution"] <= 1
    offset = result["analysis"]["fiducial"]["sinusoid"]["orbital_offset_deg"]
    assert -20 < offset < 0
    for path in (energy.SUMMARY_FILE, energy.TEMPERATURE_FILE, energy.FIGURE_FILE):
        assert Path(path).is_file() and Path(path).stat().st_size > 200
