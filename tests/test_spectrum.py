"""Reproduction checks for the archived published spectrum."""
from pathlib import Path
import numpy as np
import analyze_spectrum as spectrum

def test_spectrum_analysis_is_finite_and_reproducible():
    result = spectrum.main()
    assert result["n"] == 14
    assert result["n_tested"] == 11
    assert result["n_excluded"] == 3
    assert result["rows"]
    assert all(row["n"] == 11 for row in result["rows"])
    assert all(np.isfinite(row["chi2"]) and row["dof"] > 0 for row in result["rows"])
    for path in (spectrum.STATS_FILE, spectrum.FIGURE_FILE):
        assert Path(path).is_file() and Path(path).stat().st_size > 100
