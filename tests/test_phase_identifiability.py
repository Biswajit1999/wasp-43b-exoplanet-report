from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_phase_identifiability as audit
import verify_zenodo_provenance as provenance


def test_committed_spectra_match_audited_archive_manifest():
    checked = provenance.verify()
    assert {row["local_path"] for row in checked} == {
        "data/spectra/fiducial_combined.h5",
        "data/spectra/eureka_v1.h5",
    }


def test_two_harmonics_are_underidentified_by_four_phase_bins():
    phase = audit.load_reduction("fiducial_combined.h5")["phase"]
    design = audit.harmonic_design(phase, 2)
    assert design.shape == (4, 5)
    assert np.linalg.matrix_rank(design) == 4


def test_pre_eclipse_direction_survives_predeclared_multiverse():
    result = audit.main()
    offsets = np.asarray([row["offset_deg"] for row in result["multiverse"]])
    assert len(offsets) == 180
    assert np.all(offsets < 0)
    assert offsets.max() - offsets.min() > 1.5
    assert offsets.max() - offsets.min() < 3.0
    for path in (
        audit.MULTIVERSE_FILE,
        audit.JACKKNIFE_FILE,
        audit.SUMMARY_FILE,
        audit.FIGURE_FILE,
    ):
        assert path.is_file() and path.stat().st_size > 200
