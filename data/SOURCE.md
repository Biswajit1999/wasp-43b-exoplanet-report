# Data sources

## TESS light curve

- File: `tess2019058134432-s0009-0000000036734222-0139-s_lc.fits`
- Archive: Mikulski Archive for Space Telescopes (MAST), TESS SPOC light-curve product
- TESS sector: 9
- TIC target ID: 36734222
- MAST observation ID: 62892339
- MAST data URI: `mast:TESS/product/tess2019058134432-s0009-0000000036734222-0139-s_lc.fits`
- Exact download URL: <https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS%2Fproduct%2Ftess2019058134432-s0009-0000000036734222-0139-s_lc.fits>
- Collection DOI: [10.17909/t9-nmc8-f686](https://doi.org/10.17909/t9-nmc8-f686) (TESS 2-minute light curves, all sectors; sector 9 used here)
- Retrieved: 2026-08-15
- SHA-256: `24e025d8efb700db1b83981d7d40d5145b116f06f3b387822e2f4af1fe136374`

The FITS file is stored unmodified. The analysis reads `TIME`, `PDCSAP_FLUX`,
`PDCSAP_FLUX_ERR`, and `QUALITY`. PDCSAP flux is the SPOC light curve with common
instrumental trends removed and aperture/crowding corrections applied; this does
not make it free of residual stellar or instrumental systematics.

## System parameters

- File: `system_parameters.csv`
- Service: NASA Exoplanet Archive TAP, `pscomppars` table
- Exact query: <https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name%2Chostname%2Cra%2Cdec%2Cpl_orbper%2Cpl_tranmid%2Cpl_trandur%2Cpl_rade%2Cpl_bmasse%2Cpl_eqt%2Cpl_orbsmax%2Csy_dist%2Csy_tmag%2Cst_teff%2Cst_rad%2Cst_mass%2Cdisc_year%2Cdiscoverymethod%2Cdisc_refname%2Cdisc_pubdate%2Cdisc_facility+from+pscomppars+where+pl_name%3D%27WASP-43+b%27&format=csv>
- Retrieved: 2026-08-15

The saved row is the input actually used by `scripts/analyze_transit.py`; the
analysis does not query a changing live service at run time.


## Additional TESS sectors for robustness analysis

All are unmodified standard-cadence SPOC light curves from the same [MAST TESS collection](https://doi.org/10.17909/t9-nmc8-f686).

- Sector 9: `tess2019058134432-s0009-0000000036734222-0139-s_lc.fits` (1,848,960 bytes)
  - MAST URI: `mast:TESS/product/tess2019058134432-s0009-0000000036734222-0139-s_lc.fits`
  - SHA-256: `24e025d8efb700db1b83981d7d40d5145b116f06f3b387822e2f4af1fe136374`

## Published planetary spectrum

- Archive record: [10.5281/zenodo.10525170](https://zenodo.org/records/10525170)
- Archive file: `WASP43b_MIRI_Data.zip` (11,533,870 bytes)
- Archive checksum reported by Zenodo and independently reproduced: MD5 `d14c633b7cdc15a9c990a92cc4fc9b86`
- Data type: phase-resolved emission; instrument: JWST MIRI/LRS
- `data/spectra/fiducial_combined.h5` is byte-identical to archive member `WASP43b_MIRI_Data/2_Planetary_Spectra/fiducial_combined.h5` — 9,674 bytes; MD5 `3c0de266fee1fccc832d87954efa7ba0`; SHA-256 `3c270fd46236ec671b88c35846a95c11aef3c2030b61434b600c74d844b4a2a1`
- `data/spectra/eureka_v1.h5` is byte-identical to archive member `WASP43b_MIRI_Data/2_Planetary_Spectra/eureka_v1.h5` — 12,431 bytes; MD5 `f4b242e27f21ca6c5377a2fa71926588`; SHA-256 `c3a05fbb77638f0101986c40fbbe34e6421262069fbdbabf956728daa7cb4ba2`

The complete machine-readable mapping is in `data/zenodo_manifest.csv`.
`python scripts/verify_zenodo_provenance.py` checks the committed files offline;
pass a downloaded archive with `--archive` to verify the outer ZIP and exact
member bytes as well.  The audit was performed on 2026-09-25.

The source analysis restricts its final spectroscopic interpretation to
5–10.5 microns: the 10.6–11.8 micron shadowed-region data could not be
detrended reliably.  The repository therefore retains all 14 archived bins
for provenance and display but uses only the 11 bins centred at 5.25–10.25
microns in inferential spectrum tests and thermal summaries.
