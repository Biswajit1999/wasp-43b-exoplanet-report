"""Verify the committed spectra against the audited Zenodo archive manifest.

The default check is offline and verifies the exact byte hashes recorded while
auditing Zenodo record 10525170.  Supplying ``--archive`` additionally checks
the outer ZIP checksum and compares each local file with its archive member.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "zenodo_manifest.csv"


def digest(payload: bytes, algorithm: str) -> str:
    return hashlib.new(algorithm, payload).hexdigest()


def verify(archive: Path | None = None) -> list[dict[str, str]]:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Zenodo manifest is empty")

    archive_file = None
    if archive is not None:
        payload = archive.read_bytes()
        expected = rows[0]["archive_md5"]
        if digest(payload, "md5") != expected:
            raise ValueError("archive MD5 does not match Zenodo's checksum")
        archive_file = zipfile.ZipFile(archive)

    checked = []
    try:
        for row in rows:
            local = ROOT / row["local_path"]
            payload = local.read_bytes()
            if len(payload) != int(row["size_bytes"]):
                raise ValueError(f"size mismatch: {row['local_path']}")
            for algorithm in ("md5", "sha256"):
                expected = row[f"entry_{algorithm}"]
                if digest(payload, algorithm) != expected:
                    raise ValueError(f"{algorithm} mismatch: {row['local_path']}")
            if archive_file is not None and archive_file.read(row["archive_entry"]) != payload:
                raise ValueError(f"archive member differs: {row['archive_entry']}")
            checked.append(row)
    finally:
        if archive_file is not None:
            archive_file.close()
    return checked


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, help="optional downloaded WASP43b_MIRI_Data.zip")
    args = parser.parse_args()
    checked = verify(args.archive)
    mode = "archive + local bytes" if args.archive else "offline manifest hashes"
    print(f"Verified {len(checked)} Zenodo spectrum files ({mode}).")


if __name__ == "__main__":
    main()
