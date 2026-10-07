#!/usr/bin/env python3
"""
Export spectral entropy per spectrum from an MGF file to a TSV.

This script reuses the logic from spectral_entropy.py:
- preprocess_peaks
- spectral_entropy
- entropy_quality_label

Usage:
    python export_mgf_entropy_tsv.py \
        --mgf "/path/to/input.mgf" \
        --out "/path/to/output.tsv"

Optional:
    --top 50
    --min-intensity 0.01

Omit --top to keep all peaks after min-intensity filtering.
"""

import argparse
import csv
from pathlib import Path
from typing import Optional

from spectral_entropy import (
    preprocess_peaks,
    spectral_entropy,
    entropy_quality_label,
)


def parse_mgf(mgf_path: Path):
    """
    Parse an MGF file into a list of spectra dicts:
    {
        "scan_number": str,
        "peaks": [[mz, intensity], ...]
    }

    Priority for scan identifier:
    1. SCANS=
    2. TITLE=
    3. sequential spectrum index
    """
    spectra = []

    current_peaks = []
    current_scan = None
    current_title = None
    spectrum_index = 0
    in_ions_block = False

    with mgf_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line == "BEGIN IONS":
                in_ions_block = True
                current_peaks = []
                current_scan = None
                current_title = None
                continue

            if line == "END IONS":
                if in_ions_block:
                    spectrum_index += 1
                    scan_number = current_scan or current_title or str(spectrum_index)
                    spectra.append(
                        {
                            "scan_number": str(scan_number),
                            "peaks": current_peaks,
                        }
                    )
                in_ions_block = False
                continue

            if not in_ions_block:
                continue

            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip().upper()
                value = value.strip()

                if key == "SCANS":
                    current_scan = value
                elif key == "TITLE":
                    current_title = value
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    current_peaks.append([mz, intensity])
                except ValueError:
                    pass

    return spectra


def export_entropy_tsv(
    mgf_path: Path,
    out_path: Path,
    top_n: Optional[int] = None,
    min_intensity_frac: float = 0.01,
):
    spectra = parse_mgf(mgf_path)

    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scan_number",
                "entropy_score",
                "spectra_quality_definition",
            ]
        )

        for spectrum in spectra:
            _, intensities = preprocess_peaks(
                spectrum["peaks"],
                top_n=top_n,
                min_intensity_frac=min_intensity_frac,
            )

            entropy = spectral_entropy(intensities)
            quality_label, _ = entropy_quality_label(entropy)

            writer.writerow(
                [
                    spectrum["scan_number"],
                    f"{entropy:.6f}",
                    quality_label,
                ]
            )

    return len(spectra)


def main():
    parser = argparse.ArgumentParser(
        description="Compute spectral entropy for each spectrum in an MGF and export to TSV."
    )
    parser.add_argument(
        "--mgf",
        required=True,
        help="Path to the input MGF file.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Path to the output TSV file.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=None,
        metavar="N",
        help="Keep only the top N most intense peaks per spectrum. Omit to keep all peaks.",
    )
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.01,
        metavar="FRAC",
        help="Remove peaks below FRAC × base peak intensity (default: 0.01).",
    )

    args = parser.parse_args()

    if args.top is not None and args.top <= 0:
        raise ValueError("--top must be a positive integer when provided.")

    mgf_path = Path(args.mgf)
    out_path = Path(args.out)

    if not mgf_path.exists():
        raise FileNotFoundError(f"MGF file not found: {mgf_path}")

    count = export_entropy_tsv(
        mgf_path=mgf_path,
        out_path=out_path,
        top_n=args.top,
        min_intensity_frac=args.min_intensity,
    )

    print(f"Wrote {count} spectra to {out_path}")


if __name__ == "__main__":
    main()
