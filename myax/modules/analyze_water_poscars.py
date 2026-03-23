#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any


def _load_numpy():
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires NumPy. Activate the environment first, for example: conda activate ase"
        ) from exc

    return np


def _load_ase_read():
    try:
        from ase.io import read
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires ASE. Activate the environment first, for example: conda activate ase"
        ) from exc

    return read


def ordered_water_indices(symbols: list[str]) -> tuple[int, list[int]]:
    oxygen_indices = [i for i, symbol in enumerate(symbols) if symbol == "O"]
    hydrogen_indices = [i for i, symbol in enumerate(symbols) if symbol == "H"]

    if len(oxygen_indices) != 1 or len(hydrogen_indices) != 2:
        raise ValueError(
            "Expected exactly one O atom and two H atoms per file, "
            f"got O={len(oxygen_indices)}, H={len(hydrogen_indices)}."
        )

    return oxygen_indices[0], hydrogen_indices


def dominant_axis(unit_vector: Any) -> str:
    np = _load_numpy()

    axis_labels = ["x", "y", "z"]
    axis_index = int(np.argmax(np.abs(unit_vector)))
    sign = "+" if unit_vector[axis_index] >= 0.0 else "-"
    return f"{sign}{axis_labels[axis_index]}"


def analyze_structure(path: Path) -> dict[str, object]:
    np = _load_numpy()
    read = _load_ase_read()

    atoms = read(path, format="vasp")
    symbols = atoms.get_chemical_symbols()
    oxygen_index, hydrogen_indices = ordered_water_indices(symbols)
    h1_index, h2_index = hydrogen_indices

    oh1 = atoms.get_distance(oxygen_index, h1_index, mic=True)
    oh2 = atoms.get_distance(oxygen_index, h2_index, mic=True)
    hoh_angle = atoms.get_angle(h1_index, oxygen_index, h2_index, mic=True)

    # Use the midpoint of the two H atoms to define the geometric dipole direction.
    oh_vectors = atoms.get_distances(
        oxygen_index, hydrogen_indices, mic=True, vector=True
    )
    dipole_vector = oh_vectors.mean(axis=0)
    dipole_norm = float(np.linalg.norm(dipole_vector))
    if dipole_norm == 0.0:
        raise ValueError("Dipole vector is zero; could not determine its direction.")

    dipole_unit = dipole_vector / dipole_norm

    return {
        "file": path.name,
        "oh1": oh1,
        "oh2": oh2,
        "hoh_angle": hoh_angle,
        "dipole_vector": dipole_vector,
        "dipole_unit": dipole_unit,
        "dominant_axis": dominant_axis(dipole_unit),
    }


def analyze_structures(paths: list[Path]) -> list[dict[str, object]]:
    return [analyze_structure(path) for path in paths]


def format_report(results: list[dict[str, object]]) -> str:
    headers = [
        "File",
        "OH1 (A)",
        "OH2 (A)",
        "HOH (deg)",
        "Dipole unit vector [x, y, z]",
        "Main dir",
    ]

    rows = []
    for result in results:
        dipole_unit = result["dipole_unit"]
        rows.append(
            [
                str(result["file"]),
                f"{result['oh1']:.6f}",
                f"{result['oh2']:.6f}",
                f"{result['hoh_angle']:.6f}",
                f"[{dipole_unit[0]: .6f}, {dipole_unit[1]: .6f}, {dipole_unit[2]: .6f}]",
                str(result["dominant_axis"]),
            ]
        )

    widths = []
    for column_index, header in enumerate(headers):
        column_width = max(len(header), *(len(row[column_index]) for row in rows))
        widths.append(column_width)

    def format_row(values: list[str]) -> str:
        return " | ".join(
            value.ljust(widths[index]) for index, value in enumerate(values)
        )

    separator = "-+-".join("-" * width for width in widths)

    lines = [
        "Water molecule geometry report",
        "Dipole direction is defined geometrically as the vector from O to the midpoint of the two H atoms.",
        "",
        format_row(headers),
        separator,
    ]
    lines.extend(format_row(row) for row in rows)
    return "\n".join(lines) + "\n"


def write_report(
    input_paths: list[Path], output_path: Path = Path("water_geometry_report.txt")
) -> Path:
    results = analyze_structures(input_paths)
    report = format_report(results)
    output_path.write_text(report, encoding="utf-8")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze water geometry in VASP/POSCAR-style files."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        default=None,
        help="Input files. If omitted, all '*.vasp' files in the current directory are used.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="water_geometry_report.txt",
        help="Output text file path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.inputs:
        input_paths = [Path(path) for path in args.inputs]
    else:
        input_paths = sorted(Path.cwd().glob("*.vasp"))

    if not input_paths:
        raise SystemExit("No input files found.")

    output_path = write_report(input_paths, Path(args.output))
    print(f"Wrote report to {output_path}")


if __name__ == "__main__":
    main()
