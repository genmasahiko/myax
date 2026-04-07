#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


ORTHOGONAL_TOLERANCE = 1.0e-8
PERMUTATIONS = {
    "qe-to-vasp": [2, 0, 1],
    "vasp-to-qe": [1, 2, 0],
}


def _load_ase_io():
    try:
        from ase.io import read, write
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires ASE. Activate the environment first, for example: conda activate ase"
        ) from exc

    return read, write


def _load_numpy():
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires NumPy. Activate the environment first, for example: conda activate ase"
        ) from exc

    return np


def ensure_vasp_suffix(path: str) -> Path:
    output_path = Path(path)
    if output_path.suffix != ".vasp":
        output_path = output_path.with_suffix(".vasp")
    return output_path


def is_orthogonal(vector_a, vector_b, tolerance: float = ORTHOGONAL_TOLERANCE) -> bool:
    np = _load_numpy()
    norm_a = float(np.linalg.norm(vector_a))
    norm_b = float(np.linalg.norm(vector_b))
    if norm_a == 0.0 or norm_b == 0.0:
        return False
    return abs(float(np.dot(vector_a, vector_b))) <= tolerance * norm_a * norm_b


def classify_surface_normal(cell) -> str:
    vectors = cell[:]

    a1_is_normal = is_orthogonal(vectors[0], vectors[1]) and is_orthogonal(
        vectors[0], vectors[2]
    )
    a3_is_normal = is_orthogonal(vectors[2], vectors[0]) and is_orthogonal(
        vectors[2], vectors[1]
    )

    if a1_is_normal and a3_is_normal:
        return "ambiguous"
    if a1_is_normal:
        return "a1"
    if a3_is_normal:
        return "a3"
    return "unknown"


def format_from_normal_label(normal_label: str) -> str:
    if normal_label == "a1":
        return "vasp"
    if normal_label == "a3":
        return "qe"
    raise ValueError(f"Unsupported surface normal label: {normal_label}")


def direction_from_format(format_name: str) -> str:
    if format_name == "qe":
        return "qe-to-vasp"
    if format_name == "vasp":
        return "vasp-to-qe"
    raise ValueError(f"Unsupported format: {format_name}")


def target_format(format_name: str) -> str:
    if format_name == "qe":
        return "VASP-ESM"
    if format_name == "vasp":
        return "QE"
    raise ValueError(f"Unsupported format: {format_name}")


def display_format_name(format_name: str) -> str:
    if format_name == "qe":
        return "QE"
    if format_name == "vasp":
        return "VASP-ESM"
    raise ValueError(f"Unsupported format: {format_name}")


def resolve_format(atoms, input_func=None) -> str:
    if input_func is None:
        input_func = input

    normal_label = classify_surface_normal(atoms.cell)
    if normal_label == "a1":
        return "vasp"
    if normal_label == "a3":
        return "qe"
    if normal_label == "unknown":
        raise SystemExit(
            "Could not detect the surface-normal lattice vector automatically."
        )

    while True:
        answer = input_func(
            "Surface-normal lattice vector is ambiguous. Enter a1 for VASP-ESM or a3 for QE: "
        ).strip().lower()
        if answer in {"a1", "a3"}:
            return format_from_normal_label(answer)
        print("Please enter a1 for VASP-ESM or a3 for QE.")


def should_convert(format_name: str, input_func=None) -> bool:
    if input_func is None:
        input_func = input

    target_name = target_format(format_name)
    print(f"This POSCAR is {display_format_name(format_name)} format.")
    answer = input_func(
        f"Convert to {target_name} format? [ok/yes/y]: "
    ).strip().lower()
    return answer in {"ok", "yes", "y"}


def transform_atoms(atoms, direction: str):
    np = _load_numpy()

    if direction not in PERMUTATIONS:
        raise ValueError(f"Unsupported direction: {direction}")

    permutation = PERMUTATIONS[direction]
    transformed = atoms.copy()

    # Apply the same axis relabeling to atomic Cartesian coordinates.
    positions = transformed.get_positions()
    transformed.set_positions(positions[:, permutation])

    # Transform the cell as C' = P^T C P so the slab normal changes axis
    # while lengths, dot products, and handedness are preserved.
    cell = np.asarray(transformed.cell)
    transformed.set_cell(cell[np.ix_(permutation, permutation)], scale_atoms=False)

    return transformed


def convert_atoms(atoms, format_name: str):
    return transform_atoms(atoms, direction_from_format(format_name))


def convert_poscar(input_path: Path, output_path: Path, format_name: str) -> Path:
    read, write = _load_ase_io()
    atoms = read(input_path, format="vasp")
    converted = convert_atoms(atoms, format_name)
    write(output_path, converted, format="vasp", direct=False, sort=False)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Detect QE/VASP-ESM POSCAR axes automatically and convert after confirmation."
    )
    parser.add_argument("input_path", nargs="?", help="Input POSCAR file.")
    parser.add_argument("output_path", nargs="?", help="Output POSCAR file.")
    return parser.parse_args()


def main() -> None:
    read, _ = _load_ase_io()
    _load_numpy()
    args = parse_args()

    input_name = args.input_path or input("Input POSCAR file: ").strip()
    output_name = args.output_path or input("Output POSCAR file name: ").strip()

    if not input_name or not output_name:
        raise SystemExit("Input POSCAR and output POSCAR are required.")

    input_path = Path(input_name)
    output_path = ensure_vasp_suffix(output_name)

    if not input_path.exists():
        raise SystemExit(f"Input POSCAR file does not exist: {input_path}")

    atoms = read(input_path, format="vasp")
    format_name = resolve_format(atoms)
    if not should_convert(format_name):
        print("Conversion cancelled.")
        return

    written = convert_poscar(input_path, output_path, format_name)
    print(f"Wrote converted POSCAR file to {written}")


if __name__ == "__main__":
    main()
