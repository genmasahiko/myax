#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def _load_ase_io():
    try:
        from ase.io import read, write
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires ASE. Activate the environment first, for example: conda activate ase"
        ) from exc

    return read, write


def _load_vibrations():
    try:
        from ase.vibrations import Vibrations
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires ASE. Activate the environment first, for example: conda activate ase"
        ) from exc

    return Vibrations


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate displaced POSCAR structures for finite-difference vibrational calculations."
    )
    parser.add_argument("input_path", nargs="?", help="Input POSCAR file.")
    parser.add_argument(
        "output_root",
        nargs="?",
        help="Output parent directory for displaced structures.",
    )
    parser.add_argument(
        "--delta",
        type=float,
        help="Displacement magnitude in Angstrom. If omitted, prompt interactively.",
    )
    return parser.parse_args()


def read_atoms(input_path: Path):
    read, _ = _load_ase_io()
    return read(input_path, format="vasp")


def parse_delta(delta_value: str) -> float:
    try:
        delta = float(delta_value)
    except ValueError as exc:
        raise SystemExit("Displacement delta must be a positive number.") from exc

    if delta <= 0.0:
        raise SystemExit("Displacement delta must be greater than zero.")

    return delta


def ensure_inputs(input_name: str, output_name: str, delta: float) -> tuple[Path, Path]:
    if not input_name:
        raise SystemExit("Input POSCAR file is required.")
    if not output_name:
        raise SystemExit("Output parent directory is required.")

    input_path = Path(input_name)
    if not input_path.exists():
        raise SystemExit(f"Input POSCAR file does not exist: {input_path}")

    output_root = Path(output_name)
    if delta <= 0.0:
        raise SystemExit("Displacement delta must be greater than zero.")

    return input_path, output_root


def write_displaced_poscars(input_path: Path, output_root: Path, delta: float) -> list[Path]:
    _, write = _load_ase_io()
    Vibrations = _load_vibrations()

    atoms = read_atoms(input_path)
    vibrations = Vibrations(atoms, name="vib", delta=delta)

    written_paths: list[Path] = []

    # Write the equilibrium and each displaced structure in its own directory.
    for displacement, displaced_atoms in vibrations.iterdisplace():
        target_dir = output_root / displacement.name
        target_dir.mkdir(parents=True, exist_ok=True)

        output_path = target_dir / "POSCAR"
        write(output_path, displaced_atoms, format="vasp", direct=False, sort=False)
        written_paths.append(output_path)

    return written_paths


def main() -> None:
    args = parse_args()

    input_name = args.input_path or input("Input POSCAR file: ").strip()
    output_name = args.output_root or input("Output parent directory: ").strip()

    if args.delta is None:
        delta_value = input("Displacement delta in Angstrom: ").strip()
        delta = parse_delta(delta_value)
    else:
        delta = parse_delta(str(args.delta))

    input_path, output_root = ensure_inputs(input_name, output_name, delta)
    written_paths = write_displaced_poscars(input_path, output_root, delta)
    print(f"Wrote {len(written_paths)} POSCAR files under {output_root}")


if __name__ == "__main__":
    main()
