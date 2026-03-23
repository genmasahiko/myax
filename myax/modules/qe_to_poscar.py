#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

ASSUME_ISOLATED_PATTERN = re.compile(
    r"assume_isolated\s*=\s*['\"]?([A-Za-z0-9_+-]+)['\"]?",
    re.IGNORECASE,
)


def _load_ase_io():
    try:
        from ase.io import read, write
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires ASE. Activate the environment first, for example: conda activate ase"
        ) from exc

    return read, write


def ensure_vasp_suffix(path: str) -> Path:
    output_path = Path(path)
    if output_path.suffix != ".vasp":
        output_path = output_path.with_suffix(".vasp")
    return output_path


def uses_esm(qe_input: Path) -> bool:
    text = qe_input.read_text(errors="replace").lower()
    match = ASSUME_ISOLATED_PATTERN.search(text)
    return match is not None and match.group(1) == "esm"


def shift_half_cell_along_z(atoms) -> None:
    scaled_positions = atoms.get_scaled_positions(wrap=False)
    scaled_positions[:, 2] = (scaled_positions[:, 2] + 0.5) % 1.0
    atoms.set_scaled_positions(scaled_positions)


def read_final_structure(qe_input: Path, qe_output: Path):
    read, _ = _load_ase_io()

    input_atoms = read(qe_input, format="espresso-in")
    output_atoms = read(qe_output, format="espresso-out", index=-1)

    # Keep the cell from the input in case QE output formatting varies.
    output_atoms.set_cell(input_atoms.cell)
    output_atoms.set_pbc(input_atoms.pbc)
    return output_atoms


def write_poscar(qe_input: Path, qe_output: Path, output_path: Path) -> Path:
    _, write = _load_ase_io()
    atoms = read_final_structure(qe_input, qe_output)

    if uses_esm(qe_input):
        shift_half_cell_along_z(atoms)

    write(output_path, atoms, format="vasp", direct=False, sort=False)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Quantum ESPRESSO input/output files to a POSCAR-format .vasp file."
    )
    parser.add_argument("qe_input", nargs="?", help="Quantum ESPRESSO input file.")
    parser.add_argument("qe_output", nargs="?", help="Quantum ESPRESSO output file.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output VASP file name. The .vasp suffix is added if needed.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    qe_input = args.qe_input or input("QE input file: ").strip()
    qe_output = args.qe_output or input("QE output file: ").strip()
    output_name = args.output or input("Output POSCAR file name: ").strip()

    if not qe_input or not qe_output or not output_name:
        raise SystemExit("QE input, QE output, and output file name are required.")

    qe_input_path = Path(qe_input)
    qe_output_path = Path(qe_output)
    output_path = ensure_vasp_suffix(output_name)

    if not qe_input_path.exists():
        raise SystemExit(f"QE input file does not exist: {qe_input_path}")
    if not qe_output_path.exists():
        raise SystemExit(f"QE output file does not exist: {qe_output_path}")

    written = write_poscar(qe_input_path, qe_output_path, output_path)
    print(f"Wrote POSCAR file to {written}")


if __name__ == "__main__":
    main()
