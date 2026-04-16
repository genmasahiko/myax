#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

ATOMIC_POSITIONS_RE = re.compile(
    r"^\s*ATOMIC_POSITIONS\s*(?:[({]\s*([A-Za-z_]+)\s*[)}]|([A-Za-z_]+))?",
    re.IGNORECASE,
)
QE_CARD_RE = re.compile(
    r"^\s*(?:ATOMIC_SPECIES|ATOMIC_POSITIONS|K_POINTS|CELL_PARAMETERS|CONSTRAINTS|OCCUPATIONS)\b",
    re.IGNORECASE,
)
NAMELIST_RE = re.compile(r"^\s*&[A-Za-z_]+")


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
        description="Generate displaced input structures for finite-difference vibrational calculations."
    )
    parser.add_argument(
        "--format",
        choices=("vasp", "qe"),
        help="Input/output format. Defaults to vasp in non-interactive use.",
    )
    parser.add_argument("input_path", nargs="?", help="Input POSCAR or Quantum ESPRESSO input file.")
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


def read_atoms(input_path: Path, file_format: str = "vasp"):
    read, _ = _load_ase_io()
    ase_format = "espresso-in" if file_format == "qe" else "vasp"
    return read(input_path, format=ase_format)


def parse_delta(delta_value: str) -> float:
    try:
        delta = float(delta_value)
    except ValueError as exc:
        raise SystemExit("Displacement delta must be a positive number.") from exc

    if delta <= 0.0:
        raise SystemExit("Displacement delta must be greater than zero.")

    return delta


def ensure_inputs(
    input_name: str, output_name: str, delta: float, file_format: str = "vasp"
) -> tuple[Path, Path]:
    if not input_name:
        raise SystemExit("Input file is required.")
    if not output_name:
        raise SystemExit("Output parent directory is required.")

    input_path = Path(input_name)
    if not input_path.exists():
        label = "QE input file" if file_format == "qe" else "Input POSCAR file"
        raise SystemExit(f"{label} does not exist: {input_path}")

    output_root = Path(output_name)
    if delta <= 0.0:
        raise SystemExit("Displacement delta must be greater than zero.")

    return input_path, output_root


def normalize_format(format_name: str) -> str:
    normalized = format_name.strip().lower()
    if normalized in ("vasp", "qe"):
        return normalized
    raise SystemExit("Format must be either vasp or qe.")


def write_displaced_poscars(input_path: Path, output_root: Path, delta: float) -> list[Path]:
    _, write = _load_ase_io()
    Vibrations = _load_vibrations()

    atoms = read_atoms(input_path, "vasp")
    constraints = list(atoms.constraints)
    atoms.set_constraint([])
    vibrations = Vibrations(atoms, name="vib", delta=delta)

    written_paths: list[Path] = []

    # Write the equilibrium and each displaced structure in its own directory.
    for displacement, displaced_atoms in vibrations.iterdisplace():
        target_dir = output_root / displacement.name
        target_dir.mkdir(parents=True, exist_ok=True)

        output_path = target_dir / "POSCAR"
        displaced_atoms.set_constraint(constraints)
        write(output_path, displaced_atoms, format="vasp", direct=False, sort=False)
        written_paths.append(output_path)

    return written_paths


def write_displaced_qe_inputs(input_path: Path, output_root: Path, delta: float) -> list[Path]:
    Vibrations = _load_vibrations()

    atoms = read_atoms(input_path, "qe")
    atoms.set_constraint([])
    template = read_qe_template(input_path)
    vibrations = Vibrations(atoms, name="vib", delta=delta)

    written_paths: list[Path] = []

    # Keep the user's QE input settings and replace only atomic coordinates.
    for displacement, displaced_atoms in vibrations.iterdisplace():
        target_dir = output_root / displacement.name
        target_dir.mkdir(parents=True, exist_ok=True)

        output_path = target_dir / "in"
        output_path.write_text(render_qe_input(template, displaced_atoms), encoding="utf-8")
        written_paths.append(output_path)

    return written_paths


def write_displaced_inputs(
    input_path: Path, output_root: Path, delta: float, file_format: str
) -> list[Path]:
    if file_format == "qe":
        return write_displaced_qe_inputs(input_path, output_root, delta)
    return write_displaced_poscars(input_path, output_root, delta)


def read_qe_template(input_path: Path) -> dict[str, object]:
    lines = input_path.read_text(encoding="utf-8", errors="replace").splitlines()
    start_index = None
    unit = "angstrom"

    for i, line in enumerate(lines):
        match = ATOMIC_POSITIONS_RE.match(line)
        if match:
            start_index = i
            unit = (match.group(1) or match.group(2) or "alat").lower()
            break

    if start_index is None:
        raise SystemExit(f"ATOMIC_POSITIONS block was not found in QE input: {input_path}")
    if unit not in ("angstrom", "crystal"):
        raise SystemExit(
            f"ATOMIC_POSITIONS unit '{unit}' is not supported yet. Use angstrom or crystal."
        )

    atom_start = start_index + 1
    atom_end = atom_start
    while atom_end < len(lines):
        line = lines[atom_end]
        if not line.strip():
            break
        if QE_CARD_RE.match(line) or NAMELIST_RE.match(line):
            break
        atom_end += 1

    atom_lines = lines[atom_start:atom_end]
    if not atom_lines:
        raise SystemExit(f"ATOMIC_POSITIONS block has no atoms: {input_path}")

    return {
        "lines": lines,
        "unit": unit,
        "atom_start": atom_start,
        "atom_end": atom_end,
        "tails": [parse_qe_atom_tail(line) for line in atom_lines],
    }


def parse_qe_atom_tail(line: str) -> str:
    parts = line.split()
    if len(parts) < 4:
        raise SystemExit(f"Invalid ATOMIC_POSITIONS line: {line}")
    if len(parts) == 4:
        return ""
    return " " + " ".join(parts[4:])


def render_qe_input(template: dict[str, object], atoms) -> str:
    lines = list(template["lines"])
    atom_start = int(template["atom_start"])
    atom_end = int(template["atom_end"])
    unit = str(template["unit"])
    tails = list(template["tails"])

    if len(atoms) != len(tails):
        raise SystemExit("The number of displaced atoms does not match the QE template.")

    if unit == "crystal":
        positions = atoms.get_scaled_positions(wrap=False)
    else:
        positions = atoms.positions

    atom_lines = []
    for symbol, position, tail in zip(atoms.get_chemical_symbols(), positions, tails):
        atom_lines.append(
            f"{symbol:<2s} {position[0]:16.10f} {position[1]:16.10f} {position[2]:16.10f}{tail}"
        )

    output_lines = lines[:atom_start] + atom_lines + lines[atom_end:]
    return "\n".join(output_lines) + "\n"


def main() -> None:
    args = parse_args()

    if args.format is None and args.input_path is None:
        format_name = input("Calculation format (vasp/qe) [vasp]: ").strip() or "vasp"
    else:
        format_name = args.format or "vasp"
    file_format = normalize_format(format_name)

    input_prompt = "QE input file: " if file_format == "qe" else "Input POSCAR file: "
    input_name = args.input_path or input(input_prompt).strip()
    output_name = args.output_root or input("Output parent directory: ").strip()

    if args.delta is None:
        delta_value = input("Displacement delta in Angstrom: ").strip()
        delta = parse_delta(delta_value)
    else:
        delta = parse_delta(str(args.delta))

    input_path, output_root = ensure_inputs(input_name, output_name, delta, file_format)
    written_paths = write_displaced_inputs(input_path, output_root, delta, file_format)
    output_label = "QE input files" if file_format == "qe" else "POSCAR files"
    print(f"Wrote {len(written_paths)} {output_label} under {output_root}")


if __name__ == "__main__":
    main()
