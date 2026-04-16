#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


RY_TO_EV = 13.605693122994
BOHR_TO_ANGSTROM = 0.529177210903
RY_PER_BOHR_TO_EV_PER_ANGSTROM = RY_TO_EV / BOHR_TO_ANGSTROM

DISPLACEMENT_RE = re.compile(r"^(\d+)([xyz])([+-])$")
DIRECTION_INDEX = {"x": 0, "y": 1, "z": 2}


@dataclass(frozen=True)
class ScfResult:
    energy: float
    forces: list[list[float]]


@dataclass(frozen=True)
class ForceComparison:
    atom_index: int
    direction: str
    eq_force: float
    fd_force: float

    @property
    def diff(self) -> float:
        return self.fd_force - self.eq_force

    @property
    def abs_diff(self) -> float:
        return abs(self.diff)


def _load_ase_io():
    try:
        from ase.io import read
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This feature requires ASE. Activate the environment first, for example: conda activate ase"
        ) from exc

    return read


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare equilibrium forces with finite-difference forces from SCF energies."
    )
    parser.add_argument("root", nargs="?", help="Parent directory made by finite-difference POSCAR generation.")
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.0e-3,
        help="Allowed absolute force difference in eV/Angstrom. Default: 1.0e-3",
    )
    parser.add_argument("--csv", dest="csv_path", help="Optional CSV output path.")
    return parser.parse_args()


def find_output_file(directory: Path) -> Path:
    for filename in ("OUTCAR", "out"):
        path = directory / filename
        if path.exists():
            return path

    raise SystemExit(f"SCF output file was not found in {directory} (expected OUTCAR or out).")


def read_poscar_positions(poscar_path: Path) -> list[list[float]]:
    if not poscar_path.exists():
        raise SystemExit(f"POSCAR file was not found: {poscar_path}")

    read = _load_ase_io()
    try:
        atoms = read(poscar_path, format="vasp")
    except Exception as exc:
        raise SystemExit(f"Failed to read POSCAR with ASE: {poscar_path}: {exc}") from exc

    return atoms.positions.tolist()


def read_scf_result(output_path: Path) -> ScfResult:
    try:
        return _read_scf_result_with_ase(output_path)
    except Exception as ase_exc:
        try:
            if output_path.name == "OUTCAR":
                return _parse_vasp_outcar(output_path)
            if output_path.name == "out":
                return _parse_qe_out(output_path)
        except Exception as fallback_exc:
            raise SystemExit(
                f"Failed to read SCF result from {output_path}: ASE error: {ase_exc}; "
                f"fallback parser error: {fallback_exc}"
            ) from fallback_exc

        raise SystemExit(f"Failed to read SCF result from {output_path}: {ase_exc}") from ase_exc


def _read_scf_result_with_ase(output_path: Path) -> ScfResult:
    read = _load_ase_io()
    if output_path.name == "OUTCAR":
        atoms = read(output_path, index=-1, format="vasp-out")
    elif output_path.name == "out":
        atoms = read(output_path, index=-1, format="espresso-out")
    else:
        atoms = read(output_path, index=-1)

    return ScfResult(
        energy=float(atoms.get_potential_energy()),
        forces=atoms.get_forces().tolist(),
    )


def _parse_vasp_outcar(output_path: Path) -> ScfResult:
    energy = None
    force_blocks: list[list[list[float]]] = []
    lines = output_path.read_text(encoding="utf-8", errors="replace").splitlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        energy_match = re.search(r"free\s+energy\s+TOTEN\s+=\s+([-+0-9.Ee]+)", line)
        if energy_match:
            energy = float(energy_match.group(1))

        if "TOTAL-FORCE" in line and "eV/Angst" in line:
            block: list[list[float]] = []
            i += 2
            while i < len(lines):
                parts = lines[i].split()
                if len(parts) < 6:
                    break
                try:
                    block.append([float(parts[3]), float(parts[4]), float(parts[5])])
                except ValueError:
                    break
                i += 1
            if block:
                force_blocks.append(block)
            continue
        i += 1

    if energy is None:
        raise ValueError("TOTEN was not found.")
    if not force_blocks:
        raise ValueError("TOTAL-FORCE block was not found.")

    return ScfResult(energy=energy, forces=force_blocks[-1])


def _parse_qe_out(output_path: Path) -> ScfResult:
    energy = None
    force_blocks: list[list[list[float]]] = []
    lines = output_path.read_text(encoding="utf-8", errors="replace").splitlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        energy_match = re.search(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", line)
        if energy_match:
            energy = float(energy_match.group(1)) * RY_TO_EV

        if "Forces acting on atoms" in line:
            block: list[list[float]] = []
            i += 1
            while i < len(lines):
                force_match = re.search(
                    r"force\s+=\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
                    lines[i],
                )
                if force_match:
                    block.append(
                        [
                            float(force_match.group(1)) * RY_PER_BOHR_TO_EV_PER_ANGSTROM,
                            float(force_match.group(2)) * RY_PER_BOHR_TO_EV_PER_ANGSTROM,
                            float(force_match.group(3)) * RY_PER_BOHR_TO_EV_PER_ANGSTROM,
                        ]
                    )
                elif block:
                    break
                i += 1
            if block:
                force_blocks.append(block)
            continue
        i += 1

    if energy is None:
        raise ValueError("final total energy line was not found.")
    if not force_blocks:
        raise ValueError("force block was not found.")

    return ScfResult(energy=energy, forces=force_blocks[-1])


def collect_displacement_pairs(root: Path) -> dict[tuple[int, str], dict[str, Path]]:
    pairs: dict[tuple[int, str], dict[str, Path]] = {}

    for path in sorted(root.iterdir()):
        if not path.is_dir():
            continue

        match = DISPLACEMENT_RE.match(path.name)
        if not match:
            continue

        atom_index = int(match.group(1))
        direction = match.group(2)
        sign = match.group(3)
        pairs.setdefault((atom_index, direction), {})[sign] = path

    return pairs


def compare_forces(root: Path) -> list[ForceComparison]:
    if not root.exists():
        raise SystemExit(f"Finite-difference parent directory does not exist: {root}")

    eq_dir = root / "eq"
    if not eq_dir.is_dir():
        raise SystemExit(f"Equilibrium directory was not found: {eq_dir}")

    eq_output = find_output_file(eq_dir)
    eq_result = read_scf_result(eq_output)
    eq_positions = read_poscar_positions(eq_dir / "POSCAR")

    comparisons: list[ForceComparison] = []
    pairs = collect_displacement_pairs(root)
    if not pairs:
        raise SystemExit(f"No displacement directories like 0x- or 0x+ were found under {root}")

    for (atom_index, direction), dirs in sorted(pairs.items()):
        if "-" not in dirs or "+" not in dirs:
            raise SystemExit(f"Both plus and minus directories are required for {atom_index}{direction}.")

        minus_dir = dirs["-"]
        plus_dir = dirs["+"]
        minus_result = read_scf_result(find_output_file(minus_dir))
        plus_result = read_scf_result(find_output_file(plus_dir))

        axis = DIRECTION_INDEX[direction]
        minus_positions = read_poscar_positions(minus_dir / "POSCAR")
        plus_positions = read_poscar_positions(plus_dir / "POSCAR")
        if atom_index >= len(eq_positions):
            raise SystemExit(f"Atom index {atom_index} is outside eq/POSCAR atom list.")
        if atom_index >= len(minus_positions) or atom_index >= len(plus_positions):
            raise SystemExit(f"Atom index {atom_index} is outside displaced POSCAR atom list.")

        displacement_width = plus_positions[atom_index][axis] - minus_positions[atom_index][axis]

        if abs(displacement_width) < 1.0e-15:
            raise SystemExit(f"Displacement width is zero for {atom_index}{direction}.")

        # Force is the negative derivative of total energy with respect to displacement.
        fd_force = -(plus_result.energy - minus_result.energy) / displacement_width
        comparisons.append(
            ForceComparison(
                atom_index=atom_index,
                direction=direction,
                eq_force=_force_component(eq_result, atom_index, axis, eq_output),
                fd_force=fd_force,
            )
        )

    return comparisons


def print_report(comparisons: list[ForceComparison], threshold: float) -> None:
    print("Finite-difference force check matrices (atom rows x Cartesian columns)")
    print()
    _print_matrix("F_eq (eV/Angstrom)", comparisons, "eq_force")
    print()
    _print_matrix("F_fd from energy differences (eV/Angstrom)", comparisons, "fd_force")
    print()
    _print_matrix("F_fd - F_eq (eV/Angstrom)", comparisons, "diff")
    print()

    worst = max(comparisons, key=lambda item: item.abs_diff)
    status = "PASS" if worst.abs_diff <= threshold else "FAIL"
    print(
        f"Max abs diff: {worst.abs_diff:.8e} eV/Angstrom "
        f"at atom {worst.atom_index} ({worst.atom_index + 1}), dir {worst.direction}"
    )
    print(f"Threshold: {threshold:.8e} eV/Angstrom")
    print(f"Result: {status}")


def _print_matrix(title: str, comparisons: list[ForceComparison], value_name: str) -> None:
    matrix = _comparison_matrix(comparisons, value_name)

    print(title)
    print(f"{'atom(0)':>8s} {'atom(1)':>8s} {'x':>16s} {'y':>16s} {'z':>16s}")
    print("-" * 70)
    for atom_index in sorted(matrix):
        row = matrix[atom_index]
        print(
            f"{atom_index:8d} {atom_index + 1:8d} "
            f"{_format_matrix_value(row.get('x')):>16s} "
            f"{_format_matrix_value(row.get('y')):>16s} "
            f"{_format_matrix_value(row.get('z')):>16s}"
        )


def _comparison_matrix(comparisons: list[ForceComparison], value_name: str) -> dict[int, dict[str, float]]:
    matrix: dict[int, dict[str, float]] = {}
    for item in comparisons:
        matrix.setdefault(item.atom_index, {})[item.direction] = float(getattr(item, value_name))
    return matrix


def _format_matrix_value(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.8f}"


def _force_component(result: ScfResult, atom_index: int, axis: int, output_path: Path) -> float:
    try:
        return result.forces[atom_index][axis]
    except IndexError as exc:
        raise SystemExit(f"Force component for atom {atom_index} was not found in {output_path}.") from exc


def write_csv(csv_path: Path, comparisons: list[ForceComparison]) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "atom_index_0",
                "atom_index_1",
                "direction",
                "f_eq_ev_a",
                "f_fd_ev_a",
                "diff_ev_a",
                "abs_diff_ev_a",
            ]
        )
        for item in comparisons:
            writer.writerow(
                [
                    item.atom_index,
                    item.atom_index + 1,
                    item.direction,
                    item.eq_force,
                    item.fd_force,
                    item.diff,
                    item.abs_diff,
                ]
            )


def main() -> None:
    args = parse_args()
    root_name = args.root or input("Finite-difference parent directory: ").strip()
    if not root_name:
        raise SystemExit("Finite-difference parent directory is required.")
    if args.threshold < 0.0:
        raise SystemExit("Threshold must be non-negative.")

    comparisons = compare_forces(Path(root_name))
    print_report(comparisons, args.threshold)

    if args.csv_path:
        csv_path = Path(args.csv_path)
        write_csv(csv_path, comparisons)
        print(f"Wrote CSV report to {csv_path}")


if __name__ == "__main__":
    main()
