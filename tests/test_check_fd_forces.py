from __future__ import annotations

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from myax.modules.check_fd_forces import (
    RY_PER_BOHR_TO_EV_PER_ANGSTROM,
    RY_TO_EV,
    compare_forces,
    find_output_file,
    main,
    read_scf_result,
)


POSCAR_EQ = """H
1.0
        5.0000000000         0.0000000000         0.0000000000
        0.0000000000         5.0000000000         0.0000000000
        0.0000000000         0.0000000000         5.0000000000
   H
    1
Cartesian
     0.000000000         0.000000000         0.000000000
"""

POSCAR_MINUS = """H
1.0
        5.0000000000         0.0000000000         0.0000000000
        0.0000000000         5.0000000000         0.0000000000
        0.0000000000         0.0000000000         5.0000000000
   H
    1
Cartesian
    -0.010000000         0.000000000         0.000000000
"""

POSCAR_PLUS = """H
1.0
        5.0000000000         0.0000000000         0.0000000000
        0.0000000000         5.0000000000         0.0000000000
        0.0000000000         0.0000000000         5.0000000000
   H
    1
Cartesian
     0.010000000         0.000000000         0.000000000
"""


def vasp_outcar(energy: float, fx: float = 0.0, fy: float = 0.0, fz: float = 0.0) -> str:
    return f""" free  energy   TOTEN  =      {energy: .8f} eV
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
   0.00000000  0.00000000  0.00000000   {fx: .8f}  {fy: .8f}  {fz: .8f}
 -----------------------------------------------------------------------------------
"""


class CheckFdForcesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = Path(self.tmpdir.name)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def _make_vasp_tree(self) -> Path:
        root = self.workdir / "vib"
        for dirname, poscar, energy, force in [
            ("eq", POSCAR_EQ, -1.0, 0.2),
            ("0x-", POSCAR_MINUS, -1.0, 0.0),
            ("0x+", POSCAR_PLUS, -1.004, 0.0),
        ]:
            directory = root / dirname
            directory.mkdir(parents=True)
            (directory / "POSCAR").write_text(poscar, encoding="utf-8")
            (directory / "OUTCAR").write_text(vasp_outcar(energy, fx=force), encoding="utf-8")
        return root

    def test_compare_forces_uses_poscar_displacement_width(self) -> None:
        root = self._make_vasp_tree()

        comparisons = compare_forces(root)

        self.assertEqual(len(comparisons), 1)
        self.assertEqual(comparisons[0].atom_index, 0)
        self.assertEqual(comparisons[0].direction, "x")
        self.assertAlmostEqual(comparisons[0].eq_force, 0.2)
        self.assertAlmostEqual(comparisons[0].fd_force, 0.2)
        self.assertAlmostEqual(comparisons[0].abs_diff, 0.0)

    def test_find_output_file_prefers_vasp_then_qe_name(self) -> None:
        directory = self.workdir / "calc"
        directory.mkdir()
        (directory / "out").write_text("", encoding="utf-8")
        self.assertEqual(find_output_file(directory), directory / "out")

        (directory / "OUTCAR").write_text("", encoding="utf-8")
        self.assertEqual(find_output_file(directory), directory / "OUTCAR")

    def test_read_qe_out_fallback_converts_units(self) -> None:
        output = self.workdir / "out"
        output.write_text(
            """!    total energy              =      -2.00000000 Ry
Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.01000000    0.00000000   -0.02000000
""",
            encoding="utf-8",
        )

        result = read_scf_result(output)

        self.assertAlmostEqual(result.energy, -2.0 * RY_TO_EV)
        self.assertAlmostEqual(result.forces[0][0], 0.01 * RY_PER_BOHR_TO_EV_PER_ANGSTROM)
        self.assertAlmostEqual(result.forces[0][2], -0.02 * RY_PER_BOHR_TO_EV_PER_ANGSTROM)

    def test_main_prompts_and_writes_csv(self) -> None:
        root = self._make_vasp_tree()
        csv_path = self.workdir / "report.csv"
        answers = iter([str(root)])
        stdout = io.StringIO()

        with patch("sys.argv", ["check_fd_forces", "--csv", str(csv_path)]):
            with patch("builtins.input", side_effect=lambda prompt: next(answers)):
                with patch("sys.stdout", stdout):
                    main()

        self.assertTrue(csv_path.exists())
        self.assertIn("Result: PASS", stdout.getvalue())
        self.assertIn("atom_index_0,atom_index_1,direction", csv_path.read_text(encoding="utf-8"))

    def test_missing_pair_stops_with_clear_error(self) -> None:
        root = self.workdir / "vib"
        (root / "eq").mkdir(parents=True)
        (root / "0x+").mkdir()
        (root / "eq" / "POSCAR").write_text(POSCAR_EQ, encoding="utf-8")
        (root / "eq" / "OUTCAR").write_text(vasp_outcar(-1.0, fx=0.2), encoding="utf-8")

        with self.assertRaises(SystemExit):
            compare_forces(root)


if __name__ == "__main__":
    unittest.main()
