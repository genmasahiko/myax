from __future__ import annotations

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from myax.modules.convert_poscar_axes import (
    classify_surface_normal,
    convert_poscar,
    main,
    resolve_format,
)


QE_SAMPLE = """Al001_bc1
1.0
        5.7269001007         0.0000000000         0.0000000000
       -2.8634500503         4.9591840273         0.0000000000
        0.0000000000         0.0000000000        11.9999971390
   Al
    4
Cartesian
     0.000000000         0.000000000         5.999998569
     2.863450050         0.000000000         5.999998569
    -1.431725025         2.479592014         5.999998569
     1.431725025         2.479592014         5.999998569
"""


VASP_SAMPLE = """Al001_bc1
1.0
       11.9999971390         0.0000000000         0.0000000000
        0.0000000000         5.7269001007         0.0000000000
        0.0000000000        -2.8634500503         4.9591840273
   Al
    4
Cartesian
     5.999998569         0.000000000         0.000000000
     5.999998569         2.863450050         0.000000000
     5.999998569        -1.431725025         2.479592014
     5.999998569         1.431725025         2.479592014
"""


DIRECT_SAMPLE = """DirectCell
1.0
        2.0000000000         0.0000000000         0.0000000000
        0.5000000000         3.0000000000         0.0000000000
        0.0000000000         0.0000000000         8.0000000000
   H
    1
Direct
     0.250000000         0.500000000         0.750000000
"""


AMBIGUOUS_SAMPLE = """Cubic
1.0
        4.0000000000         0.0000000000         0.0000000000
        0.0000000000         4.0000000000         0.0000000000
        0.0000000000         0.0000000000         8.0000000000
   H
    1
Cartesian
     1.000000000         2.000000000         3.000000000
"""


UNKNOWN_SAMPLE = """Skew
1.0
        4.0000000000         0.0000000000         1.0000000000
        1.5000000000         4.0000000000         0.0000000000
        0.0000000000         2.0000000000         8.0000000000
   H
    1
Cartesian
     1.000000000         2.000000000         3.000000000
"""


class ConvertPoscarAxesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = Path(self.tmpdir.name)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def _read_vasp(self, path: Path):
        from ase.io import read

        return read(path, format="vasp")

    def _pairwise_dot_products(self, atoms) -> list[float]:
        cell = atoms.cell[:]
        return sorted(
            [
                float(cell[0] @ cell[0]),
                float(cell[1] @ cell[1]),
                float(cell[2] @ cell[2]),
                float(cell[0] @ cell[1]),
                float(cell[1] @ cell[2]),
                float(cell[2] @ cell[0]),
            ]
        )

    def _sorted_lengths(self, atoms) -> list[float]:
        return sorted(float(length) for length in atoms.cell.lengths())

    def test_classify_surface_normal_for_qe_and_vasp(self) -> None:
        qe_path = self.workdir / "qe.vasp"
        vasp_path = self.workdir / "vasp.vasp"
        qe_path.write_text(QE_SAMPLE, encoding="utf-8")
        vasp_path.write_text(VASP_SAMPLE, encoding="utf-8")

        qe_atoms = self._read_vasp(qe_path)
        vasp_atoms = self._read_vasp(vasp_path)

        self.assertEqual(classify_surface_normal(qe_atoms.cell), "a3")
        self.assertEqual(classify_surface_normal(vasp_atoms.cell), "a1")

    def test_ambiguous_detection_asks_user_for_surface_normal(self) -> None:
        ambiguous_path = self.workdir / "ambiguous.vasp"
        ambiguous_path.write_text(AMBIGUOUS_SAMPLE, encoding="utf-8")
        atoms = self._read_vasp(ambiguous_path)

        replies = iter(["a2", "a1"])
        with patch("builtins.print") as mock_print:
            format_name = resolve_format(atoms, input_func=lambda _: next(replies))

        self.assertEqual(format_name, "vasp")
        mock_print.assert_called_with("Please enter a1 for VASP-ESM or a3 for QE.")

    def test_unknown_detection_stops_without_conversion(self) -> None:
        unknown_path = self.workdir / "unknown.vasp"
        unknown_path.write_text(UNKNOWN_SAMPLE, encoding="utf-8")
        atoms = self._read_vasp(unknown_path)

        with self.assertRaises(SystemExit):
            resolve_format(atoms, input_func=lambda _: "a1")

    def test_qe_to_vasp_and_back_preserves_geometry(self) -> None:
        qe_path = self.workdir / "qe.vasp"
        vasp_path = self.workdir / "vasp.vasp"
        roundtrip_path = self.workdir / "roundtrip.vasp"
        qe_path.write_text(QE_SAMPLE, encoding="utf-8")

        convert_poscar(qe_path, vasp_path, "qe")
        convert_poscar(vasp_path, roundtrip_path, "vasp")

        qe_atoms = self._read_vasp(qe_path)
        vasp_atoms = self._read_vasp(vasp_path)
        roundtrip_atoms = self._read_vasp(roundtrip_path)

        self.assertAlmostEqual(vasp_atoms.cell[0, 0], qe_atoms.cell[2, 2], places=9)
        self.assertAlmostEqual(vasp_atoms.cell[1, 1], qe_atoms.cell[0, 0], places=9)
        self.assertAlmostEqual(vasp_atoms.cell[2, 1], qe_atoms.cell[1, 0], places=9)
        self.assertAlmostEqual(vasp_atoms.cell[2, 2], qe_atoms.cell[1, 1], places=9)
        self.assertAlmostEqual(
            vasp_atoms.positions[0, 0], qe_atoms.positions[0, 2], places=9
        )
        self.assertAlmostEqual(
            vasp_atoms.positions[2, 2], qe_atoms.positions[2, 1], places=9
        )

        self.assertTrue(
            max(
                abs(left - right)
                for left, right in zip(
                    self._sorted_lengths(qe_atoms), self._sorted_lengths(vasp_atoms)
                )
            )
            < 1.0e-10
        )
        self.assertEqual(
            self._pairwise_dot_products(qe_atoms),
            self._pairwise_dot_products(vasp_atoms),
        )
        self.assertGreater(qe_atoms.cell.handedness, 0)
        self.assertEqual(qe_atoms.cell.handedness, vasp_atoms.cell.handedness)

        self.assertTrue(
            abs(qe_atoms.get_positions() - roundtrip_atoms.get_positions()).max()
            < 1.0e-10
        )
        self.assertTrue(abs(qe_atoms.cell[:] - roundtrip_atoms.cell[:]).max() < 1.0e-10)

    def test_direct_input_is_written_as_cartesian(self) -> None:
        input_path = self.workdir / "direct.vasp"
        output_path = self.workdir / "converted.vasp"
        input_path.write_text(DIRECT_SAMPLE, encoding="utf-8")

        convert_poscar(input_path, output_path, "qe")

        output_text = output_path.read_text(encoding="utf-8")
        self.assertIn("Cartesian", output_text)
        self.assertNotIn("Direct\n", output_text)

    def test_main_confirms_before_writing_output(self) -> None:
        input_path = self.workdir / "qe.vasp"
        output_path = self.workdir / "converted.vasp"
        input_path.write_text(QE_SAMPLE, encoding="utf-8")

        answers = iter([str(input_path), str(output_path), "yes"])
        stdout = io.StringIO()

        with patch("sys.argv", ["convert_poscar_axes"]):
            with patch("builtins.input", side_effect=lambda prompt: next(answers)):
                with patch("sys.stdout", stdout):
                    main()

        self.assertTrue(output_path.exists())
        self.assertIn("Wrote converted POSCAR file", stdout.getvalue())

    def test_main_cancels_without_writing_output(self) -> None:
        input_path = self.workdir / "vasp.vasp"
        output_path = self.workdir / "cancelled.vasp"
        input_path.write_text(VASP_SAMPLE, encoding="utf-8")

        answers = iter([str(input_path), str(output_path), "no"])
        stdout = io.StringIO()

        with patch("sys.argv", ["convert_poscar_axes"]):
            with patch("builtins.input", side_effect=lambda prompt: next(answers)):
                with patch("sys.stdout", stdout):
                    main()

        self.assertFalse(output_path.exists())
        self.assertIn("Conversion cancelled.", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
