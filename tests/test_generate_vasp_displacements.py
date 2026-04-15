from __future__ import annotations

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from myax.modules.generate_vasp_displacements import (
    main,
    parse_delta,
    write_displaced_poscars,
)


POSCAR_SAMPLE = """H2
1.0
        5.0000000000         0.0000000000         0.0000000000
        0.0000000000         5.0000000000         0.0000000000
        0.0000000000         0.0000000000         5.0000000000
   H
    2
Cartesian
     0.000000000         0.000000000         0.000000000
     0.000000000         0.000000000         0.740000000
"""


class GenerateVaspDisplacementsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = Path(self.tmpdir.name)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def _read_vasp(self, path: Path):
        from ase.io import read

        return read(path, format="vasp")

    def test_write_displaced_poscars_creates_expected_directories(self) -> None:
        input_path = self.workdir / "POSCAR"
        output_root = self.workdir / "vib"
        input_path.write_text(POSCAR_SAMPLE, encoding="utf-8")

        written_paths = write_displaced_poscars(input_path, output_root, 0.01)

        expected_names = [
            "eq",
            "0x-",
            "0x+",
            "0y-",
            "0y+",
            "0z-",
            "0z+",
            "1x-",
            "1x+",
            "1y-",
            "1y+",
            "1z-",
            "1z+",
        ]

        self.assertEqual(len(written_paths), len(expected_names))
        self.assertEqual([path.parent.name for path in written_paths], expected_names)
        for name in expected_names:
            self.assertTrue((output_root / name / "POSCAR").exists())

    def test_write_displaced_poscars_writes_equilibrium_and_delta_shift(self) -> None:
        input_path = self.workdir / "POSCAR"
        output_root = self.workdir / "vib"
        input_path.write_text(POSCAR_SAMPLE, encoding="utf-8")

        write_displaced_poscars(input_path, output_root, 0.01)

        original_atoms = self._read_vasp(input_path)
        eq_atoms = self._read_vasp(output_root / "eq" / "POSCAR")
        displaced_atoms = self._read_vasp(output_root / "1z+" / "POSCAR")

        self.assertTrue(abs(original_atoms.positions - eq_atoms.positions).max() < 1.0e-12)
        self.assertTrue(abs(original_atoms.cell[:] - eq_atoms.cell[:]).max() < 1.0e-12)
        self.assertAlmostEqual(
            displaced_atoms.positions[1, 2] - original_atoms.positions[1, 2],
            0.01,
            places=12,
        )
        self.assertTrue(
            abs(displaced_atoms.positions[0] - original_atoms.positions[0]).max() < 1.0e-12
        )

    def test_parse_delta_rejects_non_positive_values(self) -> None:
        with self.assertRaises(SystemExit):
            parse_delta("0")

        with self.assertRaises(SystemExit):
            parse_delta("-0.1")

    def test_main_prompts_for_missing_values(self) -> None:
        input_path = self.workdir / "POSCAR"
        output_root = self.workdir / "vib"
        input_path.write_text(POSCAR_SAMPLE, encoding="utf-8")

        answers = iter([str(input_path), str(output_root), "0.02"])
        stdout = io.StringIO()

        with patch("sys.argv", ["generate_vasp_displacements"]):
            with patch("builtins.input", side_effect=lambda prompt: next(answers)):
                with patch("sys.stdout", stdout):
                    main()

        self.assertTrue((output_root / "eq" / "POSCAR").exists())
        self.assertIn("Wrote 13 POSCAR files", stdout.getvalue())

    def test_main_stops_when_input_file_is_missing(self) -> None:
        output_root = self.workdir / "vib"

        with patch(
            "sys.argv",
            ["generate_vasp_displacements", str(self.workdir / "missing.vasp"), str(output_root), "--delta", "0.01"],
        ):
            with self.assertRaises(SystemExit):
                main()


if __name__ == "__main__":
    unittest.main()
