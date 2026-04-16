# My Analysis programs

This is a test source code for QuantumESPRESSO pre-/post-processing.<br>
As this program is for personal use, **I do not accept any responsibility for any losses from using this program.** <br>

***

## Installation

Use `make` in current directory. <br>
The compiler is set in `make.inc`, so you should revies it before compiling. <br>

## Python launcher

Run the Python entrypoint with `bin/myax`. <br>
This launcher moves to the repository root and executes `python3 -m myax`. <br>

## Python environment

The Python tools assume the user has already activated an environment with the required packages. <br>
For your current setup, use `conda activate ase` before running the program. <br>

```bash
conda activate ase
bin/myax
```

If `ase` or `numpy` is missing, the program will stop with a short message telling the user to activate the environment first. <br>

## Available Python features

`bin/myax` includes a converter from Quantum ESPRESSO input/output files to POSCAR-format `.vasp` files. <br>
If the QE output file is omitted, the converter writes the structure directly from the QE input file. <br>
When `assume_isolated = 'esm'` is used, the converter shifts atomic positions by half a cell along the fractional z direction before writing the POSCAR file. <br>
`bin/myax` also includes a POSCAR-to-POSCAR axis converter between the QE convention (surface normal = `z`) and the VASP-ESM convention (surface normal = `x`). <br>
The converter first detects whether the input is QE or VASP-ESM, prints the detected format, and asks for confirmation before writing the converted file. <br>
If the surface-normal direction is ambiguous, it asks the user whether the normal is `a1` or `a3`. <br>
This converter preserves Cartesian geometry and handedness, so it also works when the two in-plane lattice vectors are not orthogonal. <br>
`bin/myax` also includes a displacement input generator for vibrational calculations. <br>
It uses ASE to generate `eq`, `0x-`, `0x+` style displaced structures. <br>
For VASP it writes one `POSCAR` file into each output subdirectory. <br>
For Quantum ESPRESSO it uses a user-provided input file as the template and writes one `in` file into each output subdirectory. <br>
The displacement magnitude can be passed by `--delta` or entered interactively. <br>
`bin/myax` can also check finite-difference forces after SCF calculations in those displaced directories. <br>
It reads `OUTCAR` and `POSCAR` for VASP or `out` and `in` for Quantum ESPRESSO, computes forces from total-energy differences, and compares them with the forces in `eq`. <br>
