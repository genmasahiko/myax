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
