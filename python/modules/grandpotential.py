from ase.io import read, write, espresso
import numpy as np
import os
from pathlib import Path

def generate_input(qe_input, qe_output):

    # Read the QE input file using ASE
    # params include the 'Namelist', cards include the 'Card' section
    with open(qe_input, 'r') as fin:
        params, cards = espresso.read_fortran_namelist(fin)

    old_params = params.copy()

    # Check if the calculation type is 'relax' or not
    # In the case of 'relax', extract the final coordinates
    if params['control'].get('calculation') == 'relax':
        params['control']['calculation'] = 'scf'
        
        with open(qe_output, 'r') as fout:
            lines = fout.readlines()

            for i, line in enumerate(lines):
                if "Begin final coordinates" in line:
                    new_coords = lines[i+2:i+2+params['system']['nat']+1]

    charge = 0

    # Check if Constant-mu calculation or not
    # If so, extract the final tot_charge
    if params['control'].get('lfcp'):
        params['control']['lfcp'] = False

        with open(qe_output, 'r') as fout:
            lines = fout.readlines()

            # tot_charge appears many times,
            # so we need to find the last one
            for line in reversed(lines):
                if "FCP: Total Charge =" in line:
                    charge = line.split("=")[1].strip()
                    charge = float(charge)
                    break
    print("Calculation has ended for charge =", charge)

    # Extract the number of points and step width you need
    npoints, stepw = input("\nThe number of points / step width: ").split()
    print()

    # npoints and stepw should be integers and floats, respectively
    npoints = int(npoints)
    stepw = float(stepw)

    # Make the new charges to be calculated
    charges = [ charge + (i * stepw)
               for i in range(- npoints // 2, npoints // 2 + 1) if i != 0]
    print("New charges to be calculated:", charges)

    # Check if the grandpotential directory already exists
    # If exists, ask the user if they want to overwrite it
    if os.path.exists("grandpotential"):
        print("Directory 'grandpotential' already exists.\n")
        response = input("Are you sure you want to overwrite it? (y/n): ").strip().lower()
        if response != 'y':
            print("Exiting")
            return
    os.makedirs("grandpotential", exist_ok=True)

    # To avoid overwriting the existing results,
    # Check the existing directories in grandpotential
    dirs = [int(p.name) for p in Path("grandpotential").iterdir() if p.is_dir()]
    num_dirs = max(dirs) if dirs else 0

    # Write the new input file
    # To avoid overwriting the existing results,
    # The name of the directory is (the largest number of the directory) + 1
    j = 0
    for i, charge in enumerate(charges, start=num_dirs+1):
        os.makedirs(f"grandpotential/{i}", exist_ok=True)
        filename = f"grandpotential/{i}/in"

        # Prepare the new charge 
        params['system']['tot_charge'] = charges[j]

        # Write the new input file (Namelist only)
        with open(filename, 'w') as f:
            espresso.write_fortran_namelist(f, input_data = params, binary = 'pw',)

        # This function append the 'EOF' at the end of the file (crazy)
        # In this program, we have to remove it
        with open(filename, 'r') as f:
            lines = f.readlines()
        with open(filename, 'w') as f:
            f.writelines(lines[:-1])

        # Add the new coordinates to the input file
        with open(filename, 'a') as f:
            # The new coordinates is used if the calculation was 'relax'
            for line in cards:
                if "ATOMIC_POSITIONS" in line and old_params['control']['calculation'] == 'relax':
                    for line in new_coords:
                        f.write(line)
                    break
                f.write(line + '\n')

        j += 1

    print("\nCreated input files")

def main():
    print()
    print("Making directory for calculating grandpotential\n")
    
    # Get the QE input and output files from the user
    try:
        qe_input, qe_output = input("QE input / output: ").split()
        print()
    except ValueError:
        print("Please provide both QE input and output files.")
        return

    # Check if the files exist
    if not os.path.exists(qe_input):
        print(f"QE input file {qe_input} does not exist.")
        return
    if not os.path.exists(qe_output):
        print(f"QE output file {qe_output} does not exist.")
        return

    # Make the input files for the grand potential calculation
    generate_input(qe_input, qe_output)

if __name__ == "__main__":
    main()
