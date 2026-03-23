#!/usr/bin/env python3

import importlib
import sys

modules = {
        "1": ("Grandpotential", "myax.modules.grandpotential"),
        "2": ("Analyze profiling in VASP", "myax.modules.profile_vasp"),
        "3": ("Analyze water POSCAR files", "myax.modules.analyze_water_poscars"),
        "4": ("Convert QE input/output to POSCAR", "myax.modules.qe_to_poscar")
        }

def main():
    print("Hello, This program has some features.\n")

    while True:
        for key, (description, _) in modules.items():
            print(f"    {key}. {description}")
        print("    q. Quit\n")

        choice = input("Enter your choice: ")
        print("--"*30)

        if choice == "q":
            print("Quitting the program.")
            break
        elif choice in modules.keys():
            _, module_name = modules[choice]
            try:
                module = importlib.import_module(module_name)
                module.main()
            except ModuleNotFoundError as exc:
                print(f"Missing Python dependency: {exc.name}")
                print("Activate the required environment first, for example: conda activate ase")
                sys.exit(1)
            break
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()
