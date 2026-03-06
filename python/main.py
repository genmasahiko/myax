from modules import grandpotential
from modules import profile_vasp

modules = {
        "1": ("Grandpotential", grandpotential),
        "2": ("Analyze profiling in VASP", profile_vasp)
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
            _, module = modules[choice]
            module.main()
            break
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()
