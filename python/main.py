from modules import grandpotential

modules = {
        "1": ("Grandpotential", grandpotential)
        }

def main():
    print("Hello, This program has some features.")

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
