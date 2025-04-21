// This program check which of the programs to be run
//


#include <iostream>

#include "data.h"
#include "grandpotential.h"


int main() {
    std::cout << "Hello, This program has ome features." << std::endl;

    while (true) {
        std::cout
            << "      1. Grand potential calculation" << std::endl
            << "      q. Quit" << std::endl
            << "Enter your choice: ";

        char choice;
        std::cin >> choice;
        std::cout << std::string(60, '-') << std::endl;

        switch (choice) {
            case '1': {
                gp::run();
                return 0;
            };
            case 'q': {
                return 0;
            };
            default:
                std::cout << "Invalid choice. Please try again." << std::endl;
        }
    }
}


