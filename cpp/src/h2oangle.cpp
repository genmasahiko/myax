#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "data.h"
#include "h2oangle.h"

namespace angle {
    double AngleBetween(const std::vector<double>& v1, const std::vector<double>& v2) {
        double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        double abs_v1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        double abs_v2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
        return std::acos(dot / (abs_v1 * abs_v2)) * 180.0 / M_PI; // Convert to degrees
    }

    std::vector<double> VectorMiddle(const std::vector<double>& h1, const std::vector<double>& h2) {
        return { (h1[0] + h2[0]) / 2.0, (h1[1] + h2[1]) / 2.0, (h1[2] + h2[2]) / 2.0 };
    }

    void WriteOutput(double dipole_angle, double dipole_angle_xy, double hoh_angle, std::ofstream& output_file) {
        output_file << "Dipole vs z-axis" << "    "
                    << "Dipole vs x-axis in xy-plane" << "    "
                    << "HOH angle" << " "
                    << "[degrees]" << std::endl;
        output_file << "----------------------------------------" << std::endl;
        output_file << dipole_angle << "    "
                    << dipole_angle_xy << "    "
                    << hoh_angle << std::endl;
    }

    int run() {
        std::cout << "Running H2O angle calculation module..." << std::endl;

        // Reat the poscar file
        std::cout << "Enter the POSCAR file name: ";
        std::string poscar_file_name;
        std::cin >> poscar_file_name;

        // Open the POSCAR file
        std::ifstream poscar_file(poscar_file_name);
        if (!poscar_file) {
            std::cerr << "Error opening POSCAR file: " << poscar_file_name << std::endl;
            return 1;
        }

        // Create Data instance and read the POSCAR data
        Data poscar;
        poscar = poscar.ReadPoscar(poscar_file);

        // Get the oxygen and hydrogen atoms index
        int o_index = poscar.GetAtomIndex("O");
        int h1_index = poscar.GetAtomIndex("H", 1);
        int h2_index = poscar.GetAtomIndex("H", 2);

        // Get the oxygen and hydrogen positions
        std::vector<double> o_pos = poscar.GetAtom(o_index);
        std::vector<double> h1_pos = poscar.GetAtom(h1_index);
        std::vector<double> h2_pos = poscar.GetAtom(h2_index);

        // Calculate the dipole angle against the z-axis
        std::vector<double> z_axis = { 0.0, 0.0, 1.0 };
        std::vector<double> h_middle = VectorMiddle(h1_pos, h2_pos);
        std::vector<double> dipole = { h_middle[0] - o_pos[0], h_middle[1] - o_pos[1], h_middle[2] - o_pos[2] };
        double dipole_angle = AngleBetween(dipole, z_axis);

        // Calculate the dipole angle against the x-axis in the xy-plane
        std::vector<double> dipole_inxy = { dipole[0], dipole[1], 0.0 };
        std::vector<double> x_axis = { 1.0, 0.0, 0.0 };
        double dipole_angle_xy = AngleBetween(dipole_inxy, x_axis);

        // Calculate the HOH angle
        std::vector<double> oh1 = { h1_pos[0] - o_pos[0], h1_pos[1] - o_pos[1], h1_pos[2] - o_pos[2] };
        std::vector<double> oh2 = { h2_pos[0] - o_pos[0], h2_pos[1] - o_pos[1], h2_pos[2] - o_pos[2] };
        double hoh_angle = AngleBetween(oh1, oh2);

        // Output the results
        std::string output_file_name = "angle.dat";
        if (std::ifstream(output_file_name)) {
            std::cout << "Output file already exists. Overwrite? (y/n): ";
            char choice;
            std::cin >> choice;
            if (choice != 'y') {
                std::cout << "Exiting without writing output." << std::endl;
                return 0;
            }
        }
        std::ofstream output_file(output_file_name);
        if (!output_file) {
            std::cerr << "Error opening output file: " << output_file_name << std::endl;
            return 1;
        }
        WriteOutput(dipole_angle, dipole_angle_xy, hoh_angle, output_file);

        return 0;
    }
}
