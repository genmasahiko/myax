#include "cube.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <array>

void Cube::LoadCubeFile(std::string filename) {
    std::ifstream cubeFile(filename);

    if ( !cubeFile.is_open() ) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        exit(1);
    }

    std::string line;
    int lineNum = 0;

    while (std::getline(cubeFile, line)) {
        lineNum++;

        std::istringstream iss(line);

        if (lineNum == 3) {
            iss >> nat >> origin[0] >> origin[1] >> origin[2];
        } else if (lineNum == 4) {
            iss >> grid[0] >> stepw[0][0] >> stepw[0][1] >> stepw[0][2];
        } else if (lineNum == 5) {
            iss >> grid[1] >> stepw[1][0] >> stepw[1][1] >> stepw[1][2];
        } else if (lineNum == 6) {
            iss >> grid[2] >> stepw[2][0] >> stepw[2][1] >> stepw[2][2];
        } else if (lineNum > 6 && lineNum <= 6 + nat) {
            Atom atom;
            iss >> atom.atomic_number >> atom.charge >> atom.position[0] >> atom.position[1] >> atom.position[2];
            atoms.push_back(atom);
        } else if (lineNum >= 6 + nat) {
            break;
            // to be fixed
        }
    }

    // reflesh cubeFile: This feature will be fixed
    cubeFile.clear();
    cubeFile.seekg(0, std::ios::beg);
    for ( int i = 0; i < 6 + nat; i++ ) {
        std::getline(cubeFile, line);
    }

    // Read volumatic data as 1D array
    vol_1d.resize(grid[0] * grid[1] * grid[2]);
    for (int i = 0; i < grid[0] * grid[1] * grid[2]; i++) {
        cubeFile >> vol_1d[i];
    }

    for (int i = 0; i < 3; i++) {
        alat[i] = sqrt(stepw[i][0] * stepw[i][0] + stepw[i][1] * stepw[i][1] + stepw[i][2] * stepw[i][2]);
    }

}

//void Cube::ShiftVolumaticData() {
//
//    if ( lset ) {
//
//        for ( int i = 0; i < data.grid[0]; i++ ) {
//            for ( int j = 0; j < data.grid[1]; j++ ) {
//                for ( int k = 0; k < data.grid[2]; k++ ) {
//                    volnew[i][j][k] = data.vol[ ( i + data.grid[0] / 2 ) % data.grid[0] ][ ( j + data.grid[1] / 2 ) % data.grid[1] ][ ( k + data.grid[2] / 2 ) % data.grid[2] ];
//                }
//            }
//        }
//
//    } else {
//
//        for ( int i = 0; i < data.grid[0]; i++ ) {
//            for ( int j = 0; j < data.grid[1]; j++ ) {
//                for ( int k = 0; k < data.grid[2]; k++ ) {
//                    volnew[i][j][k] = data.vol[i][j][ ( k + data.grid[2] / 2 ) % data.grid[2] ];
//                }
//            }
//        }
//
//    }
//
//    volnew.swap(data.vol);
//
//}
//

std::array<int, 3> Cube::Getgrid() { return grid; }

std::array<int, 3> Cube::Getorigin() { return origin; }

std::array<std::array<double, 3>, 3> Cube::Getstepw() { return stepw; }

std::array<double, 3> Cube::Getalat() { return alat; }

std::vector<double> Cube::Getvol_1d() { return vol_1d; }
