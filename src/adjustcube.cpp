// Purpose: Adjust the gaussian cube file to a human-readable format
// when read by VESTA.
// Particulary, if you use ESM method in Quantum ESPRESSO,
// the slab model is centered at the origin of the cell.
// This program moves the slab to the center of the cell and also shift the
// volumatic data in cube file.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

class Data {
public:
    int nat;
    double origin[3];
    int grid[3];
    double stepw[3][3];

    class Atom {
        public:
        int atomnum;
        double charge;
        double pos[3];
    };

    std::vector<Data::Atom> atoms;

    int FindAtom(int number, int order = 1) {
        int count = 0;
        for (int i = 0; i < atoms.size(); i++) {
            if (atoms[i].atomnum == number) {
                count++;
                if (count == order) {
                    return i;
                };
            };
        };
        return -1;
    };
};

void LoadCubeFile(std::string filename, Data& data); 
void ShiftVolumaticData(std::string filename, bool lset, Data data); 
void ShiftAtomicData(std::string filename, bool lset, Data& data);
void WriteCubeFile(std::string filename, Data data);

int main(int argc, char* argv[]) {

    std::string filename;

    if (argc == 1) {
        std::cout << "filename?" << std::endl;
        std::cin >> filename;
    }
    else if (argc == 2) {
        filename = argv[1];
    }
    else if (argc > 2) {
        std::cerr << "Error: Invalid number of arguments" << std::endl;
        return 1;
    }

    std::cout << "Set H2O molecule at..." << std::endl;
    std::cout << "1. Center of the surface" << std::endl;
    std::cout << "2. Do not change" << std::endl;
    int tmp;
    std::cin >> tmp;

    bool lset_center;
    if (tmp == 1) {
        lset_center = true;
    }
    else if (tmp == 2) {
        lset_center = false;
    }

    // Load cube file and store data in Data class
    Data data;
    LoadCubeFile(filename, data);

    // Shift volumatic data and write to a temporary file named "tmp.cube"
    ShiftVolumaticData(filename, lset_center, data);
    ShiftAtomicData(filename, lset_center, data);

    // Write the shifted data to the original file
    // and delete the temporary file
    // New file is named filename_shifted.cube
    
    std::string newfilename = filename.substr(0, filename.find(".cube")) + "_shifted.cube";
    WriteCubeFile(newfilename, data);

    return 0;
};

void LoadCubeFile(std::string filename, Data& data) {
    std::ifstream cubeFile(filename);

    if (cubeFile.is_open()) {
        std::string line;
        int lineNum = 0;
        while (std::getline(cubeFile, line)) {
            lineNum++;

            std::istringstream iss(line);

            if (lineNum == 3) {
                iss >> data.nat >> data.origin[0] >> 
                    data.origin[1] >> data.origin[2];
            }
            else if (lineNum == 4) {
                iss >> data.grid[0] >> data.stepw[0][0] >> 
                    data.stepw[0][1] >> data.stepw[0][2];
            }
            else if (lineNum == 5) {
                iss >> data.grid[1] >> data.stepw[1][0] >> 
                    data.stepw[1][1] >> data.stepw[1][2];
            }
            else if (lineNum == 6) {
                iss >> data.grid[2] >> data.stepw[2][0] >> 
                    data.stepw[2][1] >> data.stepw[2][2];
            }
            else if (lineNum > 6 && lineNum <= 6 + data.nat) {
                Data::Atom atom;
                iss >> atom.atomnum >> atom.charge >> 
                    atom.pos[0] >> atom.pos[1] >> atom.pos[2];
                data.atoms.push_back(atom);
            }
            else if (lineNum >= 6 + data.nat) {
                break;
            }
        }
    }
    else {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
    }
}

void ShiftVolumaticData(std::string filename, bool lset, Data data) {

    // Shift z-axis -> x-axis and y-axis

    // for z-axis
    std::ifstream cubeFile(filename);
    std::ofstream tmpFile("tmp.cube");

    if (cubeFile.is_open() && tmpFile.is_open()) {

        int HeaderLineEnd = 6;
        int VolLineStart = HeaderLineEnd + data.nat;
        std::string line;

        for ( int i = 0; i < VolLineStart; i++ ) {
            std::getline(cubeFile, line);
        }

        int NvolSurf = data.grid[0] * data.grid[1];
        std::vector<double> vol( data.grid[2] );
        std::vector<double> volnew( data.grid[2] );

        for ( int i = 0; i < NvolSurf; i++ ) {
            for ( int j = 0; j < data.grid[2]; j++ ) {

                // In c++, the number of lines read by getline() is stored in a filestream object,
                // so we can retrieve the data as follows:
                cubeFile >> vol[j];
            }

            for ( int j = 0; j < data.grid[2]; j++ ) {
                
                volnew[j] = vol[ ( j + data.grid[2] / 2 ) % data.grid[2] ];

                // In fact, the volumatic data format is like 0.xxxxxE-xx (precision 5)
                // But c++ cannot write the data in this format;
                // c++ has to write the data like x.xxxxxE-xx
                // So, I have to set the precision to 4.
                tmpFile << std::right << std::scientific << std::setprecision(4) << std::setw(13) << volnew[j];

                // Break the line every 6 data and at the end of one z-block.
                if ( ( j + 1 ) % 6 == 0 || ( j + 1 ) == data.grid[2] ) {
                    tmpFile << std::endl;
                }
            }
        }

        cubeFile.close();
        tmpFile.close();

    } else {

        std::cerr << "Error: Cannot open file: " << filename << std::endl;
    };

    // for x-axis and y-axis
    if ( lset ) { 

        std::ifstream tmpFile("tmp.cube");

        if ( !tmpFile.is_open() ) {
            std::cerr << "Error: Cannot open file: tmp.cube or tmp2.cube" << std::endl;
        }

        int nblock = data.grid[0] * data.grid[1];
        int nline = ( data.grid[2] / 6 + 1 ) * nblock;

        std::cout << nblock << std::endl;

        std::string line;
        std::vector<std::string> block( nblock );

        for ( int i = 0; i < nline; i++ ) {
            std::getline(tmpFile, line);
            block[ i / ( data.grid[2] / 6 + 1 ) ] += line + "\n";
        };

        tmpFile.close();
        std::ofstream tmpFile2("tmp.cube");

        for ( int i = 0; i < nblock; i++ ) {
            tmpFile2 << block[ ( i - 10000 + nblock / 2 ) % nblock ];
        };

        tmpFile2.close();

    };

}

void ShiftAtomicData(std::string filename, bool lset, Data& data) {

    std::ifstream cubeFile(filename);

    // Shift toward c_axis

    double alat_c = data.grid[2] * data.stepw[2][2];
    double shift_c = alat_c / 2;

    for ( int i = 0; i < data.nat; i++ ) {
        data.atoms[i].pos[2] = fmod( data.atoms[i].pos[2] + shift_c, alat_c );
    }

    // Shift toward a_axis and b_axis

    double alat_surf = data.grid[0] * data.stepw[0][0];
    double shift_a = alat_surf * 0.5;
    double shift_bx = data.grid[1] * data.stepw[1][0] * 0.5;
    double shift_by = data.grid[1] * data.stepw[1][1] * 0.5;

    double slope_b = data.stepw[1][1] / data.stepw[1][0];
    double hypo = sqrt( pow( data.stepw[1][0], 2 ) + pow( data.stepw[1][1], 2 ) );
    double y2b = hypo / data.stepw[1][1];

    if ( lset ) {
        for ( int i = 0; i < data.nat; i++ ) {

            // for a_axis
            //data.atoms[i].pos[0] = fmod( data.atoms[i].pos[0] + shift_a, alat_surf );
            data.atoms[i].pos[0] = data.atoms[i].pos[0] + shift_a;

            if ( data.atoms[i].pos[0] > alat_surf + data.atoms[i].pos[1] / slope_b ) {
                data.atoms[i].pos[0] -= alat_surf;
            };

            // for b_axis
            data.atoms[i].pos[0] = data.atoms[i].pos[0] + shift_bx;
            data.atoms[i].pos[1] = data.atoms[i].pos[1] + shift_by;

            double b_coord = data.atoms[i].pos[1] * y2b;

            if ( b_coord > alat_surf ) {
                data.atoms[i].pos[0] -= shift_bx*2;
                data.atoms[i].pos[1] -= shift_by*2;
            };
        };
    }
}

void WriteCubeFile(std::string filename, Data data) {
    std::ofstream cubeFile(filename);

    if (cubeFile.is_open()) {
        cubeFile << "CUBE FILE" << std::endl;
        cubeFile << "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << std::endl;
        cubeFile << std::setw(5) << data.nat << std::setw(12) << data.origin[0] << 
            std::setw(12) << data.origin[1] << std::setw(12) << data.origin[2] << std::endl;
        cubeFile << std::setw(5) << data.grid[0] << std::setw(12) << data.stepw[0][0] << 
            std::setw(12) << data.stepw[0][1] << std::setw(12) << data.stepw[0][2] << std::endl;
        cubeFile << std::setw(5) << data.grid[1] << std::setw(12) << data.stepw[1][0] << 
            std::setw(12) << data.stepw[1][1] << std::setw(12) << data.stepw[1][2] << std::endl;
        cubeFile << std::setw(5) << data.grid[2] << std::setw(12) << data.stepw[2][0] << 
            std::setw(12) << data.stepw[2][1] << std::setw(12) << data.stepw[2][2] << std::endl;

        for ( int i = 0; i < data.nat; i++ ) {
            cubeFile << std::setw(5) << data.atoms[i].atomnum << std::setw(12) 
                << data.atoms[i].charge << std::setw(12) << data.atoms[i].pos[0] 
                << std::setw(12) << data.atoms[i].pos[1] << std::setw(12)
                << data.atoms[i].pos[2] << std::endl;
        }
    }
    else {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
    }

    std::ifstream tmpFile("tmp.cube");

    if ( tmpFile.is_open() ) {
        std::string line;
        while ( std::getline(tmpFile, line) ) {
            cubeFile << line << std::endl;
        }
    }
    else {
        std::cerr << "Error: Cannot open file: tmp.cube" << std::endl;
    }

    cubeFile.close();
    tmpFile.close();

    //std::remove("tmp.cube");
}
