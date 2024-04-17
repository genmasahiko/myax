#include "data.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

void Data::SetNat(int n) { 
    nat = n;
}

void Data::SetNtyp(int n) { 
    ntyp = n;
}

void Data::SetAtom( std::string symbol, std::vector<double> pos, std::vector<int> ifpos, std::vector<double> force) {
    Atom new_atom;
    new_atom.symbol = symbol;
    new_atom.pos = pos;
    new_atom.ifpos = ifpos;
    new_atom.force = force;

    atoms.push_back(new_atom);
}

int Data::GetNat() { 
    return nat;
}

std::vector<double> Data::GetForce(int i) { 
    return atoms[i].force;
}

std::vector<int> Data::GetIfpos(int i) { 
    return atoms[i].ifpos;
}

Data Data::ReadOutfile(std::ifstream &file) {

    Data out;
    std::string line;
    bool lfinish = false;

    int nat, ntyp;
    std::vector< std::string > atom_symbol;
    std::vector< std::vector<double> > atom_pos;
    std::vector< std::vector<std::string> > atom_ifpos_str;
    std::vector< std::vector<int> > atom_ifpos;
    std::vector< std::vector<double> > atom_force;


    while (std::getline(file, line)) {

        if (line.find("number of atoms/cell      =") != std::string::npos) {
            line = line.substr( line.find("=") + 1 );
            nat = std::stoi(line);

            atom_symbol.resize(nat);
            atom_pos.resize(nat);
            atom_ifpos_str.resize(nat);
            atom_ifpos.resize(nat);
            atom_force.resize(nat);
            for (int i = 0; i < nat; i++) {
                atom_pos[i].resize(3);
                atom_ifpos_str[i].resize(3);
                atom_ifpos[i].resize(3);
                atom_force[i].resize(3);
            }
        }


        if (line.find("number of atomic types    =") != std::string::npos) {
            line = line.substr( line.find("=") + 1 );
            ntyp = std::stoi(line);
        }


        if (line.find("Forces acting on atoms") != std::string::npos) {
            std::getline(file, line); // skip blank line
            for (int i = 0; i < nat; i++) {
                std::getline(file, line);
                line = line.substr( line.find("=") + 1 );
                std::istringstream iss(line);

                iss >> atom_force[i][0] >> atom_force[i][1] >> atom_force[i][2];
            }
        }

        if (line.find("ATOMIC_POSITIONS") != std::string::npos) {
            for (int i = 0; i < nat; i++) {
                std::getline(file, line);
                std::istringstream iss(line);

                iss
                    >> atom_symbol[i]
                    >> atom_pos[i][0] >> atom_pos[i][1] >> atom_pos[i][2]
                    >> atom_ifpos_str[i][0] >> atom_ifpos_str[i][1] >> atom_ifpos_str[i][2];

                for (int j = 0; j < 3; j++) {
                    if (atom_ifpos_str[i][j] == "0") {
                        atom_ifpos[i][j] = 0;
                    } else if ( atom_ifpos_str[i][j] == "1" || atom_ifpos_str[i][j] == "" ) {
                        atom_ifpos[i][j] = 1;
                    }
                }
            }
        }

        if (line.find("End of BFGS Geometry Optimization") != std::string::npos || 
            line.find("End of damped dynamics calculation") != std::string::npos) {

            lfinish = true;

        }
    }

    out.SetNat(nat);
    out.SetNtyp(ntyp);
    for (int i = 0; i < nat; i++) {
        out.SetAtom(atom_symbol[i], atom_pos[i], atom_ifpos[i], atom_force[i]);
    };

    if (!lfinish) {
        std::cerr << "Error: Optimization did not finish!" << std::endl;
    };


    return out;
}
