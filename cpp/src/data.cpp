#include "data.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <any>

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

void Data::SetCalculation(std::string calc) {
    calculation = calc;
}

void Data::SetRestartmode(std::string mode) {
    restartmode = mode;
}

void Data::SetParam( std::unordered_map<std::string, std::any> param ) {
    params = param;
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

std::string Data::GetCalculation() {
    return calculation;
}

std::string Data::GetRestartmode() {
    return restartmode;
}

// It will be revised and used.
//
//template <typename T> 
//T Data::GetParam( std::string key ) {
//    T value;
//    if ( decltype(value) == params[key].type() ) {
//        return std::any_cast<T>(params[key]);
//    } else {
//        std::cerr << "Error: The type of the parameter is not " << typeid(T).name() << std::endl;
//    }
//}

std::string GetStringParam( std::string line ) {
    int start = line.find("'") + 1;
    int end = line.find("'", start);
    return line.substr(start, end - start);
}

Data Data::ReadInfile(std::ifstream &file) {

    Data in;
    std::string line;

    std::unordered_map<std::string, std::any> param;
    std::string calculation;
    std::string restartmode;

    while (std::getline(file, line)) {
        if (line.find("calculation") != std::string::npos) {
            param["calculation"] = GetStringParam(line);
            break;
        }

        if (line.find("BEGIN_PATH_INPUT") != std::string::npos ||
            line.find("begin_path_input") != std::string::npos) {

            calculation = "neb";
            break;
        }
    }

    file.clear();
    file.seekg(0, std::ios::beg);

    if ( calculation == "neb" ) {
        while (std::getline(file, line)) {

            if ( line.find("restart_mode") != std::string::npos ) {
                restartmode = GetStringParam(line);
            }

        }
    } else {
        std::cerr << "Error: Function ReadInfile do not ready for calculation except NEB."
            << std::endl;
    }

    //in.SetCalculation(calculation);
    in.SetParam(param);
    in.SetRestartmode(restartmode);

    file.clear();
    file.seekg(0, std::ios::beg);

    return in;
}

Data Data::ReadOutfile(std::ifstream &file) {

    Data out;
    std::string line;
    std::string program;
    bool lfinish = false;

    int nat = 0;
    int ntyp = 0;
    std::string calculation = "scf";
    std::vector< std::string > atom_symbol;
    std::vector< std::vector<double> > atom_pos;
    std::vector< std::vector<std::string> > atom_ifpos_str;
    std::vector< std::vector<int> > atom_ifpos;
    std::vector< std::vector<double> > atom_force;

    while (std::getline(file, line)) {
        if (line.find("PWSCF") != std::string::npos) {
            program = "pw";
        } else if (line.find("NEB") != std::string::npos) {
            program = "neb";
            calculation = "neb";
        }
    }

    if ( program == "pw" ) {

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
                calculation = "relax";

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
        }

        if (!lfinish && calculation == "relax") {
            std::cerr << "Error: Optimization did not finish!" << std::endl;
        }

    } else if ( program == "neb" ) {

    }

    file.clear();
    file.seekg(0, std::ios::beg);

    return out;
}

