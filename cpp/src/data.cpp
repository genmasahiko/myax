#include "data.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <any>
#include <variant>
#include <array>
#include <iomanip>

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

void Data::SetParam( std::unordered_map< std::string, std::unordered_map< std::string, std::string > > param ) {
    params = param;
}

void Data::SetPseudo( std::string symbol, float mass, std::string name) {
    Pseudo new_pseudo(symbol, mass, name);
    pseudos_.push_back(new_pseudo);
}

void Data::SetKpoints( const std::string& style, const std::array<int, 3>& nk, const std::array<int, 3>& sk) {
    kpoints.style = style;
    kpoints.nk = nk;
    kpoints.sk = sk;
}

void Data::SetCell( const std::string& style, const std::vector< std::vector<double> >& v) {
    cell.style = style;
    cell.v = v;
}

void Data::SetParam( const std::string& section, const std::string& key, const std::string& value ) {
    params[section][key] = value;
}

void Data::SetKpoints4bands(
        const std::string& style,
        const int& nks,
        const std::vector< std::vector<float> >& xk,
        const std::vector<float>& wk
        ) {
    kpoints.style = style;
    kpoints.nks = nks;
    kpoints.xk = xk;
    kpoints.wk = wk;
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

std::string Data::GetParam( const std::string& section, const std::string& key ) {
    if (params.find(section) != params.end() && params[section].find(key) != params[section].end()) {
        return params[section][key];
    } else {
        std::cerr << "Error: Parameter not found: " << section << " -> " << key << std::endl;
        return "";
    }
}

std::string GetStringParam( std::string line ) {
    int start = line.find("'") + 1;
    int end = line.find("'", start);
    return line.substr(start, end - start);
}

std::string Trim( std::string str ) {
    str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
    str.erase(str.find_last_not_of(" \t\n\r\f\v") + 1);
    return str;
}

std::string ExtractInnerString( std::string str ) {
    // Find the first '(', '{', or '['
    size_t start = str.find_first_of("({[");
    if (start == std::string::npos) {
        return ""; // No opening bracket found
    }

    // Find the matching ')', '}', or ']'
    size_t end = str.find_first_of(")}]", start);
    if (end == std::string::npos) {
        return ""; // No closing bracket found
    }
    // Extract the substring between the brackets
    std::string inner = str.substr(start + 1, end - start - 1);
    return inner;
}

Data Data::ReadInfile(std::ifstream &file) {

    Data in;
    std::string line;

    //using ParamType = std::variant< int, double, std::string >;

    std::unordered_map<std::string, std::unordered_map<std::string, std::string> > param;
    std::string calculation;
    std::string restartmode;


    while (std::getline(file, line)) {
        if (line.find("calculation") != std::string::npos) {
            param["control"]["calculation"] = GetStringParam(line);
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
    }

    std::string section;

    while (std::getline(file, line)) {
        line = Trim(line);
        if ( line.empty() || line == "/" ) continue;

        if ( line.find("&") != std::string::npos ) {
            section = line.substr(1);

        } else if (line.find("ATOMIC_SPECIES") != std::string::npos) {
            for ( int i = 0; i < std::stoi(param["system"]["ntyp"]) ; i++ ) {
                std::getline(file, line);
                std::istringstream iss(line);
                std::string symbol;
                float mass;
                std::string name;

                iss >> symbol >> mass >> name;

                in.SetPseudo(symbol, mass, name);
            }

        } else if (line.find("K_POINTS") != std::string::npos) {
            std::istringstream iss(line);
            std::string dum, style;
            iss >> dum >> style;

            style = ExtractInnerString(style);

            if ( style == "automatic" ) {
                std::getline(file, line);
                std::istringstream iss(line);
                std::array<int, 3> nk, sk;

                iss >> nk[0] >> nk[1] >> nk[2];
                iss >> sk[0] >> sk[1] >> sk[2];

                in.SetKpoints(style, nk, sk);
            }

        } else if (line.find("CELL_PARAMETERS") != std::string::npos) {
            std::istringstream iss(line);
            std::string dum, style;
            iss >> dum >> style;

            style = ExtractInnerString(style);

            std::vector< std::vector<double> > v;

            for (int i = 0; i < 3; i++) {
                std::getline(file, line);
                std::istringstream iss(line);
                std::vector<double> vec(3);
                iss >> vec[0] >> vec[1] >> vec[2];
                v.push_back(vec);
            }

            in.SetCell(style, v);

        } else if (line.find("ATOMIC_POSITIONS") != std::string::npos) {
            std::istringstream iss(line);
            std::string dum, style;
            iss >> dum >> style;

            style = ExtractInnerString(style);

            int nat = std::stoi(param["system"]["nat"]);

            for (int i = 0; i < nat; i++) {
                std::getline(file, line);
                std::istringstream iss(line);
                std::string symbol;
                std::vector<double> pos(3);
                std::vector<int> ifpos(3);

                iss >> symbol >> pos[0] >> pos[1] >> pos[2] >> ifpos[0] >> ifpos[1] >> ifpos[2];

                in.SetAtom(symbol, pos, ifpos);
            }

        } else {
            if (line.find("=") == std::string::npos) {
                continue;
            }
            std::string key = Trim( line.substr(0, line.find("=")) );
            std::string value = Trim( line.substr(line.find("=") + 1) );
            param[section][key] = value;
        }
        
    }

    in.SetParam(param);

    file.clear();
    file.seekg(0, std::ios::beg);

    return in;
};

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

void Data::WriteBandInfile (std::ofstream& file) {

    std::vector< std::string > section_order = {
        "control",
        "system",
        "electrons",
    };
    
    for ( const auto& section : section_order ) {
        file << "&" << section << std::endl;
        for ( const auto& [key, value] : params[section] ) {
            file << " " << key << " = " << value << std::endl;
        }
        file << "/" << std::endl;
    }

    file << "\nATOMIC_SPECIES" << std::endl;
    for ( const auto& pseudo : pseudos_ ) {
        file << std::setw(3) << std::left << pseudo.symbol
             << std::setw(5) << std::left << pseudo.mass
             << std::left << pseudo.name
             << std::endl;
    }

    file << "\nK_POINTS " << "(" << kpoints.style  << ")" << std::endl;
    file << kpoints.nks << std::endl;
    for ( int i = 0; i < kpoints.nks; i++ ) {
        file << std::fixed << std::setprecision(5) << std::setw(10) << std::right << kpoints.xk[i][0]
             << std::fixed << std::setprecision(5) << std::setw(10) << std::right << kpoints.xk[i][1]
             << std::fixed << std::setprecision(5) << std::setw(10) << std::right << kpoints.xk[i][2]
             << std::fixed << std::setprecision(5) << std::setw(10) << std::right << kpoints.wk[i]
             << std::endl;
    }

    file << "\nCELL_PARAMETERS " << "(" << cell.style  << ")" << std::endl;
    for ( int i = 0; i < 3; i++ ) {
        file << std::fixed << std::setprecision(9) << std::setw(14) << std::right << cell.v[i][0]
             << std::fixed << std::setprecision(9) << std::setw(14) << std::right << cell.v[i][1]
             << std::fixed << std::setprecision(9) << std::setw(14) << std::right << cell.v[i][2]
             << std::endl;
    }

    std::string nat = params["system"]["nat"];

    file << "\nATOMIC_POSITIONS " << "(" << "angstrom"  << ")" << std::endl;
    for ( int i = 0; i < std::stoi(nat); i++ ) {
        file << std::setw(3) << std::left << atoms[i].symbol
             << std::fixed << std::setprecision(9) << std::setw(14) << std::right << atoms[i].pos[0]
             << std::fixed << std::setprecision(9) << std::setw(14) << std::right << atoms[i].pos[1]
             << std::fixed << std::setprecision(9) << std::setw(14) << std::right << atoms[i].pos[2]
             << std::endl;
    }
    
}
