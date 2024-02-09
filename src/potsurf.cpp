#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include <iomanip>

class Indata {
public:
    std::vector<std::vector<double>> latvec;
    int surf_size;

    class Atom {
        public:
        std::string name;
        std::vector<double> pos;
    };

    std::vector<Indata::Atom> atoms;

    int FindAtom(std::string name, int order = 1) {
        int count = 0;
        for (int i = 0; i < atoms.size(); i++) {
            if (atoms[i].name == name) {
                count++;
                if (count == order) {
                    return i;
                };
            };
        };
        return -1;
    };

    Indata() : surf_size(0) {};

};

int LoadInitialData(std::string filename, Indata &data);
int Shifth2o(int nearest, int index, std::string filename, Indata &data);
int Rotateh2o(std::string filename, Indata &data);

int main() {

    std::cout << "Initial strucure QE input file: " << std::endl;
    std::string filename;
    std::cin >> filename;

    std::cout << "There are three options" << std::endl;
    std::cout << "1: Toward 1st nearest" << std::endl;
    std::cout << "2: Toward 2nd nearest" << std::endl;
    std::cout << "3: Rotate" << std::endl;
    int option;
    std::cin >> option;

    Indata data;
    LoadInitialData(filename, data);

    if ( option < 3 ) {

        for (int i = 0; i < 6; i++) {

            Shifth2o(option, i, filename, data);

        };

    } else if ( option == 3 ) {

        Rotateh2o(filename, data);

    };

    return 0;

}

int LoadInitialData(std::string filename, Indata &data) {

    std::ifstream infile(filename);
    std::string line;

    if ( !infile.is_open() ) {
        std::cout << "Error: file not exist" << std::endl;
        return 1;
    };

    while (std::getline(infile, line)) {

        if (line.find("nat") != std::string::npos) {

            std::istringstream iss(line);
            std::string temp;
            int nat;

            iss >> temp >> temp >> nat;

            data.surf_size = sqrt( (nat - 3) / 4 );
        };

        if (line.find("CELL_PARAM") != std::string::npos) {

            for (int i = 0; i < 3; i++) {
                std::getline(infile, line);
                std::istringstream iss(line);

                double a, b, c;
                iss >> a >> b >> c;

                std::vector<double> temp;
                temp.push_back(a);
                temp.push_back(b);
                temp.push_back(c);

                data.latvec.push_back(temp);

            };

        };

        if (line.find("ATOMIC_POSITIONS") != std::string::npos) {

            while (std::getline(infile, line)) {

                if ( line.find("H") != std::string::npos || line.find("O") != std::string::npos ) {

                    std::istringstream iss(line);
                    Indata::Atom atom;
                    double x, y, z;

                    iss >> atom.name >> x >> y >> z;

                    atom.pos.push_back(x);
                    atom.pos.push_back(y);
                    atom.pos.push_back(z);

                    data.atoms.push_back(atom);

                };
            };
        };
    };

    infile.close();

    return 0;
};

int Shifth2o(int nearest, int index, std::string filename, Indata &data) {

    double const pi = 3.14159265359;
    std::string line;
    std::ifstream infile(filename);

    if ( !infile.is_open() ) {
        std::cout << "Error: file not exist" << std::endl;
        return 1;
    };

    std::string dirname = std::to_string(index) + "ori";
    std::string command1 = "mkdir " + dirname;
    std::system(command1.c_str());

    for (int i = 0; i < 11; i++) {

        std::string command2 = "mkdir " + dirname + "/" + std::to_string(i);
        std::system(command2.c_str());

    };

    for (int i = 0; i < 11; i++) {

        std::ofstream outfile( dirname + "/" + std::to_string(i) + "/in" );
        int j = 0;

        while (std::getline(infile, line)) {

            outfile << line << std::endl;

            if ( line.find("ATOMIC_POSITIONS") != std::string::npos ) {
                break;
            };

        };

        if ( nearest == 1 ) {

            while (std::getline(infile, line)) {

                if ( line.find("H") != std::string::npos || line.find("O") != std::string::npos ) {

                    double xnew = data.atoms[j].pos[0] + data.latvec[0][0] / data.surf_size * i * 0.1 * cos( index * pi / 3 );
                    double ynew = data.atoms[j].pos[1] + data.latvec[0][0] / data.surf_size * i * 0.1 * sin( index * pi / 3 );

                    outfile << std::left << std::setw(6) << data.atoms[j].name
                            << std::right << std::setw(14) << std::fixed << std::setprecision(9) << xnew
                            << std::right << std::setw(14) << std::fixed << std::setprecision(9) << ynew
                            << std::right << std::setw(14) << std::fixed << std::setprecision(9) << data.atoms[j].pos[2] 
                            << std::right << std::setw(5) << "0"
                            << std::right << std::setw(4) << "0"
                            << std::right << std::setw(4) << "1"
                            << std::endl;

                    j++;

                } else {

                    outfile << line << std::endl;

                };

            };

        } else if ( nearest == 2 ) {

            while ( std::getline(infile, line) ) {

                if ( line.find("H") != std::string::npos || line.find("O") != std::string::npos ) {

                    double xnew = data.atoms[j].pos[0] + sqrt(3) * data.latvec[0][0] / data.surf_size * i * 0.1 * cos( (2*index + 1) * pi / 6 );
                    double ynew = data.atoms[j].pos[1] + sqrt(3) * data.latvec[0][0] / data.surf_size * i * 0.1 * sin( (2*index + 1) * pi / 6 );

                    outfile << std::left << std::setw(6) << data.atoms[j].name
                            << std::right << std::setw(14) << std::fixed << std::setprecision(9) << xnew
                            << std::right << std::setw(14) << std::fixed << std::setprecision(9) << ynew
                            << std::right << std::setw(14) << std::fixed << std::setprecision(9) << data.atoms[j].pos[2] 
                            << std::right << std::setw(5) << "0"
                            << std::right << std::setw(4) << "0"
                            << std::right << std::setw(4) << "1"
                            << std::endl;

                    j++;

                } else {

                    outfile << line << std::endl;

                };

            };

        };

        infile.clear(); infile.seekg(0, std::ios::beg);

    };

    return 0;

};

int Rotateh2o(std::string filename, Indata &data) {

    double const pi = 3.14159265359;
    std::string line;
    std::ifstream infile(filename);

    if ( !infile.is_open() ) {
        std::cout << "Error: file not exist" << std::endl;
        return 1;
    };

    std::string dirname = "rotate";
    std::string command1 = "mkdir " + dirname;
    std::system(command1.c_str());

    for (int i = 0; i < 11; i++) {

        std::string command2 = "mkdir " + dirname + "/" + std::to_string(i);
        std::system(command2.c_str());

    };

    for (int i = 0; i < 10; i++) {

        std::ofstream outfile( dirname + "/" + std::to_string(i) + "/in" );

        while (std::getline(infile, line)) {

            outfile << line << std::endl;

            if ( line.find("ATOMIC_POSITIONS") != std::string::npos ) {
                break;
            };

        };

        int j = 0;
        int oindex = data.FindAtom("O");

        while ( std::getline(infile, line) ) {

            if ( line.find("H") != std::string::npos ) {

                j++;
                int hindex = data.FindAtom("H", j);

                double xnew = ( data.atoms[hindex].pos[0] - data.atoms[oindex].pos[0] ) * cos( i * pi * 0.2 ) 
                            - ( data.atoms[hindex].pos[1] - data.atoms[oindex].pos[1] ) * sin( i * pi * 0.2 ) 
                            + data.atoms[oindex].pos[0];

                double ynew = ( data.atoms[hindex].pos[0] - data.atoms[oindex].pos[0] ) * sin( i * pi * 0.2 ) 
                            + ( data.atoms[hindex].pos[1] - data.atoms[oindex].pos[1] ) * cos( i * pi * 0.2 ) 
                            + data.atoms[oindex].pos[1];

                outfile << std::left << std::setw(6) << data.atoms[j].name
                        << std::right << std::setw(14) << std::fixed << std::setprecision(9) << xnew
                        << std::right << std::setw(14) << std::fixed << std::setprecision(9) << ynew
                        << std::right << std::setw(14) << std::fixed << std::setprecision(9) << data.atoms[j].pos[2] 
                        << std::right << std::setw(5) << "0"
                        << std::right << std::setw(4) << "0"
                        << std::right << std::setw(4) << "1"
                        << std::endl;

            } else {

                outfile << line << std::endl;

            };

        };

        infile.clear(); infile.seekg(0, std::ios::beg);

    };

    return 0;

};
