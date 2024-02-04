#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <filesystem>
#include <cmath>
#include <cstdlib>

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

    Indata() : surf_size(0) {};

};

int LoadInitialData(std::string filename, Indata &data);
int Shifth2o(int index, std::string filename, Indata &data);

int main() {

    std::cout << "Initial strucure QE input file: " << std::endl;
    std::string filename;
    std::cin >> filename;

    Indata data;

    LoadInitialData(filename, data);

    for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
            std::cout << data.latvec[i][j] << " ";
        };
        std::cout << std::endl;
    };

    for (int i = 0; i < 6; i++) {
        Shifth2o(i, filename, data);
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

int Shifth2o(int index, std::string filename, Indata &data) {

    double const pi = 3.14159265359;
    std::string line;
    std::ifstream infile(filename);

    if ( !infile.is_open() ) {
        std::cout << "Error: file not exist" << std::endl;
        return 1;
    };

    std::string dirname = std::to_string(index) + "ori";

    if ( !std::filesystem::exists( dirname ) ) {

        for (int i = 0; i < 11; i++) {

            std::filesystem::create_directory( dirname );
            std::filesystem::create_directory( dirname + "/" + std::to_string(i) );

        };

    }

    for (int i = 0; i < 11; i++) {

        std::ofstream outfile( dirname + "/" + std::to_string(i) + "/in" );
        int j = 0;

        while (std::getline(infile, line)) {

            outfile << line << std::endl;

            if ( line.find("ATOMIC_POSITIONS") != std::string::npos ) {
                break;
            };

        };

        while (std::getline(infile, line)) {

            if ( line.find("H") != std::string::npos || line.find("O") != std::string::npos ) {

                outfile << std::left << std::setw(6) << data.atoms[j].name
                        << std::right << std::setw(14) << std::fixed << std::setprecision(9) << data.atoms[j].pos[0] + data.latvec[0][0] / data.surf_size * i * 0.1 * cos( index * pi / 3 )
                        << std::right << std::setw(14) << std::fixed << std::setprecision(9) << data.atoms[j].pos[1] + data.latvec[0][0] / data.surf_size * i * 0.1 * sin( index * pi / 3 )
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

        infile.clear(); infile.seekg(0, std::ios::beg);

    };

    return 0;

}
