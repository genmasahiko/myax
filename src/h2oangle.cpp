#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <unordered_map>

struct Atom {
    std::string type = "";
    int count = 0;
    double x = 0.0, y = 0.0, z = 0.0;
};

std::vector<Atom> readPOSCAR(const std::string& filename) {
    std::vector<Atom> atoms;
    std::ifstream file(filename);
    std::string line;

    std::vector<std::string> atomTypes;
    std::vector<int> atomNums;

    if (file.is_open()) {
        for (int i = 0; i < 5; ++i) {
            std::getline(file, line);
        }

        std::getline(file, line);
        std::istringstream iss(line);

        std::string atomType;
        while (iss >> atomType) {
            atomTypes.push_back(atomType);
        }

        std::getline(file, line);
        iss.clear();
        iss.str(line);

        int atomNum;
        while (iss >> atomNum) {
            atomNums.push_back(atomNum);
        }

        std::getline(file, line);

        int i = 0;
        int j = 0;
        while ( std::getline(file, line) ) {
            std::istringstream iss(line);

            Atom atom;
            iss >> atom.x >> atom.y >> atom.z;

            atom.type = atomTypes[j];
            atom.count = i + 1;

            atoms.push_back(atom);

            if (i == atomNums[j] - 1) {
                i = 0;
                j++;
            } else {
                i++;
            }

        }
    }

    return atoms;
}

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const std::vector<double>& a) {
    return sqrt(dotProduct(a, a));
}

double angleBetween(const std::vector<double>& a, const std::vector<double>& b) {
    double dot = dotProduct(a, b);
    double norms = norm(a) * norm(b);
    return acos(dot / norms) * 180.0 / M_PI; // ラジアンから度へ変換
}

Atom midlle_point(const Atom& atom1, const Atom& atom2) {
    Atom midpoint;

    midpoint.type = "M";
    midpoint.x = ( atom1.x + atom2.x ) * 0.5;
    midpoint.y = ( atom1.y + atom2.y ) * 0.5;
    midpoint.z = ( atom1.z + atom2.z ) * 0.5;

    return midpoint;
}

int main( int argc, char* argv[]) {
    std::string filename;

    if ( argc > 1 ) {
        filename = argv[1];
    } else {
        std::cout << "Enter the name of the POSCAR file: ";
        std::cin >> filename;
    }

    std::vector<Atom> atoms = readPOSCAR(filename);

    Atom hydrogen1, hydrogen2, oxygen;

    for (const auto& atom : atoms) {
        if (atom.type == "O") {
            oxygen = atom;
        } else if (atom.type == "H") {
            if (atom.count == 1) {
                hydrogen1 = atom;
            } else {
                hydrogen2 = atom;
            }
        }
    }

    Atom midpoint = midlle_point(hydrogen1, hydrogen2);

    std::vector<double> h2oDipoleDirection = { midpoint.x - oxygen.x, midpoint.y - oxygen.y, midpoint.z - oxygen.z };
    std::vector<double> zAxis = { 0.0, 0.0, 1.0 };

    double angle = angleBetween(h2oDipoleDirection, zAxis);
    std::cout << angle << std::endl;

    return 0;
}

