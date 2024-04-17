#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include "data.h"

int main() {

    std::cout << "This program checks the convergence of relax" << std::endl;
    std::cout << "Output filename?: ";
    std::string filename;
    std::cin >> filename;

    std::ifstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return 1;
    }

    Data out;
    out = out.ReadOutfile(outfile);

    int nat = out.GetNat();
    double forc_conv_thr = 1.0e-4;
    std::vector< std::vector<int> > conv_flag( nat, std::vector<int>(3) );

    for ( int i = 0; i < nat; i++ ) {
        std::vector<double> force = out.GetForce(i);
        for (int j = 0; j < 3; j++) {
            if (fabs(force[j]) > forc_conv_thr && out.GetIfpos(i)[j] == 1) {
                conv_flag[i][j] = 1;
            } else {
                conv_flag[i][j] = 0;
            }
        };
    };

    std::cout << "Convergence not achieved for the following atoms:" << std::endl;

    for ( int i = 0; i < nat; i++ ) {
        if ( conv_flag[i][0] == 1 || conv_flag[i][1] == 1 || conv_flag[i][2] == 1 ) {
            std::cout << "Atom: " << i+1;
            if ( conv_flag[i][0] == 1 ) std::cout << " X";
            if ( conv_flag[i][1] == 1 ) std::cout << " Y";
            if ( conv_flag[i][2] == 1 ) std::cout << " Z";
            std::cout << std::endl;
        }
    }

    return 0;
    
}

