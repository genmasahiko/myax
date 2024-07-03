#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>



#include "data.h"


int main() {

    std::string infilename;
    std::string outfilename;
    std::cout << "Enter the input and output file names:" << std::endl;

    std::cin >> infilename >> outfilename;

    std::ifstream infile(infilename);
    std::ifstream outfile(outfilename);

    Data in;
    Data out;
    in = in.ReadInfile(infile);
    out = out.ReadOutfile(outfile);

    // Change filenames
    std::string old_infilename = infilename + ".0001";
    std::string old_outfilename = outfilename + ".0001";
    int count = 1;
    while ( system( ("ls " + old_infilename + " > /dev/null 2>&1" ).c_str() ) == 0 ) {
        count++;
        std::ostringstream oss;
        oss << "." << std::setw(4) << std::setfill('0') << count;
        old_infilename = infilename + oss.str();
        old_outfilename = outfilename + oss.str();
    }
    system( ("mv " + infilename + " " + old_infilename).c_str() );
    system( ("mv " + outfilename + " " + old_outfilename).c_str() );

    std::cout << in.GetParam("calculation") << std::endl;

    // if ( in.GetParam("calculation") == "neb" ) {
    //     std::cout << "The calculation is NEB." << std::endl;

    //     std::ofstream new_infile(infilename);
    //     std::string line;
    //     
    //     while ( std::getline(infile, line) ) {
    //         new_infile << line << std::endl;

    //         if ( line.find("&PATH") != std::string::npos ||
    //              line.find("&path") != std::string::npos ) {
    //             if ( in.GetRestartmode() != "restart" ) {
    //                 new_infile << "  restart_mode = 'restart'" << std::endl;
    //             }
    //         }
    //     }
    // }

    return 0;

}
