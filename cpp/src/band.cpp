#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include "band.h"
#include "data.h"

void trim_dash(std::string& str) {
    auto end = std::remove(str.begin(), str.end(), '\'');
    str.erase(end, str.end());
}

namespace band {
    int run() {
        std::cout << "Running band module..." << std::endl;

        // Read the QE input and output file name
        std::cout << "Enter the pw.x input / output: ";
        std::string input_file_name, output_file_name;
        std::cin >> input_file_name >> output_file_name;

        // Open the input and output files
        std::ifstream input_file(input_file_name);
        std::ifstream output_file(output_file_name);
        if ( !input_file || !output_file ) {
            std::cerr << "Something wrong. Check the file you entered." << std::endl;
            return 1;
        }

        // Create Data instance and read the data
        Data pw_in;
        Data pw_out;

        pw_in = pw_in.ReadInfile( input_file );
        pw_out = pw_out.ReadOutfile( output_file );

        Data band_in = pw_in;

        band_in.SetParam("control", "calculation", "'bands'");
        
        // Set the kpoints
        // ------------------------------------------------------
        // This feature is not completely implemented yet
        // We will make it possible to set the kpoints interactively by the user
        // Here, the surface is supposed to be hexagonal,
        // and the kpoints path are set to G -> K -> M -> G
        std::string kpoints_style = "tpiba_b";
        int nks = 4;
        std::vector< std::vector<float> > xk;
        std::vector<float> wk;
        xk.resize(nks);
        wk.resize(nks);

        xk[0] = { 0.0f, 0.0f, 0.0f };
        xk[1] = { 1.0f/2.0f, -1.0f/2.0f/std::sqrt(3.0f), 0.0f };
        xk[2] = { 2.0f/3.0f, 0.0f, 0.0f };
        xk[3] = { 0.0f, 0.0f, 0.0f };

        wk[0] = 10.0f;
        wk[1] = 10.0f;
        wk[2] = 10.0f;
        wk[3] = 10.0f;

        band_in.SetKpoints4bands(kpoints_style, nks, xk, wk);
        // ------------------------------------------------------
        
        // Make the directory for the band calculation
        std::string band_dir = "band";

        if ( !std::filesystem::exists(band_dir) ) {
            std::filesystem::create_directory(band_dir);
        } else {
            std::cout << "Caution: The directory named 'band' has already existed." << std::endl;
            std::cout << "Do you want to overwrite it? (y/n): ";
            char answer;
            std::cin >> answer;
            if ( answer != 'y' && answer != 'Y' ) {
                std::cerr << "Exiting..." << std::endl;
                return 1;
            }
        }

        // Copy the outdir from current directory to band_dir
        std::string outdir = pw_in.GetParam("control", "outdir");
        std::cout << outdir << std::endl;
        trim_dash(outdir);
        
        // Check if the outdir exists in ./
        if ( !std::filesystem::exists(outdir) ) {
            std::cerr << "Error: The directory named '" << outdir << "' does not exist." << std::endl;
            return 1;
        }

        // Check if the outdir exists in ./band/
        // If it exists, remove it and copy the new one
        if ( std::filesystem::exists(band_dir + "/" + outdir) ) {
            std::filesystem::remove_all(band_dir + "/" + outdir);
        }
        std::filesystem::copy(outdir, band_dir + "/" + outdir, std::filesystem::copy_options::recursive);

        // Write the input file for the band calculation
        std::string band_input_file_name = band_dir + "/in";

        std::ofstream band_input_file(band_input_file_name);

        band_in.WriteBandInfile(band_input_file);

        return 0;
    }
}
