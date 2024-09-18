#include "cube.h"

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>

double VolumeElement ( Cube &cube ) {
    double ve = 0.0;
    std::array<std::array<double, 3>, 3> stepw = cube.Getstepw();
    
    // Calculate the scalar triple product of the stepw vectors
    ve = stepw[0][0] * ( stepw[1][1] * stepw[2][2] - stepw[1][2] * stepw[2][1] )
        - stepw[0][1] * ( stepw[1][0] * stepw[2][2] - stepw[1][2] * stepw[2][0] )
        + stepw[0][2] * ( stepw[1][0] * stepw[2][1] - stepw[1][1] * stepw[2][0] );

    return ve;
}

double Integral ( Cube &cube ) {
    std::array<int, 3> grid = cube.Getgrid();
    std::array<double, 3> alat = cube.Getalat();
    std::vector<double> vol_1d = cube.Getvol_1d();

    double integrate = 0.0;
    double ve = VolumeElement( cube );
    int numgrid = grid[0] * grid[1] * grid[2];

    for ( int i = 0; i < numgrid; i++ ) {
        integrate += std::abs(vol_1d[i]);
    }

    integrate = integrate * ve ;//* alat[0] * alat[1] * alat[2] ;/// numgrid;

    return integrate;
}

int main(int argc, char *argv[]) {
    Cube cube;
    double integrate;

    std::string filename;

    if ( argc == 2 ) {
        filename = argv[1];
    } else {
        std::cout << "Enter the File name to be integrated: ";
        std::cin >> filename;
    }

    cube.LoadCubeFile( filename );
    integrate = Integral( cube );

    std::cout << "The integration of the cube file is: " << integrate << std::endl;

    return 0;
}
