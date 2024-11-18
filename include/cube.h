#ifndef CUBE_H
#define CUBE_H

#include <string>
#include <vector>
#include <array>
#include "atoms.h"

class Cube : protected Atoms {
public:
    void LoadCubeFile(std::string filename);
    void WriteCubeFile(std::string filename);
    void ShiftVolumaticData();
    void ShiftAtomicPosition();

    std::array<int, 3> Getgrid();
    std::array<int, 3> Getorigin();
    std::array<std::array<double, 3>, 3> Getstepw();
    std::array<double, 3> Getalat();
    std::vector<double> Getvol_1d();
protected:
    std::array<int, 3> grid;
    std::array<int, 3> origin;
    std::array<std::array<double, 3>, 3> stepw;
    std::array<double, 3> alat;
    std::vector<double> vol_1d;
};

#endif
