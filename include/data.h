#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <any>
#include <unordered_map>

#pragma once

// This class is used to store the data from the output file of the QuantumESPRESSO.
//
// class Data has the following attributes:
// - nat:   number of atoms in the cell
// - ntyp:  number of atom types
// - atoms: a vector of class Atom
//
// class Atom has the following attributes:
// - symbol:    the symbol of the atom
// - pos:       3D vector defining the position of the atom
// - ifpos:     3D vector defining whether the each component of the position is fixed or not
// - force:     3D vector defining the force acting on the atom
//

class Data {
public:

    // Following three functions are mainly used in the ReadOutfile
    void SetNat(int n); 
    void SetNtyp(int n);
    void SetAtom( std::string symbol, std::vector<double> pos, std::vector<int> ifpos, std::vector<double> force);
    void SetCalculation(std::string calc);
    void SetRestartmode(std::string mode);
    void SetParam( std::unordered_map<std::string, std::any> param );

    // Use these functions to get the data
    int GetNat(); 
    std::vector<double> GetForce(int i);
    std::vector<int> GetIfpos(int i);
    std::string GetCalculation();
    std::string GetRestartmode();

    template <typename T>
    T GetParam( std::string key );

    // Read the io file of QuantumESPRESSO
    Data ReadInfile( std::ifstream &file );
    Data ReadOutfile( std::ifstream &file );

private:

    class Atom {
    public:
        std::string symbol;
        std::vector<double> pos;
        std::vector<int> ifpos;
        std::vector<double> force;
    };

    std::vector<Data::Atom> atoms;

    int nat;
    int ntyp;

    std::string calculation;
    std::string restartmode;

    std::unordered_map<std::string, std::any> params;

};
