#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <any>
#include <unordered_map>
#include <variant>
#include <array>

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
    // using ParamType = std::variant< int, double, std::string >;

    // Following three functions are mainly used in the ReadOutfile
    void SetNat(int n); 
    void SetNtyp(int n);
    void SetAtom( std::string symbol, std::vector<double> pos, std::vector<int> ifpos, std::vector<double> force = {});
    void SetCalculation(std::string calc);
    void SetRestartmode(std::string mode);
    void SetParam( std::unordered_map<std::string, std::unordered_map<std::string, std::string>> param );
    void SetPseudo( std::string symbol, float mass, std::string name );
    void SetKpoints( const std::string& style, const std::array<int, 3>& nk, const std::array<int, 3>& sk );
    void SetCell( const std::string& style, const std::vector< std::vector<double> >& v );

    void SetParam( const std::string& section, const std::string& key, const std::string& value );

    void SetKpoints4bands(
            const std::string& style,
            const int& nks,
            const std::vector< std::vector<float> >& xk,
            const std::vector<float>& wk
            );

    // Use these functions to get the data
    int GetNat(); 
    std::vector<double> GetForce(int i);
    std::vector<int> GetIfpos(int i);
    std::string GetCalculation();
    std::string GetRestartmode();

    std::string GetParam( const std::string& section, const std::string& key );

    // Read the io file of QuantumESPRESSO
    Data ReadInfile( std::ifstream &file );
    Data ReadOutfile( std::ifstream &file );

    void WriteBandInfile( std::ofstream &file );

private:

    class Atom {
    public:
        std::string symbol;
        std::vector<double> pos;
        std::vector<int> ifpos;
        std::vector<double> force;
    };

    std::vector<Atom> atoms;

    int nat;
    int ntyp;

    std::string calculation;
    std::string restartmode;

    class Pseudo {
    public:
        std::string symbol;
        float mass;
        std::string name;
        
        Pseudo( std::string s, float m, std::string n ):
            symbol( std::move(s) ), mass( std::move(m) ), name( std::move(n) ) {}
    };

    std::vector<Pseudo> pseudos_;

    class Kpoints {
    public:
        std::string style;
        std::array<int, 3> nk;
        std::array<int, 3> sk;

        int nks;
        std::vector< std::vector<float> > xk;
        std::vector<float> wk;
    };

    Kpoints kpoints;

    class Cell {
    public:
        std::string style;
        std::vector< std::vector<double> > v;
    };

    Cell cell;
    
    std::unordered_map<std::string, std::unordered_map< std::string, std::string > > params;

};
