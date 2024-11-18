#ifndef ATOMS_H
#define ATOMS_H

class Atoms {
public:

protected:
    class Atom {
    public:
        double position[3];
        double charge;
        int atomic_number;
    };
    int nat;
    std::vector<Atom> atoms;
};

#endif
