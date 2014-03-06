#ifndef ATOM_H
#define ATOM_H

#include <vector>

using namespace std;

class Atom
{
public:
    Atom();
    const int AtomIndex();
private:
    vector < double > position();
    vector < double > velocity();
    vector < double > travel(int t);
    vector < double > intitial_position();

};

#endif // ATOM_H
