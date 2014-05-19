#ifndef GENERATENANOPOROUSSYSTEM_H
#define GENERATENANOPOROUSSYSTEM_H

#include <atom.h>
#include <stdlib.h>  // rand()
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class GenerateNanoPorousSystem
{
public:
    GenerateNanoPorousSystem(vector <Atom >&atoms, double &R0, double &R1, double &Lx, double &Ly, double &Lz, int &nSpheres, int &N);

private:
    void spheres(double &R0, double &R1, int &nSpheres);
    void create_pores(vector <Atom> &atoms, int &nSpheres, int &N);
    vector < vector < double > > spherePos;
    vector <double> sphereRad;
    vector <double> L;

};

#endif // GENERATENANOPOROUSSYSTEM_H
