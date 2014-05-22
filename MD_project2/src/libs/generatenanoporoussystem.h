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
    int numberOfFreeParticles();
    double density();
    double volume();

private:
    void cylinder(vector <Atom> &atoms, double &R);
    void spheres(double &R0, double &R1, int &nSpheres);
    void create_pores(vector <Atom> &atoms, int &nSpheres);
    vector < vector < double > > spherePos;
    vector <double> sphereRad;
    vector <double> L;
    int NumberOfFreeParticles;
    int counter;
    int Nparticles;
    double dens;
    double vol;


};

#endif // GENERATENANOPOROUSSYSTEM_H
