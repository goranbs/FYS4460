#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include <vector>

using namespace std;

class InitialState
{
public:
    InitialState();
    void two_particles(vector < vector < double > > &R, vector < vector < double > > &V, double Lx, double Ly, double Lz);

private:
    double boltzmann_randomnr();

};
#endif // INITIALSTATE_H
