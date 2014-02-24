#include "initialstate.h"
#include <vector>
#include <random>
#include <cmath>

using namespace std;

/* Class that generates the initial system state of MD-system
 * Initial positions and velocities
 */

InitialState::InitialState(){
}
void InitialState::two_particles(vector < vector < double > > &R, vector < vector < double > > &V, double Lx, double Ly, double Lz){
    R[0][0] = Lx/3.0;
    R[0][1] = 0;
    R[0][2] = 0;
    R[1][0] = 2*Lx/3.0;
    R[1][1] = 0;
    R[1][2] = 0;


    V[0][0] = -1.0;
    V[0][1] = 0;
    V[0][2] = 0;
    V[1][0] = 1.0;
    V[1][1] = 0;
    V[1][2] = 0;


}

double InitialState::boltzmann_randomnr(){
    /****************************************************************************
     *  Random number generator - Bolzmann distribution
     *                Box-Muller transform
     */
    double pi = 4*atan(1);
    double U1 = rand()/float(RAND_MAX);
    double U2 = rand()/float(RAND_MAX);
    double val = sqrt(-2*log(U1))*cos(2*pi*U2);
    return val;
}
