#include <iostream>
#include <iomanip>
#include <time.h>
#include <fstream>
#include <lib.h>
#include <vector>
#include <cmath>
using namespace std;

double pi = 4*atan(1);

void load_initial_state(){

    double N;
    ifstream ifile;
    ifile.open('state000.xyz');
    if(!file){
        cout << "Error loading the file" << endl;
    }

    N >> ifile.read();
    print N
    vector < vector < double > > R (N,3);
    vector < vector < double > > V (N,3);








/* *********************************************************************************************
 * Now we should give the particles some initial velocity.
 * -First we have to calculate the
 */

//    // r_ij is holding vectors of distances between every particle in the grid:
//    vector < vector < double > > r_ij;
//    // velocity should contain the velocity of every particle in the grid:
//    vector < vector < double > > velocity;

//    double rx_ij,ry_ij,rz_ij, rx,ry,rz, r2,r2i,r6i,r12i,fix,fiy,fiz,fjx,fjy,fjz;
//    // f_ij contains the force felt from each particle i on particle j.
//    vector < vector < double > > f_ij (N_atoms*N_atoms*16, vector < double > (3,0));
//    int force = 0;
//    for (int i = 0; i < N_atoms*4; ++i) {
//        rx = atom[i][0]; ry = atom[i][1]; rz = atom[i][2];
//        for (int j = 0; j < N_atoms*4; ++j) {
//            r_ij[force][0] = rx - atom[j][0];
//            r_ij[force][1] = ry - atom[j][1];
//            r_ij[force][2] = rz - atom[j][2];

//            rx_ij =  rx - atom[j][0];ry_ij =  ry - atom[j][1];ry_ij =  rz - atom[j][2];
//            r2 =rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij;
//            r2i = 1.0/r2;
//            r6i = r2i*r2i*r2i;
//            r12i = r6*r6;

//            //f_ij = 24*(2*r12i - r6i)*r2i*rx*r_ij;
//            f_ij[force][0] = 24*(2*r12i - r6i)*r2i*rx_ij;
//            f_ij[force][1] = 24*(2*r12i - r6i)*r2i*ry_ij;
//            f_ij[force][2] = 24*(2*r12i - r6i)*r2i*rz_ij;

//            fix = fix + f_ij[force][0];
//            fiy = fiy + f_ij[force][1];
//            fiz = fiz + f_ij[force][2];
//            fjx = fjx - f_ij[force][0];
//            fjy = fjy - f_ij[force][1];
//            fjz = fjz - f_ij[force][2];

//            force = force + 1;
//        }
//    }

}

int main(){
    load_initial_state();

return 0;
}
