#include <iostream>
#include <iomanip>
#include <time.h>
#include <fstream>
#include <lib.h>
#include <vector>
#include <cmath>
using namespace std;

double pi = 4*atan(1);

double random_number(){
    /****************************************************************************
     *  Random number generator - Bolzmann distribution
     *                Box-Muller transform
     */
    double U1 = rand()/float(RAND_MAX);
    double U2 = rand()/float(RAND_MAX);
    double val = sqrt(-2*log(U1))*cos(2*pi*U2);
    return val;
}

void initialize(double **R, double *V,int N){

    /* Initialize system:
     * - initial particle positions   - fcc
     * - initial particle velocities  - Boltzmann distibution
     */

    double x,y,z,b,N,Nx,Ny,Nz,Lx,Ly,Lz,sigma,T,kB;
    int natom,N_atoms;

    Nx = 3;    // number of origins (from where to place the four atoms) within a box
    Ny = 3;
    Nz = 3;

    b = 1.0;     // 5.26;    // Ångstrøm [Å]

    Lx = Nx*b;   // length of box along x-axis
    Ly = Ny*b;
    Lz = Nz*b;

    N_atoms = Nx*Ny*Nz;  // number of origins within a box
    N = 4*N_atoms;       // number of atoms in one box. 4 atoms per origin.

    kB = 1.480*pow(10,-23); // Bolzmann constant
    T = 100.0;              // Kelvin, initial temperature
    //sigma = sqrt(kB*T/m);   // standarddeviation in temp. from Boltzmann distribution

    // create a vector within a vector, using including the <vector> library.
//    vector < vector < double > > R; // positions
//    vector < vector < double > > V; // velocities

    /*****************************************************************************
     *               Initial positions and velocities
     */

    natom = 0;
    for (int ix=0; ix < Nx; ++ix) {
        for (int iy=0; iy < Ny; ++iy){
            for (int iz=0; iz < Nz; ++iz){
                vector < double > r (3); // temp position vector
                vector < double > v (3); // temp velocity vector

                x = ix*b; y=iy*b;z=iz*b;

                r[0] = x; r[1] = y; r[2] = z;
                R.push_back(r); // [x,y,z]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);

                r[0] = x+0.5; r[1] = y+0.5; r[2] = z;
                R.push_back(r); // [x+0.5,y+0.5,z]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);

                r[0] = x+0.5; r[1] = y; r[2] = z+0.5;
                R.push_back(r); // [x+0.5,y,z+0.5]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);

                r[0] = x; r[1] = y+0.5; r[2] = z+0.5;
                R.push_back(r); // [x,y+0.5,z+0.5]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);


                natom = natom+4;
            }
        }
    }


    /*********************************************************************************************
     *                     Write initial conditions to file
     * string filename = "state000.xyz", open file and write initial postions and velocities;
     */

    double sum_v = 0;
    double quad_sum_v = 0;
    double mean_v;
    double std_dev;

    for (int i=0;i<N;++i){
        cout << "Ar" << " " << R[i][0] << " " << R[i][1] << " " << R[i][2] << " " << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " << endl;
        sum_v = sum_v + (V[i][0] + V[i][1] + V[i][2]);
        quad_sum_v = quad_sum_v + (V[i][0]*V[i][0] + V[i][1]*V[i][1] + V[i][2]*V[i][2]);
    }

    mean_v = sum_v/N;
    std_dev = sqrt(quad_sum_v/N - mean_v*mean_v);
    cout << "----------------------------------------------------------" << endl;
    cout  << "mean velocity = " << mean_v << " with standard deviation = " << std_dev << endl;

    ofstream myfile;
    myfile.open("state000.xyz");
    myfile << N << endl;
    myfile << "initial state fcc lattice of Argon gass" << endl;
    for (int i=0;i<N;++i){
        myfile << "Ar" << " " << R[i][0] << " " << R[i][1] << " " << R[i][2] << " " << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " << endl;
    }
    myfile.close();

}

void integrator(double **R,double **V,int N){

    /* *********************************************************************************************
     * Integrator uses the stable Verlet algorithm to calculate the motion of the particles
     *  using the Lenny-Jones potential to find the force between evry particle.
     */

    double rx,ry,rz,r2,r2i,r6i,r12i,fix,fiy,fiz,fjx,fjy,fjz;

    // f_ij contains the force felt from each particle i on particle j.
    vector < vector < double > > f_ij (N,3);

    int force = 0;
    for (int i = 0; i < N; ++i) {
        rx = R[i][0]; ry = R[i][1]; rz = R[i][2];
        vector < double > f (3);
        for (int j = 0; j < N; ++j) {
            if ( j != i) {
                vector < double > r_ij (3);
                vector < double > fij (3);

                r_ij[0] = rx - R[j][0];
                r_ij[1] = ry - R[j][1];
                r_ij[2] = rz - R[j][2];

                r2 = r_ij[0]*rx_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
                r2i = 1.0/r2;
                r6i = r2i*r2i*r2i;
                r12i = r6*r6;

                //f_ij = 24*(2*r12i - r6i)*r2i*rx*r_ij;
                f_ij[force][0] = 24*(2*r12i - r6i)*r2i*rx_ij;
                f_ij[force][1] = 24*(2*r12i - r6i)*r2i*ry_ij;
                f_ij[force][2] = 24*(2*r12i - r6i)*r2i*rz_ij;

                fij[0] = fix + f_ij[force][0];
                fij[1] = fiy + f_ij[force][1];
                fij[2] = fiz + f_ij[force][2];
                //fjx = fjx - f_ij[force][0];
                //fjy = fjy - f_ij[force][1];
                //fjz = fjz - f_ij[force][2];

                f = f + f.push_back(fij);

                force = force + 1;
            } // end if
        } // end for j

        f_ij.push_back(f);
    } // end for i

    for (int i = 0; i < N; ++i) {
        vector < double > v (3);
        vector < double > r (3);

        v1[0] = V[i][0] + f_ij(i)/2m*dt // dimentions!


    }

}

int main(){

    // create a vector within a vector, using including the <vector> library.
    vector < vector < double > > R; // positions
    vector < vector < double > > V; // velocities
    int N;

    initialize(R,V,N);

    integrator(R,V,R);

    return 0;
}

