#include <iostream>
#include <iomanip>
#include <time.h>
#include <fstream>
#include <lib.h>
#include <vector>
using namespace std;

int main(int argc, char* argv[]){

    double x,y,z,b,bx,by,bz,Nx,Ny,Nz,Lx,Ly,Lz;

    int natom,N_atoms;
    Nx = 3;               // number of atoms along axis
    Ny = 3;
    Nz = 3;
    Lx = 3;               // length of axis
    Ly = 3;
    Lz = 3;
    N_atoms = Nx*Ny*Nz;
    bx = Lx/(double)Nx;
    by = Ly/(double)Ny;
    bz = Lz/(double)Nz;
    cout << bx << " " << by << " " << bz << endl;
    b = 5.26; // Ångstrøm [Å]

    // create a vector within a vector, using including the <vector> library.
    vector < vector < double > > atom (N_atoms*4, vector < double > (3,0));

    natom = 0;
    for (int ix=0; ix < Nx; ++ix) {
        for (int iy=0; iy < Ny; ++iy){
            for (int iz=0; iz < Nz; ++iz){

                x = ix*bx; y=iy*by;z=iz*bz;

                // four atoms has to be placed for every grid point
                atom[natom][0]=x;
                atom[natom][1]=y;
                atom[natom][2]=z;
                atom[natom+1][0]=x+0.5;
                atom[natom+1][1]=y;
                atom[natom+1][2]=z+0.5;
                atom[natom+2][0]=x;
                atom[natom+2][1]=y+0.5;
                atom[natom+2][2]=z+0.5;
                atom[natom+3][0]=x+0.5;
                atom[natom+3][1]=y+0.5;
                atom[natom+3][2]=z;

                natom = natom+4;

            }
        }
    }

    cout << "-------------------------------------------" << endl;
    for (int i = 0; i < N_atoms*4; ++i) {
        cout << " " << atom[i][0] << " " << atom[i][1] << " " << atom[i][2] << " " << endl;
    }
    cout << "-------------------------------------------" << endl;

    //string filename = "state000.xyz", open file and write initial state;
    ofstream myfile;
    myfile.open("state000.xyz");
    myfile << N_atoms*4 << endl;
    myfile << "initial state fcc lattice of Argon gass" << endl;
    for (int i=0;i<N_atoms*4;++i){
        myfile << " " << atom[i][0] << " " << atom[i][1] << " " << atom[i][2] << " " << endl;
    }
    myfile.close();

    /* Now we should give the particles some initial velocity, and then calculate the new positions
     * accordingly. What to do with the boundaries then?
     */

    vector < vector < double > > velocity (N_atoms*4, vector < double > (3,0) );
    double rx_ij,ry_ij,rz_ij, rx,ry,rz, r2,r2i,r6i,r12i;
    vector < vector < double > > f_ij (N_atoms*Natoms*16, vector < double > (3,0));
    int force = 0;
    for (int i = 0; i < N_atoms*4; ++i) {
        rx = atom[i][0]; ry = atom[i][1]; rz = atom[i][2];
        for (int j = 0; j < N_atoms*4; ++j) {
            rx_ij =  rx - atom[j][0];ry_ij =  ry - atom[j][1];ry_ij =  rz - atom[j][2];
            r2 = rx_ij*rx_ij + ry_ij*ry_ij + rz_ij*rz_ij;
            r2i = 1.0/r2;
            r6i = r2i*r2i*r2i;
            r12i = r6*r6;
            f_ij[force][0] = 24*(2*r12i - r6i)*r2i*rx;
            f_ij[force][1] = 24*(2*r12i - r6i)*r2i*ry;
            f_ij[force][2] = 24*(2*r12i - r6i)*r2i*rz;
            fi = fi + f_ij;
            fj = fj - f_ij;

            force = force + 1;
        }
    }

    return 0;
}

