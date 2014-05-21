#include "generatenanoporoussystem.h"

GenerateNanoPorousSystem::GenerateNanoPorousSystem(vector < Atom > &atoms, double &R0, double &R1, double &Lx, double &Ly, double &Lz, int &nSpheres, int &N){
    /* set up the nanoporous system, given
     *   -the atoms in the system
     *   -the range of uniformly distributed random radiuses. The radius defines a sphere that acts as a pore within the matrix.
     *   -the system volume V = Lx*Ly*Lz
     *   -the spheres should be centered at random positions within the volume V.
     *   -nSpheres is the number of spheres in the system
     *   -N is the number of atoms in atoms.
     */
    //nSpheres = 1;
    spherePos = vector < vector <double> > (nSpheres, vector < double > (3,0));
    sphereRad = vector <double> (nSpheres,0);
    L = vector <double> (3,0);
    L[0] = Lx;
    L[1] = Ly;
    L[2] = Lz;

    // create spheres:
    spheres(R0,R1,nSpheres); // create spheres
    create_pores(atoms,nSpheres,N);

    //create cylinder:
    //cylinder(atoms,R1);
}

void GenerateNanoPorousSystem::cylinder(vector < Atom > &atoms, double &R){
    // create cylinder that goes throough the system.
    // make all surrounding particles part of the matrix
    // 1) remove particles in the cylinder, or
    // 2) make particles in cylinder part of fluid.

    vector < double > r (3,0);
    vector < double > Rc (2,0);
    double Rx,Ri;
    Rc = {L[0]/2.0, L[1]/2.0};

    auto it = atoms.begin();
    while (it != atoms.end()) {
        r = it->position();
        Ri = 0.0;
        Rx = 0.0;
        for (int cor = 0; cor < 2; ++cor) {
            Rx = Rc[cor] - r[cor];
            if (Rx < -L[cor]/2.0) {
                Rx = Rx + L[cor];
            }
            else if (Rx > L[cor]/2.0) {
                Rx = Rx - L[cor];
            }
            Ri = Ri + Rx*Rx;
        }
        Ri = sqrt(Ri);
        if (Ri > R) {                            // all atoms outside the cylinder
            //it = atoms.erase(it);
            //it->setIs_matrix(true);
            ++it;
        }
        else{  // all atoms within cylinder
            it = atoms.erase(it);
            //it->setIs_matrix(true);
            //++it;
        }
    }
}

void GenerateNanoPorousSystem::spheres(double &R0, double &R1, int &nSpheres){
    double range = R1 - R0;
    for (int n = 0; n < nSpheres; ++n) {
        spherePos[n][0] =  (rand()/float(RAND_MAX))*L[0];
        spherePos[n][1] =  (rand()/float(RAND_MAX))*L[1];
        spherePos[n][2] =  (rand()/float(RAND_MAX))*L[2];
        sphereRad[n] = R0 + (range*rand()/float(RAND_MAX));

    }
//    spherePos[0][0] =  0;//L[0]/2.0;
//    spherePos[0][1] =  0;//L[1]/2.0;
//    spherePos[0][2] =  0;//L[2]/2.0;
//    sphereRad[0] = R0 + (range*rand()/float(RAND_MAX));
}

void GenerateNanoPorousSystem::create_pores(vector <Atom > &atoms, int &nSpheres, int &N){
    vector < double > ri (3,0);

    double Ri, Rx;
    double SphereRad;

    for (int n = 0; n<nSpheres; n++){                   // for all generated spheres
        auto it = atoms.begin();
        while (it != atoms.end()) {                     // for all atoms in system
            ri = it->position();                        // position of atom i in vector of atoms.
            //it->setIs_matrix(true);                   // all atoms part of matrix
            Ri = 0;

            for (int cor = 0; cor < 3; ++cor) {
                Rx = spherePos[n][cor] - ri[cor];       // atoms distance from sphere center.
                // periodic boundary conditions:
                if (Rx < -L[cor]/2.0) {
                    Rx = Rx + L[cor];
                }
                else if (Rx > L[cor]/2.0) {
                    Rx = Rx - L[cor];
                }
                Ri = Ri + Rx*Rx;
            }

            Ri = sqrt(Ri);
            SphereRad = sphereRad[n];
            if (Ri <= sphereRad[n]) {                      // if distance from center of sphere, less than sphere radius:
                //it = atoms.erase(it);                    // erase atom from vector of atoms.
                // esase atom to just have a look at the matrix.
                //cout << "object erased!!" << endl;
                //it->setIs_matrix(false);                   // atom part of fluid
                it->setIs_matrix(true);                      // atom part of matrix
                ++it;

            }
            else{  // all atoms outside the spheres
                it = atoms.erase(it);
                //++it;
            }
        }
    }
}

