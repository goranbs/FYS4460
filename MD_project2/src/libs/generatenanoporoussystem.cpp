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


    // create shpere holes:
    //spheres(R0,R1,nSpheres); // create spheres
    //create_pores(atoms,nSpheres,N);

    //create cylinder:
    cylinder(atoms,R1);
}

void GenerateNanoPorousSystem::cylinder(vector <Atom> &atoms, double &R){

    // create cylinder that goes throough the system.
    // make all surrounding particles part of the matrix
    // 1) remove particles in the cylinder, or
    // 2) make particles in cylinder part of fluid.


    vector <double> r (3,0);
    double Rx,Ri;
    double random;

    NumberOfFreeParticles = 0;


    auto it = atoms.begin();
    while (it != atoms.end()) {
        r = it->position();

        Rx = 0;
        Ri = 0;
        for (int cor = 0; cor < 2; ++cor) {
            Rx = L[cor]/2.0 - r[cor];

            // minimum image convention:
            if (Rx < -(L[cor]/2.0)) { Rx = Rx + L[cor];}
            else if (Rx > (L[cor]/2.0)) { Rx = Rx - L[cor];}

            Ri += Rx*Rx;
        }
        Ri = sqrt(Ri);

        if (Ri > R) {
            //it = atoms.erase(it);
            it->setIs_matrix(true);
            ++it;
        }
        else{
            //it = atoms.erase(it);
            // we should decrease its density to its half! 50/50 chance of keeping the atom should do the trick.
            random = rand()/float(RAND_MAX); // random number in range [0,1].
            if (random < 0.8) {
                it = atoms.erase(it);

            }
            else{
                it->setIs_matrix(false);
                NumberOfFreeParticles += 1;
                ++it;
            }

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
    vector <int> indexes;
    double Ri, Rx;
    double SphereRad;

    NumberOfFreeParticles = 0;

    for (int n = 0; n<nSpheres; n++){                   // for all generated spheres
        auto it = atoms.begin();
        while (it != atoms.end()) {                     // for all atoms in system

            ri = it->position();                        // position of atom i in vector of atoms
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

            if (Ri <= SphereRad) {                       // if distance from center of sphere, less than sphere radius:
                it = atoms.erase(it);                    // erase atom from vector of atoms.
                // esase atom to just have a look at the matrix.
                //cout << "object erased!!" << endl;   
            }

            else {                                       // if distance from center of sphere is greater than sphere radius:
                it->setIs_matrix(true);
                NumberOfFreeParticles += 1;
                ++it;

            }

            //If we want to use the spheres as the matrix, we need to keep track on which particles that are within the spheres.
            /*if (Ri <= SphereRad) {
                indexes.push_back(*it);
            }
            for (int i = 0; i < indexes.size(); ++i) {
                atoms[indexes(i)].setIs_matrix(true);   // the atoms which are within the spheres is part of the matrix.
            }*/
        }
    }
}

int GenerateNanoPorousSystem::numberOfFreeParticles(){
    return NumberOfFreeParticles;
}
double GenerateNanoPorousSystem::density(){
    // density of the fluid
    return dens;
}
double GenerateNanoPorousSystem::volume(){
    // volume of the fluid
    return vol;
}
