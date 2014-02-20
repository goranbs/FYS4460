#include <iostream>
#include <iomanip>
#include <time.h>
#include <fstream>
//#include <lib.h>
#include <vector>
#include <cmath>
#include <random>
#include <unitconverter.h>
#include <list>

using namespace std;

double pi = 4*atan(1);

// Check out the constants, do they fit with those in the project text?
//double b = 5.260;                 // Ångstrøm [Å]
double b = 20.0;                 // Ångsrøm [Å]
double mA = 39.948;               // mass of Argon [amu]
double kB = 1.480*pow(10,-23);    // Bolzmann constant [eV/K]
double eps = 0.1*1.0303;          // Energy constant [eV]
double Temp = 100.0;              // Kelvin, initial temperature
double sigma = 3.405;             // Ångstrøm, scalefactor - Leonard-Jones
/********************************************************************************
 * Conversion factors, so that we get out units that we would like to use
 */

UnitConverter leng;
double length = leng.from_aangstrom(b);  // unitless length

UnitConverter nrj;
double E = nrj.from_energy(eps);         // unitless energy

UnitConverter mass;
double m = mass.from_amu(mA);            // unitless mass

double Time_0 = sigma*sqrt(mA/eps);      // Conversion factor time
UnitConverter time_0;
double Time = time_0.from_time(Time_0);  // unitless time

double F_0 = mA*sigma/(Time_0*Time_0);   // Conversion factor force

double T_0 = 119.74; // [K] eps/kB;      // Conversion factor temperature

double velocity = sigma/Time_0;          // Conversion factor velocity


/********************************************************************************
 *                    Random Number Generator
 * ******************************************************************************
 */
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


/********************************************************************************
 *                    Initialize Box List
 * ******************************************************************************
 */
void initialize_box_list(double Lcx, double Lcy, double Lcz , int nx, int ny, int nz, vector < vector < double > > &R, vector < list < int > > &box_list){
    /* function that creates the box_list
     */

    int ix,iy,iz,box_index;
    for (int p = 0; p < R.size(); ++p) {
        // Create x,y and z index for particle p.
        ix = floor(R[p][0]/Lcx);
        iy = floor(R[p][1]/Lcy);
        iz = floor(R[p][2]/Lcz);
        box_index = ix*ny*nz + iy*nz + iz;     // cubic to linear transform, (nested lists)
        box_list[box_index].push_back(p);      // put particle p into box number box_index
    }

}

/********************************************************************************
 *                    Update Box List
 * ******************************************************************************
 */
void update_box_list(double Lcx, double Lcy, double Lcz, int nx, int ny, int nz, vector < vector < double > > &R,vector < list < int > > &box_list){
    /* Function that updates the box list :-) - that is, it puts the particle into the box it belongs to:-)
     */

    for(int i = 0; i < box_list.size(); ++i){
        box_list[i].clear(); // empty the box list!
    }

    int ix,iy,iz,box_index;
    for (int p = 0; p < R.size(); ++p) {
        // Create x,y and z index for particle p.
        ix = floor(R[p][0]/Lcx);
        iy = floor(R[p][1]/Lcy);
        iz = floor(R[p][2]/Lcz);
        box_index = ix*ny*nz + iy*nz + iz;   // cubic to linear transform, (nested lists)
        box_list[box_index].push_back(p);        // put particle p into box number box_index
    }


}

/********************************************************************************
 *                    Initialize state
 * ******************************************************************************
 */
void initialize(vector < vector < double > > &V, vector < vector < double > > &R, int &N, int Nx, int Ny, int Nz){

    /*   Initialize system:
     * - initial particle positions   - fcc
     * - initial particle velocities  - Boltzmann distibution
     */

    double x,y,z;
    int natom,N_atoms;

    double Lx,Ly,Lz;
    Lx = Nx*length;
    Ly = Ny*length;
    Lz = Nz*length;
    cout << "Nx Ny Nz" << endl;
    cout << Nx << " " << Ny << " " << Nz << endl;

    N_atoms = Nx*Ny*Nz;  // number of origins within a box
    N = 4*N_atoms;       // number of atoms in one box. 4 atoms per origin.

    /*****************************************************************************
     *               Initial positions and velocities
     */


    natom = 0;
    for (int ix=0; ix < Nx; ++ix) {
        for (int iy=0; iy < Ny; ++iy){
            for (int iz=0; iz < Nz; ++iz){
                vector < double > r (3); // temp position vector
                vector < double > v (3); // temp velocity vector

                x = ix*length; y=iy*length;z=iz*length;

                r[0] = x; r[1] = y; r[2] = z;
                R.push_back(r); // [x,y,z]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);

                r[0] = x+0.5*length; r[1] = y+0.5*length; r[2] = z;
                R.push_back(r); // [x+0.5,y+0.5,z]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);

                r[0] = x+0.5*length; r[1] = y; r[2] = z+0.5*length;
                R.push_back(r); // [x+0.5,y,z+0.5]
                v[0] = random_number();v[1] = random_number();v[2] = random_number();
                V.push_back(v);

                r[0] = x; r[1] = y+0.5*length; r[2] = z+0.5*length;
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

    double sum_v_x = 0;
    double sum_v_y = 0;
    double sum_v_z = 0;
    double mean_x,mean_y,mean_z;
    double quad_sum_x = 0;
    double quad_sum_y = 0;
    double quad_sum_z = 0;
    double std_dev_x;
    double std_dev_y;
    double std_dev_z;

    for (int i=0;i<N;++i){
        //cout << "Ar" << " " << R[i][0] << " " << R[i][1] << " " << R[i][2] << " " << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " << endl;
        sum_v_x = sum_v_x + V[i][0];
        sum_v_y = sum_v_y + V[i][1];
        sum_v_z = sum_v_z + V[i][2];
        quad_sum_x = quad_sum_x + V[i][0]*V[i][0];
        quad_sum_y = quad_sum_y + V[i][1]*V[i][1];
        quad_sum_z= quad_sum_z + V[i][2]*V[i][2];
    }

    mean_x = sum_v_x/float(N);
    mean_y = sum_v_y/float(N);
    mean_z = sum_v_z/float(N);
    quad_sum_x = quad_sum_x/float(N);
    quad_sum_y = quad_sum_y/float(N);
    quad_sum_z = quad_sum_z/float(N);

    std_dev_x = sqrt(quad_sum_x - mean_x*mean_x);
    std_dev_y = sqrt(quad_sum_y - mean_y*mean_y);
    std_dev_z = sqrt(quad_sum_z - mean_z*mean_z);

    cout << "----------------------------------------------------------" << endl;
    cout << " mean speed x: " << mean_x << " stdev= " << std_dev_x << endl;
    cout << " mean speed y: " << mean_y << " stdev= " << std_dev_y << endl;
    cout << " mean speed z: " << mean_z << " stdev= " << std_dev_z << endl;
    cout  << "---------------------------------------------------------" << endl;

//    ofstream myfile;
//    myfile.open("state000.txt");
//    myfile << N << endl;
//    myfile << "initial_state_fcc_lattice_of_Argon_gass" << " " << 0.0 << endl;
//    for (int i=0;i<N;++i){
//        myfile << "Ar" << " " << R[i][0] << " " << R[i][1] << " " << R[i][2] << " " << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << endl;
//    }
//    myfile.close();

}

void Lennard_Jones(vector < vector < double > > &F, vector < vector < double > > &R, vector < double > &U, int &N, double Lx, double Ly, double Lz, int N_cells_x, int N_cells_y, int N_cells_z, vector < list < int > > &box_list){
    /* The Lenny-Jones potential updates the forces F on particle i in position R.
     */


    double rx,ry,rz,r2,r2i,r6i,r12i;
    int force = 0;
    for (int box_x = 0; box_x < N_cells_x; ++box_x){
        for (int box_y = 0; box_y < N_cells_y; ++box_y){
            for (int box_z = 0; box_z < N_cells_z; ++box_z){

                int box_index = box_x*N_cells_x*N_cells_y + box_y*N_cells_z + box_z;

                for (auto it = box_list[box_index].begin(); it != box_list[box_index].end(); ++it){
                    int ai = *it;  // ai = atom_index in box_list

                    vector < double > f (3) ; // sums up the total force on particle ai.
                    rx = R[ai][0]; ry = R[ai][1]; rz = R[ai][2]; // x,y,z position of atom ai in box box_index.
                    double Ui = 0; // sums up the potential energy for particle ai.

                    for (int ix=-1; ix<=1; ++ix){ // box number -1,0,1
                        for (int iy=-1; iy<=1; ++iy){
                            for (int iz=-1; iz<=1; ++iz){
                                box_x = box_x + ix;
                                box_y = box_y + iy;
                                box_z = box_z + iz;


                                if (box_x < 0){ box_x = box_x + N_cells_x;}
                                else if (box_x >= N_cells_x){ box_x = box_x - N_cells_x;}

                                if (box_y < 0){ box_y = box_y + N_cells_y; }
                                else if (box_y >= N_cells_y){box_y = box_y - N_cells_y;}

                                if (box_z < 0){ box_z = box_z + N_cells_z;}
                                else if (box_z >= N_cells_z){box_z = box_z - N_cells_z;}

                                // neighbour_index cubic to linear transform:
                                int neighbour = box_x*N_cells_y*N_cells_z + box_y*N_cells_z + box_z;


                                if (neighbour == box_index){
                                    auto iterator = it;
                                    for (auto it2 = ++ iterator; it2 != box_list[neighbour].end(); ++it2){
                                        int aj = *it2;

                                        vector < double > fij (3);
                                        vector < double > r_ij (3);

                                        //cout << "Hi there, Im in auto it2" << endl;
                                        r_ij[0]=  rx - R[aj][0];
                                        r_ij[1] = ry - R[aj][1];
                                        r_ij[2] = rz - R[aj][2];

                                        // Periodic boundary conditions
                                        if (r_ij[0] > Lx/2){r_ij[0] = - Lx + r_ij[0];}
                                        else if (r_ij[0] < -Lx/2){r_ij[0] = Lx + r_ij[0];}

                                        if (r_ij[1] > Ly/2){r_ij[1] = -Ly + r_ij[1];}
                                        else if (r_ij[1] < -Ly/2){r_ij[1] = Ly + r_ij[1];}

                                        if (r_ij[2] > Lz/2){r_ij[2] = -Lz + r_ij[2];}
                                        else if (r_ij[2] < -Lz/2){r_ij[2] = Lz + r_ij[2];}

                                        r2 = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
                                        r2i = 1.0/r2;
                                        r6i = r2i*r2i*r2i;
                                        r12i = r6i*r6i;

                                        fij[0] = 24*(2*r12i - r6i)*r2i*r_ij[0];  // force from j on i.
                                        fij[1] = 24*(2*r12i - r6i)*r2i*r_ij[1];
                                        fij[2] = 24*(2*r12i - r6i)*r2i*r_ij[2];

                                        f[0] = f[0] + fij[0]; // adding up the forces on particle i in x direction.
                                        f[1] = f[1] + fij[1];
                                        f[2] = f[2] + fij[2];


                                        Ui = Ui + 4*(r12i - r6i); // sum up the potential energy for particle i.

                                        force = force + 1;

                                    }
                                }
                                else { // neighbour != box_index, that is, we are looking at the surrounding boxes!
                                    for (auto it3 = box_list[neighbour].begin(); it3 != box_list[neighbour].end(); ++it3){
                                        int aj = *it3;

                                        //cout << "Hi there, im in auto it3" << endl;

                                        vector < double > fij (3);
                                        vector < double > r_ij {
                                            rx - R[aj][0],
                                            ry - R[aj][1],
                                            rz - R[aj][2]
                                        };

                                        //cout << "Hi there, Im in auto it2" << endl;
                                        r_ij[0]= rx - R[aj][0];
                                        r_ij[1] = ry - R[aj][1];
                                        r_ij[2] = rz - R[aj][2];

                                        // Periodic boundary conditions
                                        if (r_ij[0] > Lx/2){r_ij[0] = - Lx + r_ij[0];}
                                        else if (r_ij[0] < -Lx/2){r_ij[0] = Lx + r_ij[0];}

                                        if (r_ij[1] > Ly/2){r_ij[1] = -Ly + r_ij[1];}
                                        else if (r_ij[1] < -Ly/2){r_ij[1] = Ly + r_ij[1];}

                                        if (r_ij[2] > Lz/2){r_ij[2] = -Lz + r_ij[2];}
                                        else if (r_ij[2] < -Lz/2){r_ij[2] = Lz + r_ij[2];}

                                        r2 = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
                                        r2i = 1.0/r2;
                                        r6i = r2i*r2i*r2i;
                                        r12i = r6i*r6i;

                                        fij[0] = 24*(2*r12i - r6i)*r2i*r_ij[0];  // force from j on i.
                                        fij[1] = 24*(2*r12i - r6i)*r2i*r_ij[1];
                                        fij[2] = 24*(2*r12i - r6i)*r2i*r_ij[2];

                                        f[0] = f[0] + fij[0]; // adding up the forces on particle i in x direction.
                                        f[1] = f[1] + fij[1];
                                        f[2] = f[2] + fij[2];


                                        Ui = Ui + 4*(r12i - r6i); // sum up the potential energy for particle i.

                                        force = force + 1;

                                    }
                                }// end if/else.
                            }
                        }
                    }
                    F[ai] = f;   // total force in [x,y,z] on particle i
                    U.push_back(Ui); // total potential energy
                }
            }
        }
    }



//    double rx,ry,rz,r2,r2i,r6i,r12i;
//    //vector < vector < vector < double > > > mirror (N, vector < vector < double > > (26, vector < double > (3,0)));

//    int force = 0;
//    for (int i = 0; i < N; ++i) {
//        rx = R[i][0]; ry = R[i][1]; rz = R[i][2];
//        vector < double > f (3);   //  sums up to total force on i
//        double Ui = 0; // potential energy
//        for (int j = 0; j < N; ++j) {
//            if ( j != i) {
//                vector < double > r_ij (3);
//                vector < double > fij (3);

//                r_ij[0] = rx - R[j][0];  // x distance between particle i and j.
//                r_ij[1] = ry - R[j][1];  // y
//                r_ij[2] = rz - R[j][2];  // z


//                // Periodic boundary conditions
//                if (r_ij[0] > Lx/2){
//                    r_ij[0] = - Lx + r_ij[0];
//                }
//                else if (r_ij[0] < -Lx/2){
//                    r_ij[0] = Lx + r_ij[0];
//                }
//                if (r_ij[1] > Ly/2){
//                    r_ij[1] = -Ly + r_ij[1];
//                }
//                else if (r_ij[1] < -Ly/2){
//                    r_ij[1] = Ly + r_ij[1];
//                }
//                if (r_ij[2] > Lz/2){
//                    r_ij[2] = -Lz + r_ij[2];
//                }
//                else if (r_ij[2] < -Lz/2){
//                    r_ij[2] = Lz + r_ij[2];
//                }

//                r2 = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
//                r2i = 1.0/r2;
//                r6i = r2i*r2i*r2i;
//                r12i = r6i*r6i;

//                fij[0] = 24*(2*r12i - r6i)*r2i*r_ij[0];  // force from j on i.
//                fij[1] = 24*(2*r12i - r6i)*r2i*r_ij[1];
//                fij[2] = 24*(2*r12i - r6i)*r2i*r_ij[2];

//                f[0] = f[0] + fij[0]; // adding up the forces on particle i in x direction.
//                f[1] = f[1] + fij[1];
//                f[2] = f[2] + fij[2];


//                Ui = Ui + 4*(r12i - r6i); // sum up the potential energy for particle i.

//                force = force + 1;
//            } // end if
//        } // end for j

//        F[i] = f;   // total force in [x,y,z] on particle i
//        U.push_back(Ui); // total potential energy
//    } // end for i


    //cout << "------Forces on particle i------------" << endl;
//    for (int i = 0; i < N; ++i) {
//        cout << i << " " << F[i][0] << " "  << F[i][1] << " " << F[i][2] << endl;
//    }

} // returns the total unitless potential energy


void integrator(vector < vector < double > > &V,vector < vector < double > > &R,int &N, double Lx, double Ly, double Lz, int N_cells_x,int N_cells_y,int N_cells_z, double Lcx, double Lcy, double Lcz, vector < list < int > > &box_list){

    /* *********************************************************************************************
     * Integrator uses the stable Verlet algorithm to calculate the motion of the particles
     *  using the Lenny-Jones potential to find the force between evry particle.
     *  Creates a new state file: stateXXX.xyz for every loop iteration while t < tmax.
     */

    vector < vector < double > > F (R.size(), vector < double > (3,0)); // Vector that holds the forces on particle i.
    vector < double >  U (R.size(),0); // Vector that holds the potential energy for particle i.
    vector < double > E_system;
    vector < double > Temperature;
    vector < double > Ekin;

    double dt = 0.02;
    int tmax = 61;
    double Ek, E_kin, E_tot, E_tot_system;     // E_tot = E_kin + U
    double E_mean_system, E_quad, E_stdev;
    E_mean_system = 0;
    E_quad = 0;

    vector < double > Time_vec ;

    for (int t=0;t<tmax;++t){
        char filename [20];
        sprintf(filename, "state%03d.txt", t);
        ofstream myfile;
        myfile.open(filename);
        myfile << N << endl;
        myfile << filename << "time: " << t*dt << endl;
        for (int boxnr = 0; boxnr<box_list.size(); ++boxnr){
            for (auto it=box_list[boxnr].begin(); it != box_list[boxnr].end(); ++it){
                int index = *it;
                if (floor(index) != index){
                    cout << "somethings wrong!!!" << boxnr << " " << index << endl;
                }
                myfile << "Ar" << " " << R[index][0] << " " << R[index][1] << " " << R[index][2] << " " << V[index][0] << " " << V[index][1] << " " << V[index][2] << " " <<  F[index][0] << " " << F[index][1] << " " << F[index][2] << " " << boxnr << endl;
            }
        }
        myfile.close();
        E_tot_system = 0;
        for (int i = 0; i < N; ++i) {      // update velocity and positions from the forces acting on the particles
            V[i][0] = V[i][0] + F[i][0]*dt/(2*m);   // Calculate V[i] at (t + dt/2)
            V[i][1] = V[i][1] + F[i][1]*dt/(2*m);
            V[i][2] = V[i][2] + F[i][2]*dt/(2*m);

            R[i][0] = R[i][0] + V[i][0]*dt;         // Calculate R[i] at (t + dt)
            R[i][1] = R[i][1] + V[i][1]*dt;
            R[i][2] = R[i][2] + V[i][2]*dt;

            // Adjust positions after periodic boundary conditions
            if (R[i][0] > Lx){R[i][0] = R[i][0] - Lx;}
            else if (R[i][0] < 0){R[i][0] = R[i][0] + Lx;}

            if (R[i][1] > Ly){R[i][1] = R[i][1] - Ly;}
            else if (R[i][1] < 0){R[i][1] = R[i][1] + Ly;}

            if (R[i][2] > Lz){R[i][2] = R[i][2] - Lz;}
            else if (R[i][2] < 0){R[i][2] = R[i][2] + Lz;}
        }

        update_box_list(Lcx,Lcy,Lcz,N_cells_x,N_cells_y,N_cells_z,R,box_list);
        Ek = 0;
        Lennard_Jones(F,R,U,N,Lx,Ly,Lz,N_cells_x,N_cells_y,N_cells_z,box_list);              // calculate the force at time (t+dt) using the new positions.
        for (int i = 0; i < N; ++i) {
            V[i][0] = V[i][0] + F[i][0]*dt/(2*m);   // then find the velocities at time (t+dt)
            V[i][1] = V[i][1] + F[i][1]*dt/(2*m);
            V[i][1] = V[i][2] + F[i][2]*dt/(2*m);

            E_kin = 0.5*m*sqrt(V[i][0]*V[i][0] + V[i][0]*V[i][0] + V[i][0]*V[i][0]);
            E_tot = E_kin + U[i]; // total energy for particle i
            Ek += E_kin;          // total kinetic energy
            E_tot_system += E_tot;
            // Write to file:
            //myfile << "Ar" << " " << R[i][0] << " " << R[i][1] << " " << R[i][2] << " " << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " <<  F[i][0] << " " << F[i][1] << " " << F[i][2] << " " << E_tot << " " << endl;
        }


        Time_vec.push_back(t*dt/Time_0);           // [fs]
        E_system.push_back(E_tot_system);          // Energy of the system.
        Ekin.push_back(Ek);
        double tempi = 2*Ek/(3.0*N);
        Temperature.push_back(tempi);     // Temperature
        cout << "t= " << t << " E= " << E_tot_system << " T= " << tempi << endl;
        //cout << "Total enegy of the system= " << E_tot_system << " at time t= " << t*dt/Time_0 << endl;
    }

    // Write temperatures to file temperatures.txt
    char tempr [20];
    sprintf(tempr,"temperatures.txt");
    ofstream ofile;
    ofile.open(tempr);
    ofile << "Temperature of system as a function of time" << endl;
    ofile << "[Temprature,K] [Time,fs] [E_kin]" << endl;
    for (int i = 0; i < N; ++i) {
        ofile << Temperature[i]*T_0 << " " << i*dt/Time_0 << " " << Ekin[i]*eps << endl;
    }
    ofile.close();

    // calculate the total energy of the system and its standard deviation:
    for (int index = 0; index < N; ++index) {
        E_mean_system = E_mean_system + E_system[index];
        E_quad = E_quad + E_system[index]*E_system[index];
    }
    E_mean_system = E_mean_system/float(N);
    E_quad = E_quad/float(N);
    //cout << E_mean_system << endl;
    //cout << E_quad << endl;
    E_stdev = sqrt(E_quad - E_mean_system*E_mean_system);
    cout << "Total energy of the system= " << E_mean_system << "+-" << E_stdev << endl;
    cout << "fluctuation is " << 100*E_stdev/E_mean_system << "%" << endl;
}




int main(){

    cout << "scaled mass:   " << m << endl;
    cout << "scaled time:   " << 0.02 << endl;
    cout << "scaled length: " << length << endl;
    cout << "scaled Force:  " << eps*sigma/F_0 << endl;  // ? eps, kB ?
    cout << "scaled Energy: " << eps/E << endl;
    cout << "scaled Temp.:  " << T_0 << endl;

    // create a vector within a vector, using including the <vector> library.
    vector < vector < double > > R; // positions
    vector < vector < double > > V; // velocities
    int N;
    double Nx, Ny, Nz; // number of origins
    double Lx,Ly,Lz;   // lattice length
    Nx = 5;
    Ny = 5;
    Nz = 5;
    Lx = Nx*length;
    Ly = Ny*length;
    Lz = Nz*length;

    // Cells
    double Lcx,Lcy,Lcz; // length of cell
    double N_cells_x,N_cells_y,N_cells_z;        // number of cells
    N_cells_x = N_cells_y = N_cells_z = 3.0;
    int N_boxes;
    N_boxes = int(N_cells_x*N_cells_y*N_cells_z);

    Lcx = Lx/N_cells_x;
    Lcy = Ly/N_cells_y;
    Lcz = Lz/N_cells_z;
    // endsure that the total length of the cells equals the length dimentions of the box.
    N_cells_x = Lx/Lcx;
    N_cells_y = Ly/Lcy;
    N_cells_z = Lz/Lcz;
    Lcx = Lx/N_cells_x;
    Lcy = Ly/N_cells_y;
    Lcz = Lz/N_cells_z;

    vector < list < int > > box_list(N_cells_x*N_cells_y*N_cells_z);

    initialize(V,R,N,Nx,Ny,Nz);

    initialize_box_list(Lcx,Lcy,Lcz,N_cells_x,N_cells_y,N_cells_z,R,box_list);

    integrator(V,R,N,Lx,Ly,Lz,N_cells_x,N_cells_y,N_cells_z,Lcx,Lcy,Lcz,box_list);

    return 0;
}

//r_ij[0] = rx - R[aj][0];
//r_ij[1] = ry - R[aj][1];
//r_ij[2] = rz - R[aj][2];

// Periodic boundary conditions
//if (r_ij[0] > Lx/2){r_ij[0] = - Lx + r_ij[0];}
//else if (r_ij[0] < -Lx/2){r_ij[0] = Lx + r_ij[0];}

//if (r_ij[1] > Ly/2){r_ij[1] = -Ly + r_ij[1];}
//else if (r_ij[1] < -Ly/2){r_ij[1] = Ly + r_ij[1];}

//if (r_ij[2] > Lz/2){r_ij[2] = -Lz + r_ij[2];}
//else if (r_ij[2] < -Lz/2){r_ij[2] = Lz + r_ij[2];}

//r2 = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
//r2i = 1.0/r2;
//r6i = r2i*r2i*r2i;
//r12i = r6i*r6i;

//fij[0] = 24*(2*r12i - r6i)*r2i*r_ij[0];  // force from j on i.
//fij[1] = 24*(2*r12i - r6i)*r2i*r_ij[1];
//fij[2] = 24*(2*r12i - r6i)*r2i*r_ij[2];

//f[0] = f[0] + fij[0]; // adding up the forces on particle i in x direction.
//f[1] = f[1] + fij[1];
//f[2] = f[2] + fij[2];


//Ui = Ui + 4*(r12i - r6i); // sum up the potential energy for particle i.

//force = force + 1;
