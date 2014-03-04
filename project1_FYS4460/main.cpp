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
#include <algorithm>        //
#include "initialstate.h"
//#include <armadillo>

using namespace std;
//using namespace arma;

double pi = 4*atan(1);

/**********************************************************************************************************
 *             FUNCTION DECLARATIONS
 * ********************************************************************************************************
 */
void write_to_file(vector<vector<double> > R, vector<vector<double> > V, vector<vector<double> > F, vector<list<int> > box_list, string filename, int N, double t);

void test_2particles(double Lx, double Ly, double Lz, int N_cells_x, int N_cells_y, int N_cells_z, vector<list<int> > box_list);

double calculate_forces(vector<vector<double> > &R, vector<vector<double> > &F, int i, int j, double Lx, double Ly, double Lz, double Pi);

/**********************************************************************************************************
 *             CONSTANTS
 **********************************************************************************************************
 */
// Check out the constants, do they fit with those in the project text?
double b = 5.260;                 // Ångstrøm [Å]
//double b = 20.0;                // Ångsrøm [Å]
double mA = 39.948;               // mass of Argon [amu]
double kB = 1.480*pow(10,-23);    // Bolzmann constant [eV/K]
double eps = 0.01*1.0318;         // Energy constant [eV]
double Temp = 100.0;              // Kelvin, initial temperature
double sigma = 3.405;             // Ångstrøm, scalefactor - Leonard-Jones
double P_0 = 1.60217657*pow(10,11.0)/(b*b*b); // [N/m^2] Pressure
/**********************************************************************************************************
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
     * Standard deviation of the distribution std = sqrt(kB*T/m)
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
        if (box_index > box_list.size()){
            cout <<"Tried to access box nr: " << box_index << endl;
        }
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

    for (int p=0; p<N;++p){
        // get rid of drift:
        V[p][0] = V[p][0] - mean_x;
        V[p][1] = V[p][1] - mean_y;
        V[p][2] = V[p][2] - mean_z;
    }

}

/* ******************************************************************************************************
 *                                          LENNARD-JONES
 *                                         BOX-calculation
 ** ******************************************************************************************************
 */
void Lennard_Jones(vector < vector < double > > &F, vector < vector < double > > &R, vector < double > &U, int N, double Lx, double Ly, double Lz, int N_cells_x, int N_cells_y, int N_cells_z, vector < list < int > > &box_list, double &P_sum){
    /* The Lenny-Jones potential updates the forces F on particle i in position R.
     * One box-calculation - calulating the contribution from every particle in the system.
     *

        for (int i = 0; i < R.size(); ++i) {
            double Ui = 0; // potential energy
            for (int j = i; j < R.size(); ++j) {
                if ( j != i) {

                    Ui += calculate_forces(R,F,i,j,Lx,Ly,Lz);

                } // end if
            } // end for j

            U[i] = Ui; // Total potential energy of particle i.
        } // end for i

} // end Lennard-Jones
*/
/*******************************************************************************************************************
 * forces with boxes - currently not working :-)
 */


    vector < vector < int > > list_of_neighbours_to_calculate_for_each_box (N_cells_x*N_cells_y*N_cells_z);
    vector < int >  list_of_visited_boxes;
    for (int box_x = 0; box_x < N_cells_x; ++box_x){
        for (int box_y = 0; box_y < N_cells_y; ++box_y){
            for (int box_z = 0; box_z < N_cells_z; ++box_z){

                int main_box_index = box_x*N_cells_y*N_cells_z + box_y*N_cells_z + box_z;

                    for (int dx=-1; dx<=1; ++dx){
                        for (int dy=-1; dy<=1; ++dy){
                            for (int dz=-1; dz<=1; ++dz){
                                //neighbour boxes:
                                int neighbour_x = (box_x + dx + N_cells_x) % N_cells_x;
                                int neighbour_y = (box_y + dy + N_cells_y) % N_cells_y;
                                int neighbour_z = (box_z + dz + N_cells_z) % N_cells_z;

                                /*
                                if (box_x < 0){ box_x = box_x + N_cells_x;}
                                else if (box_x >= N_cells_x){ box_x = box_x - N_cells_x;}

                                if (box_y < 0){ box_y = box_y + N_cells_y; }
                                else if (box_y >= N_cells_y){box_y = box_y - N_cells_y;}

                                if (box_z < 0){ box_z = box_z + N_cells_z;}
                                else if (box_z >= N_cells_z){box_z = box_z - N_cells_z;}
                                */

                                int neighbour_index = neighbour_x*N_cells_y*N_cells_z + neighbour_y*N_cells_z + neighbour_z;

                                if ( find(list_of_visited_boxes.begin(), list_of_visited_boxes.end(), neighbour_index) == list_of_visited_boxes.end()){
                                    list_of_neighbours_to_calculate_for_each_box[main_box_index].push_back(neighbour_index);

                                }
                            }
                        }
                    }
                    list_of_visited_boxes.push_back(main_box_index);

            }
        }
    }

    /* the force calculation
     * We need to calulate the forces between every particle in the main box and its neighbours
     */
    for (int x = 0; x < N_cells_x; ++x) {
        for (int y = 0; y < N_cells_y; ++y) {
            for (int z = 0; z < N_cells_z; ++z) {

                int main_box_index = N_cells_x*N_cells_y*x + N_cells_y*y + z;

                for (auto it = box_list[main_box_index].begin(); it!= box_list[main_box_index].end(); ++it){
                    // iterates over all particles in box_list[main_index]
                    int i = *it;
                    double u = 0;

                    for (int neighbour_index : list_of_neighbours_to_calculate_for_each_box[main_box_index]){
                        // for every neighbour:



                        for (auto it2 = box_list[neighbour_index].begin(); it2 != box_list[neighbour_index].end(); ++it2){
                            // for every particle in neighbour box:

                            int j = *it2;

                            if (i != j){

                                u += calculate_forces(R,F,i,j,Lx,Ly,Lz,P_sum);

                            }

                        } // particle in neighbour box
                    } // neighbour box
                    U[i] = u;
                } // main box
            }
        }
    }

}// end Lennard-Jones


/********************************************************************************************************
 *                              CALCULATE FORCES
 * ******************************************************************************************************
 */

double calculate_forces(vector < vector < double > > &R, vector < vector < double > > &F, int i, int j, double Lx,double Ly,double Lz, double Pi){
    /* Takes positionvector R for particle i and j. Where F is the total Force acting on the particle.
     * The function calculate_forces returns the potential energy for particle i felt from particle j.
     * The total potential energy for particle i is then the sum of potentials from all other particles.
     */

    double r2,r2i,r6i,r12i,Ui;
    vector < double > fij (3);
    vector < double > r_ij (3);

    r_ij[0]=  R[i][0] - R[j][0];
    r_ij[1] = R[i][1] - R[j][1];
    r_ij[2] = R[i][2] - R[j][2];

    // Periodic boundary conditions - Minimum Image Convention:
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

    F[i][0] = F[i][0] + fij[0]; // adding up the forces on particle i in x direction.
    F[i][1] = F[i][1] + fij[1];
    F[i][2] = F[i][2] + fij[2];

    F[j][0] = F[j][0] - fij[0]; // adding up the forces on particle j in x direction.
    F[j][1] = F[j][1] - fij[1];
    F[j][2] = F[j][2] - fij[2];

    Ui = 4*(r12i - r6i); // the potential energy for particle i in j's presence.

    if (i < j){ // just to be sure :-)
         Pi += (F[i][0]*r_ij[0] + F[i][1]*r_ij[1] + F[i][2]*r_ij[2]); // contribution to the system pressure
    }
    return Ui;
}


/***********************************************************************************************************************
 *                                  INTEGRATOR
 * *********************************************************************************************************************
 */
void integrator(vector < vector < double > > &V,vector < vector < double > > &R,int N, double Lx, double Ly, double Lz, int N_cells_x,int N_cells_y,int N_cells_z, double Lcx, double Lcy, double Lcz, vector < list < int > > &box_list){

    /* *********************************************************************************************
     * Integrator uses the stable Verlet algorithm to calculate the motion of the particles
     *  using the Lenny-Jones potential to find the force between evry particle.
     *  Creates a new state file: stateXXX.xyz for every loop iteration while t < tmax.
     */

    vector < vector < double > > F (R.size(), vector < double > (3,0)); // Vector that holds the forces on particle i.
    vector < double >  U (R.size(),0);         // Vector that holds the potential energy for particle i.
    vector < double > E_system;
    vector < double > Temperature;
    vector < double > Ekin;
    vector < double > Epot;

    double dt = 0.02;
    int tmax = 50;
    double Ek, Ep;
    double E_mean_system, E_quad, E_stdev;
    E_mean_system = 0;
    E_quad = 0;

    vector < double > Time_vec ;
    vector < double > Pressure ;

    for (int t=0;t<tmax;++t){
        char filename [20];
        sprintf(filename, "state%04d.txt", t);

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
        Ek = 0;
        Ep = 0;
        double P_sum = 0;

        for (int p = 0; p < R.size(); ++p) {
            // clear the force matrix:
            F[p][0] = 0;
            F[p][1] = 0;
            F[p][2] = 0;
            //And the potential energy!!
            U[p] = 0;
            //cout << "clearing force and potential energy..." << endl;
        }

        update_box_list(Lcx,Lcy,Lcz,N_cells_x,N_cells_y,N_cells_z,R,box_list);      // Update box-list
        Lennard_Jones(F,R,U,N,Lx,Ly,Lz,N_cells_x,N_cells_y,N_cells_z,box_list,P_sum);  // calculate the force at time (t+dt) using the new positions.
        write_to_file(R,V,F,box_list,filename,N,t*dt*Time_0);                    // write to file

        for (int i = 0; i < N; ++i) {
            V[i][0] = V[i][0] + F[i][0]*dt/(2*m);   // then find the velocities at time (t+dt)
            V[i][1] = V[i][1] + F[i][1]*dt/(2*m);
            V[i][2] = V[i][2] + F[i][2]*dt/(2*m);

            Ek += 0.5*m*(V[i][0]*V[i][0] + V[i][1]*V[i][1] + V[i][2]*V[i][2]); // total kinetic energy
            Ep += U[i];
        }


        Time_vec.push_back(t*dt/Time_0);           // [fs]
        E_system.push_back(Ek+Ep);                 // Energy of the system.
        Ekin.push_back(Ek);
        Epot.push_back(Ep);
        double tempi = 2*Ek/(3.0*N);
        Temperature.push_back(tempi);              // Temperature
        double Press = (N*tempi + P_sum/3);        // Pressure
        Pressure.push_back(Press);                 // Pressure
        cout << "t= " << t << " E= " << E_system[t] << "  Ekin= " << Ekin[t] << "  U= " << Epot[t] << "  T= " << tempi << " P= " << Press << endl;
        //cout << "Total enegy of the system= " << E_tot_system << " at time t= " << t*dt/Time_0 << endl;
    }

    // Write temperatures to file temperatures.txt
    ofstream ofile("temperatures.txt");
    ofile << "Temperature of system, Kinetic, Potential Energy and Pressure as a function of time," << endl;
    ofile << "[Temprature,K] [Time,fs] [E_kin,eV] [Epot,eV] [P,N/Å^2]" << endl;
    for (int t = 0; t < tmax; ++t) {
        ofile << Temperature[t]*T_0 << " " << t*dt*Time_0 << " " << Ekin[t]*eps <<  " "  << Epot[t]*eps <<  " " << Pressure[t]*P_0 << endl;
    }
    ofile.close();

    // Calculate the total energy of the system and its standard deviation:
    for (int index = 0; index < N; ++index) {
        E_mean_system = E_mean_system + E_system[index];
        E_quad = E_quad + E_system[index]*E_system[index];
        //cout << E_system[index] << endl;
    }
    E_mean_system = E_mean_system/float(N);
    E_quad = E_quad/float(N);
    //cout << E_mean_system << endl;
    //cout << E_quad << endl;
    E_stdev = sqrt(E_quad - E_mean_system*E_mean_system);
    cout << "Total energy of the system= " << E_mean_system << " stdev= " << E_stdev << endl;
    cout << "fluctuation is " << 100*E_stdev/E_mean_system << "%" << endl;

}


//void ReadInitailState(string filename, vector < vector < double > > &R, vector < vector < double > > &V){
/* Read inital state from filename and return positions R, velocities V
 */

//    ifstream myfile;
//    myfile.open(filename.c_str());

//    string firstline, secondline;
//    getline(myfile,firstline);
//    getlinge(myfile,secondline);
//    secondline.c_str(); // split on whitespace
//    double t0 = secondline[2];
//    double rx,ry,rz,vx,vy,vz,fx,fy,fz;
//    int it = 0;
//    while(!myfile.eof()){
//        myfile >> rx;
//        myfile >> ry;
//        myfile >> rz;
//        myfile >> vx;
//        myfile >> vy;
//        myfile >> vz;
//        myfile >> fx;
//        myfile >> fy;
//        myfile >> fz;
//        R[it][0] = rx;
//        R[it][1] = ry;
//        R[it][2] = rz;
//        V[it][0] = vx;
//        V[it][1] = vy;
//        V[it][2] = vz;
//        F[it][0] = fx;
//        F[it][1] = fy;
//        F[it][2] = fz;
//        it += it;
//    }

//}

/*******************************************************************************************************************************
 *                                        Test with 2 particles
 * *****************************************************************************************************************************
 */

void test_2particles(double Lx,double Ly,double Lz,int N_cells_x,int N_cells_y,int N_cells_z, vector < list < int > > box_list){

    vector < vector < double > > R (2,vector < double > (3,0));
    vector < vector < double > > V (2,vector < double > (3,0));

    InitialState twoparts;
    twoparts.two_particles(R,V,Lx,Ly,Lz);

    double Lcx = Lx/N_cells_x;
    double Lcy = Ly/N_cells_y;
    double Lcz = Lz/N_cells_z;
    // endsure that the total length of the cells equals the length dimentions of the box.
    N_cells_x = Lx/Lcx;
    N_cells_y = Ly/Lcy;
    N_cells_z = Lz/Lcz;
    Lcx = Lx/N_cells_x;
    Lcy = Ly/N_cells_y;
    Lcz = Lz/N_cells_z;

    initialize_box_list(Lcx,Lcy,Lcz,3,3,3,R,box_list);

    clock_t time3;
    time3 = clock();
    integrator(V,R,2,Lx,Ly,Lz,N_cells_x,N_cells_y,N_cells_z,Lcx,Lcy,Lcz,box_list);
    time3 = clock() - time3;

    cout << "two particles" << endl;
    cout << "Integrator used time= " << float(time3)/CLOCKS_PER_SEC << " seconds" << endl;

}

/***************************************************************************************************************************
 *                                    Write to File
 ***************************************************************************************************************************
 */
void write_to_file(vector < vector < double > > R, vector < vector < double > > V, vector < vector < double > > F, vector < list < int > > box_list,string filename, int N, double t){
    /* Writes positions R[p], velocities V[p] and forcec F[p] to file where p denotest particle p, to file filename.
     * every particle belongs to a box in the system. The box is numbered boxnumber.
     */

    ofstream myfile;
    myfile.open(filename);
    myfile << N << endl;
    myfile << filename << "time: " << t << endl;

    for (int box_nr = 0; box_nr < box_list.size(); ++box_nr){

        for (auto it = box_list[box_nr].begin(); it != box_list[box_nr].end(); ++it){
            int p = *it; // particle p in boxnr

            myfile << "Ar" << " " << R[p][0] << " " << R[p][1] << " " << R[p][2] << " " << V[p][0] << " " << V[p][1] << " " << V[p][2] << " " <<  F[p][0] << " " << F[p][1] << " " << F[p][2] << " " << box_nr << endl;
        }
    }
    myfile.close();
}


/*******************************************************************************************************************
 *                                                MAIN
 * *****************************************************************************************************************
 */
int main(){

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

    /*******************************************************************************************
     *                 unittesting with only two particles
     */

    //test_2particles(Lx,Ly,Lz,N_cells_x,N_cells_y,N_cells_z, box_list);

    /*******************************************************************************************
     *  using initialize to generate intial state
     */


    clock_t time1, time2;
    time1 = clock();

    initialize(V,R,N,Nx,Ny,Nz);

    time1 = clock() - time1;

    initialize_box_list(Lcx,Lcy,Lcz,N_cells_x,N_cells_y,N_cells_z,R,box_list);

    time2 = clock();
    integrator(V,R,N,Lx,Ly,Lz,N_cells_x,N_cells_y,N_cells_z,Lcx,Lcy,Lcz,box_list);
    time2 = clock() - time2;

    cout << "____________________________________________________________________________________________" << endl;
    cout << "Unitless; mass= " << m << " Energy= " << E << " Temperature= " << Temp/T_0 << " Length= " << length << endl;
    cout << "____________________________________________________________________________________________" << endl;
    cout << "Initialize used time= " << float(time1)/CLOCKS_PER_SEC << " seconds" << endl;
    cout << "Integrator used time= " << float(time2)/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}

