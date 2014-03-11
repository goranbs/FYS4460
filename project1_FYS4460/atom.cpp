#include "atom.h"

/*****************************************************
 *                   Class Atom
 * contains information about atom i such as
 * - Atomnumber/index    i          Atom.index(i) = i
 * - Position            r(x,y,x)   Atom.position(i) = r(x,y,z) for atom i
 * - Velocity            v(x,y,z)   Atom.velocity(i) = v(x,y,z) for atom i
 * - Initial position    ri(x,y,z)  Atom.initial(i)  = ri(x,y,z) for atom i
 * - Travel              r^2(t)     Atom.travel(i)   = r(t) - r_initial
 */


Atom::Atom(vector < double > r_, vector < double > v_, vector < double > f_, double u_){
    // construct the Atom object that holds the r,v,f,u,n_crossings,dist and initial position r0.

    dist = vector <double> (3);
    n_crossings = vector <double> (3);
    N = 0;
    u = u_;
    r = r_;
    v = v_;
    f = f_;
    r0 = r;

}

void Atom::cross_boundary(int i, int j, int k){
    /* i = 1, if the atom has crossed the boundary in positive x-direction, i = -1 for negative direction
     * and i is zero if it has not crossed.
     * Same applies for j for y-direction, and k for z. */
    n_crossings[0] += i;
    n_crossings[1] += j;
    n_crossings[2] += k;
}

vector < double > Atom::return_n_crossings(){
    return n_crossings;
}

vector < double > Atom::return_distance_traveled(){
    dist[0] = dist[0];
    dist[1] = dist[1];
    dist[2] = dist[2];
    return dist;
}

void Atom::distance_traveled(vector<double> r){
    dist[0] += r[0] - r0[0];
    dist[1] += r[1] - r0[1];
    dist[2] += r[2] - r0[2];
}

void Atom::update_position(vector <double> r_){
    for (int i = 0; i < 3; ++i) {
        r[i] = r_[i];
    }
    distance_traveled(r);
}
void Atom::update_velocity(vector <double> v_){
    for (int i = 0; i < 3; ++i) {
        v[i] = v_[i];
    }
}
void Atom::update_force(vector <double> f_){
    for (int i = 0; i < 3; ++i) {
        f[i] = f_[i];
    }
}
void Atom::update_potential(double u_){
    u = u_;
}


vector < double > Atom::position(){
    return r;
}
vector < double > Atom::velocity(){
    return v;
}
vector < double > Atom::force(){
    return f;
}
vector < double > Atom::return_initial_position(){
    return r0;
}




