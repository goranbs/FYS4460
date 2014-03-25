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
    // construct the Atom object that holds the r, v, f, u, n_crossings and initial position r0.

    //dist = vector <double> (3);
    n_crossings = vector <double> (3);
    //N = 0;
    u = u_;
    r = r_;
    v = v_;
    f = f_;
    r0 = r_;

}

void Atom::cross_boundary(int i, int j, int k){
    /* i = 1, if the atom has crossed the boundary in positive x-direction, i = -1 for negative direction
     * and i is zero if it has not crossed.
     * Same applies for j for y-direction, and k for z. */
    n_crossings[0] += i;
    n_crossings[1] += j;
    n_crossings[2] += k;
}

vector < double > Atom::return_n_crossings(){           // number of crossings out of system oundaries
    return n_crossings;
}


void Atom::update_position(vector <double> &r_){        // R. reset postion
    for (int i = 0; i < 3; ++i) {
        r[i] = r_[i];
    }
}
void Atom::update_velocity(vector <double> &v_){        // V. reset velocity
    for (int i = 0; i < 3; ++i) {
        v[i] = v_[i];
    }
}
void Atom::update_force(vector <double> &f_){           // F. reset force
    for (int i = 0; i < 3; ++i) {
        f[i] = f_[i];
    }
}
void Atom::subtract_force(vector<double> &f_){          // F. subtract force from existing
    for (int i = 0; i < 3; ++i) {
        f[i] -= f_[i];
    }
}

void Atom::add_force(vector <double> &f_){              // F. add force to existing
    for (int i = 0; i < 3; ++i) {
        f[i] += f_[i];
    }
}
void Atom::clear_force(){                               // F. clear force vector. Force set to zero
    for (int i = 0; i < 3; ++i) {
        f[i]= 0;
    }
}
void Atom::update_potential(double &u_){                // P. reset potential
    u = u_;
}
void Atom::add_potential(double &u_){                    // P. add potential to existing
    u += u_;
}

void Atom::clear_potential(){                           // P. clear potential. Potential set to zero.
    u = 0;
}


const vector < double > Atom::position(){               // R. return position
    return r;
}
const vector < double > Atom::velocity(){               // V. return velocity
    return v;
}
const vector < double > Atom::force(){                  // F. return force
    return f;
}
const double Atom::potential(){                         // P. return potential
    return u;
}
const vector<double> Atom::return_initial_position(){   // R0. return initial position
    return r0;
}




