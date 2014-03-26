#ifndef ATOM_H
#define ATOM_H

#include <vector>

using namespace std;

class Atom{

public:
    Atom(vector <double> r_, vector <double> v_, vector <double> f_, double u_);
    void update_position(const vector <double> &r_);
    void update_velocity(const vector <double> &v_);
    void update_force(const vector <double> &f_);
    void subtract_force(const vector <double> &f_);
    void add_force(const vector <double> &f_);
    void update_potential(double &u_);
    void add_potential(double &u_);
    void clear_force();
    void clear_potential();
    inline const vector<double> &position();
    inline const vector <double> &velocity();
    inline const vector <double> &force();
    inline const double potential() const;
    inline const vector <double> &return_initial_position();
    void cross_boundary(int i, int j , int k);
    inline const vector<double> &return_n_crossings() const;


private:
    vector <double> r;            // position of atom
    vector <double> v;            // velocity
    vector <double> f;            // total force felt from other particles
    vector <double> r0;           // initial position
    vector <double> n_crossings ; // number of crossings out of the system (+1 pos dir, -1 neg dir.)
    //vector <double> dist;         // distance traveled.
    double u;                     // total potential
    //int N;

};

inline const vector < double > &Atom::position(){               // R. return position
    return r;
}
inline const vector < double > &Atom::velocity(){               // V. return velocity
    return v;
}
inline const vector < double > &Atom::force(){                  // F. return force
    return f;
}
inline const double Atom::potential() const {                   // P. return potential
    return u;
}
inline const vector<double> &Atom::return_initial_position(){   // R0. return initial position
    return r0;
}
inline const vector < double > &Atom::return_n_crossings() const{           // number of crossings out of system oundaries
    return n_crossings;
}


#endif // ATOM_H
