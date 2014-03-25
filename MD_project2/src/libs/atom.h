#ifndef ATOM_H
#define ATOM_H

#include <vector>

using namespace std;

class Atom{

public:
    Atom(vector <double> r_, vector <double> v_, vector <double> f_, double u_);
    void update_position(vector <double> &r_);
    void update_velocity(vector <double> &v_);
    void update_force(vector <double> &f_);
    void add_force(vector <double> &f_);
    void update_potential(double &u_);
    void add_potential(double &u_);
    void clear_force();
    void clear_potential();
    const vector <double> position();
    const vector <double> velocity();
    const vector <double> force();
    const double potential();
    const vector <double> return_initial_position();
    double get_potential();
    void cross_boundary(int i, int j , int k);
    vector <double> return_n_crossings();


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

#endif // ATOM_H
