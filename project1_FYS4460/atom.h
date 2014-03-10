#ifndef ATOM_H
#define ATOM_H

#include <vector>

using namespace std;

class Atom{

public:
    Atom(vector <double> r_, vector <double> v_, vector <double> f_, double u_);
    void update_position(vector <double> r_);
    void update_velocity(vector <double> v_);
    void update_force(vector <double> f_);
    void update_potential(double u_);
    vector <double> position();
    vector <double> velocity();
    vector <double> force();
    vector <double> get_initial_position();
    double get_potential();
    void cross_boundary(int i, int j , int k);
    vector <double> return_n_crossings();
    vector <double> return_distance_traveled();

private:
    void distance_traveled(vector <double> r);
    vector <double> r;            // position of atom
    vector <double> v;            // velocity
    vector <double> f;            // total force felt from other particles
    vector <double> r0;     // initial position
    double u;                     // total potential
    vector <double> n_crossings ;
    vector <double> dist;
    int N;

};

#endif // ATOM_H
