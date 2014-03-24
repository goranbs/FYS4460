#include <unittest++/UnitTest++.h>
#include <atom.h>
#include <unitconverter.h>
#include <iostream>

TEST(Atom) {
    Atom argon({1,2,3},{1,2,3},{0.1,0.1,0.1},0.0);
    vector <double> i_pos = argon.return_initial_position();
    double truevalue = false;
    if (i_pos[0] == 1){
        if (i_pos[1] == 2){
            if (i_pos[2] == 3){
                truevalue = true;
            }
        }
    }
    CHECK(truevalue == true);
}

TEST(Atom2){
    int N = 10;
    vector <vector <double> > r (N,vector <double> (3,0.0));
    vector <vector <double> > r2 (N,vector <double> (3,0.0));
    vector <vector <double> > v (N,vector <double> (3,0.0));
    vector <vector <double> > f (N,vector <double> (3,0.0));
    vector <vector <double> > r0 (N,vector <double> (3,0.0));
    vector <Atom> atoms;
    vector <double> v2 (3,0.0);
    vector <double> f2 (3,0.0);
    bool grandtruth,R_truevalue,F_truevalue,V_truevalue, u_true;
    grandtruth = true;
    u_true = true;
    R_truevalue = true;
    V_truevalue = true;
    F_truevalue = true;
    double potential = 0;
    double u0 = 0.1;
    int t_max = 10;



    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            r[i][k] = 1.1;
            v[i][k] = 3,14;
            f[i][k] = 666;
        }
        Atom argon(r[i],v[i],f[i],u0);
        atoms.push_back(argon);
    }
    for (int time = 0; time < t_max; ++time) {

        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                r2[i][k] = 0.1*time;
                v[i][k] = 0.5*time;
                f[i][k] = 6*time;
            }
            potential = u0*time;
            atoms[i].update_position(r[i]);
            atoms[i].update_velocity(v[i]);
            atoms[i].update_force(f[i]);
            atoms[i].update_potential(potential);
        }
    }
    for (int i = 0; i < N; ++i) {
        r0[i] = atoms[i].return_initial_position();
    }

    cout << "Before test: "<< endl;
    cout << "Grand truth= " << grandtruth << endl;
    cout << "R_truevalue= " << R_truevalue << endl;
    cout << "F_truevalue= " << F_truevalue << endl;
    cout << "V_truevalue= " << V_truevalue << endl;
    cout << "u_truevalue= " << u_true << endl;

    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (r[i][k]!=r0[i][k]){
                R_truevalue = false;
                continue;
            }
            f2 = atoms[i].position();
            if (f[i][k] != f2[k]){
                F_truevalue = false;
                continue;
            }
            v2 = atoms[i].velocity();
            if (v[i][k] != v2[k]){
                V_truevalue = false;
                continue;
            }

            if (atoms[i].potential() != u0*t_max){
                u_true = false;
                continue;
            }
        }
    }
    cout << "After test: " << endl;
    cout << "R_truevalue= " << R_truevalue << endl;
    cout << "F_truevalue= " << F_truevalue << endl;
    cout << "V_truevalue= " << V_truevalue << endl;
    cout << "u_truevalue= " << u_true << endl;

    if (R_truevalue != false) {
        grandtruth = false;
        if (V_truevalue != false) {
            grandtruth = false;
            if (F_truevalue != false){
                grandtruth = false;
                if (u_true != false) {
                    grandtruth = false;

                }
            }
        }
    }

    cout << "Grand truth= " << grandtruth << endl;
    CHECK(grandtruth == true);
}

TEST(ATOM_Force){
    vector <double> f1 (3,0.0);
    vector <double> f2 (3,0.0);
    vector <double> f3 (3,0.0);
    bool unlike, like;
    unlike = true;
    like = true;
    for (int k = 0; k < 3; ++k) {
        f1[k] = 0.11;
        f2[k] = 0.12;
    }
    Atom argon({1,1,1},{2,2,2},f1,0.0);
    argon.update_force(f2);
    f3 = argon.force();
    for (int k = 0; k < 3; ++k) {
        if (f1[k] != f3[k]) {
            like = false;
        }
    }
    //argon.reset_force();
    //f1 = argon.force();
    for (int k = 0; k < 3; ++k) {
        cout << f1[k] << "== 0.11" << endl;
        cout << f2[k] << "== 0.12" << endl;
    }

    CHECK(unlike != like);

}


int main() {

    return UnitTest::RunAllTests();
}
