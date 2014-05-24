#include <unittest++/UnitTest++.h>
#include <atom.h>
#include <unitconverter.h>
#include <iostream>
#include <list>
#include <readinitialstatefile.h>
#include <generatenanoporoussystem.h>
#include <string>
#include <unitconverter.h>

TEST(Atom) {
    vector<double> r = {1,2,3};
    vector<double> v = {1,2,3};
    vector<double> f = {0.1,0.1,0.1};
    Atom argon(r,v,f,0.0);
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
    grandtruth = false;
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
                f[i][k] = 7*time;
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

    /*
    cout << "Before test: "<< endl;
    cout << "Grand truth= " << grandtruth << endl;
    cout << "R_truevalue= " << R_truevalue << endl;
    cout << "F_truevalue= " << F_truevalue << endl;
    cout << "V_truevalue= " << V_truevalue << endl;
    cout << "u_truevalue= " << u_true << endl;
    */

    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (r[i][k]!=r0[i][k]){
                R_truevalue = false;
                char r1_name [10];
                char r0_name [10];
                sprintf(r1_name, "r[%d][%d]",i,k);
                sprintf(r0_name, "r0[%d][%d]",i,k);
                cout << "r1[i][k] = " << r1_name << " = " << r[i][k] << endl;
                cout << "r0[k] = " << r0_name << " = " << r0[i][k] << endl;
                continue;
            }
            f2 = atoms[i].force();
            if (f[i][k] != f2[k]){
                F_truevalue = false;
                char f1_name [10];
                char f2_name [10];
                sprintf(f1_name, "f1[%d][%d]",i,k);
                sprintf(f2_name, "f2[%d]",k);
                cout << "f1[i][k] = " << f1_name << " = " << f[i][k] << endl;
                cout << "f2[k] = " << f2_name << " = " << f2[k] << endl;
                continue;
            }
            v2 = atoms[i].velocity();
            if (v[i][k] != v2[k]){
                V_truevalue = false;
                char v1_name [10];
                char v2_name [10];
                sprintf(v1_name, "v1[%d][%d]",i,k);
                sprintf(v2_name, "v2[%d]",k);
                cout << "v1[i][k] = " << v1_name << " = " << v[i][k] << endl;
                cout << "v2[k] = " << v2_name << " = " << v2[k] << endl;
                continue;
            }

            if (atoms[i].potential() != u0*(t_max-1)){
                u_true = false;
                cout << "error in the update of the potential value" << endl;
                cout << i << " u0*time = " << u0*t_max << " atoms[i].potential() = " << atoms[i].potential() << endl;
                continue;
            }
        }
    }
    /*
    cout << "After test: " << endl;
    cout << "R_truevalue= " << R_truevalue << endl;
    cout << "F_truevalue= " << F_truevalue << endl;
    cout << "V_truevalue= " << V_truevalue << endl;
    cout << "u_truevalue= " << u_true << endl;
    */
    if (R_truevalue != false) {
        grandtruth = false;
        //cout << "R is true" << endl;
        if (V_truevalue != false) {
            grandtruth = false;
            //cout << "V is true" << endl;
            if (F_truevalue != false){
                grandtruth = false;
                //cout << "F is true" << endl;
                if (u_true != false) {

                    grandtruth = true;

                }
            }
        }
    }

    //cout << "Grand truth= " << grandtruth << endl;
    CHECK(grandtruth == true);
}

TEST(ATOM_Force){
    vector <double> f1 (3,0.0);
    vector <double> f2 (3,0.0);
    vector <double> f3 (3,0.0);
    bool unlike, like;
    unlike = true;  // true
    like = true;    // this should be false, f1 != f3
    for (int k = 0; k < 3; ++k) {
        f1[k] = 0.11;
        f2[k] = 0.12;
    }
    vector<double> r = {1,1,1};
    vector<double> v = {2,2,2};
    Atom argon(r,v,f1,0.0);
    argon.update_force(f2);
    f3 = argon.force();
    for (int k = 0; k < 3; ++k) {
        if (f1[k] != f3[k]) {
            like = false;
        }
    }
    //argon.reset_force();
    //f1 = argon.force();
    /*
    for (int k = 0; k < 3; ++k) {
        cout << f1[k] << "== 0.11" << endl;
        cout << f2[k] << "== 0.12" << endl;
    }
    */

    CHECK(unlike != like);
}

/*
TEST(test2particles){
    double L,T_bath;
    int N;
    vector < list < int > > box_list;

    N = 3;
    L = 1.0;
    T_bath = 0.851;

    test_2particles(L,L,L,N,N,N,box_list,T_bath);
    cout << "test_2particles. See output file" << endl;
    CHECK(2==2);
}
*/
/*
TEST(generateNanoPorousSystem_1){
    // 1) Load thermalized system file
     // 2) create nanoporous system using GenerateNanoPorousSystem
     // 3) create outputfile to be visualized.
     //

    double b = 5.72; // [Å]
    double R0,R1,Lx,Ly,Lz;
    int nSpheres,N;

    R0 = 20; // [Å]
    R1 = 30; // [Å]
    UnitConverter R_0;
    UnitConverter R_1;
    R0 = R_0.from_aangstrom(R0); // [MD units]
    R1 = R_1.from_aangstrom(R1); // [MD units]

    string filename;
    filename = "state0499.txt";
    ReadInitialStateFile initialstate(filename);
    vector <Atom> atoms = initialstate.ReturnInitialState();

    GenerateNanoPorousSystem porosSyst(atoms,R0,R1,Lx,Ly,Lz,nSpheres,N);

}
*/
int main() {

    return UnitTest::RunAllTests();
}
