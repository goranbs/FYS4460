/*
Inductory molecular dynamic model of Argon gass
*/

using namespace std;

#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>

using namespace arma;

#define pi 4*atan(1)

ofstream ofile;

int main(int argc,char* argv[]){
    /*
     * Create the initial Argon grid and write
     * data to file.
     */

    cout << "Hello World!" << endl;

    double **r;
    int N;
    r = new double*[N];
    N = 9;                   // number of atoms in cube grid
    r = new double[N][3];    // r holds positions of the atoms

    for (int i=0;i<N;i++){
        for (int x=0;x<N;x++){
            for (int y=0;y<N;y++){
                for (int z=0;z<N;z++){
                    r[i] = [x,y,Z];
                }
            }
        }
    }


    for (int i=0;i<N;i++){
    cout << r[i] << " " ;
                    }

    return 0;
}
