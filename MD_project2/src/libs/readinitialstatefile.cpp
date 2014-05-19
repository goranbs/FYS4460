#include "readinitialstatefile.h"

ReadInitialStateFile::ReadInitialStateFile(string &filename){
    double N;

    ifstream myfile;
    myfile.open(filename.c_str());
    if (myfile.is_open()) {

        string firstline,secondline;
        //getline(myfile,firstline);
        myfile >> N;

        getline(myfile,firstline);
        getline(myfile,secondline);

        vector <double> r (3);
        vector <double> v (3);
        vector <double> f (3);
        vector <double> f_ (3);
        string atomType;
        int box_index;

        while(!myfile.eof()){
            myfile >> atomType;

            for (int i = 0; i < 3; i++) myfile >> r[i];
            for (int i = 0; i < 3; i++) myfile >> v[i];
            for (int i = 0; i < 3; i++) myfile >> f[i];
            myfile >> box_index;

            if (f[0] == f_[0]){
                if (f[1] == f_[1]){
                    if (f[2] == f_[2]){
                        break;
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                f_[i] = f[i];
            }
            Atom argon(r,v,f,0.0);
            atoms.push_back(argon);

            //  myfile.ignore(999, '\n'); // ignore characters untill u find \n, max 999 characters.
        }

        myfile.close();
    }
    else{
        cout << "Error opening file: " << filename << endl;
    }

}

vector <Atom> ReadInitialStateFile::ReturnInitialState(){
    return atoms;
}
