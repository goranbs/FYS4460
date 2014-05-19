#ifndef READINITIALSTATEFILE_H
#define READINITIALSTATEFILE_H

#include <atom.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

class ReadInitialStateFile
{
public:
    ReadInitialStateFile(string &filename);
    vector <Atom> ReturnInitialState();
private:
    vector <Atom> atoms;
};

#endif // READINITIALSTATEFILE_H
