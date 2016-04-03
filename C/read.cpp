#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "read.h"

void ReadInput(int* mode, int* maxit, int* nsteps, int* nprint, int* isbinary, \
    double* inistep, double* adjust, double* tol, double* dt, double* damp, double* penalty, \
    int* materialtype, double materialprops[5], double gravity[3]) {
    // Read in the user inputs
    std::ifstream infile("input.txt");
    if (!infile.good()){
        std::cerr << "Couldn't open file: input.txt\n";
    }
    std::vector<std::string> inputs;
    std::string line;
    std::vector<int> ints;
    std::vector<double> dous;
    while (std::getline(infile, line)) {
        inputs.push_back(line);
    }
    int k = 0;
    for (int i = 0; i < inputs.size(); ++i) {
        for (int j = 0; j < inputs[i].size(); ++j) {
            if (inputs[i][j] == ':') {
                k++;
                std::string val;
                for (int m = j + 1; m < inputs[i].size(); ++m) val += inputs[i][m];
                if ( k == 1 || k == 2 || k == 5 || k == 9 || k == 11 || k == 13) {
                    ints.push_back(atoi(val.c_str()));
                } else {
                    dous.push_back(atof(val.c_str()));
                }
            }
        }
    }
    // Assign values
    *mode = ints[0];
    *isbinary = ints[1];
    *tol = dous[0];
    *penalty = dous[1];
    gravity[0] = dous[2];
    gravity[1] = dous[3];
    gravity[2] = dous[4];
    *maxit = ints[2];
    *inistep = dous[5];
    *adjust = dous[6];
    *nsteps = ints[3];
    *dt = dous[7];
    *nprint = ints[4];
    *damp = dous[8];
    *materialtype = ints[5];
    materialprops[0] = dous[9];
    materialprops[1] = dous[10];
    materialprops[2] = dous[11];
    materialprops[3] = dous[12];
    materialprops[4] = dous[13];
}

void ReadMesh(int* nsd, int* nn, int* nel, int* nen, double coords[*nn][*nsd], int connect[*nel][*nen], \
    int* bc_size, int bc_num[*bc_size][2], double bc_val[*bc_size], int* load_size, int* load_type, \
    int load_num[*load_size][2], double load_val[*load_size][*load_type], int share[*nn]) {
        // Coordinates
        std::ifstream incoords("coords.txt");
        if (!incoords.good()){
            std::cerr << "Couldn't open file: coords.txt\n";
        }
        incoords >> *nsd;
        incoords >> *nn;
        for (int i = 0; i < *nn; ++i) {
            for (int j = 0; j < *nsd; ++j) {
                incoords >> coords[i][j];
            }
        }
        // Connectivity
        std::ifstream inconnect("connect.txt");
        if (!inconnect.good()){
            std::cerr << "Couldn't open file: connect.txt\n";
        }
        inconnect >> *nel;
        inconnect >> *nen;
        for (int i = 0; i < *nel; ++i) {
            for (int j = 0; j < *nen; ++j) {
                inconnect >> connect[i][j];
            }
        }
        // Essential Boundary Conditions
        std::ifstream inbc("bc.txt");
        if (!inbc.good()){
            std::cerr << "Couldn't open file: bc.txt\n";
        }
        inbc >> *bc_size;
        for (int i = 0; i < *bc_size; ++i) {
            inbc >> bc_num[i][0] >> bc_num[i][1] >> bc_val[i];
        }
        // Natural Boundary Conditions, either traction or pressure
        std::ifstream inload("load.txt");
        if (!inload.good()){
            std::cerr << "Couldn't open file: load.txt\n";
        }
        inload >> *load_size >> *load_type;
        (*load_type) -= 2;
        for (int i = 0; i < *load_size; ++i) {
            for (int j = 0; j < *load_type; ++j) {
                if (j < 2) inload >> load_num[i][j];
                else inload >> load_val[i][j-2];
            }
        }
        for (int i = 0; i < *nel; ++i) {
            for (int j = 0; j < *nen; ++j) {
                share[connect[i][j]-1]++;
            }
        }
}
