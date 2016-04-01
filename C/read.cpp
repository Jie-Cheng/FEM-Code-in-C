#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include "global_variables.h"
#include "read.h"

void ReadInput() {
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
    mode = ints[0];
    isbinary = ints[1];
    tol = dous[0];
    penalty = dous[1];
    gravity[0] = dous[2];
    gravity[1] = dous[3];
    gravity[2] = dous[4];
    maxit = ints[2];
    inistep = dous[5];
    adjust = dous[6];
    nsteps = ints[3];
    dt = dous[7];
    nprint = ints[4];
    damp = dous[8];
    materialtype = ints[5];
    materialprops[0] = dous[9];
    materialprops[1] = dous[10];
    materialprops[2] = dous[11];
    materialprops[3] = dous[12];
    materialprops[4] = dous[13];
}

void ReadMesh() {
        // Coordinates
        std::ifstream incoords("coords.txt");
        if (!incoords.good()){
            std::cerr << "Couldn't open file: coords.txt\n";
        }
        incoords >> nsd;
        incoords >> nn;
        coords = (double**)calloc(nn, sizeof(double*));
        for (int i = 0; i < nn; ++i) {
            coords[i] = (double*)calloc(nsd, sizeof(double));
            for (int j = 0; j < nsd; ++j) {
                incoords >> coords[i][j];
            }
        }
        // Connectivity
        std::ifstream inconnect("connect.txt");
        if (!inconnect.good()){
            std::cerr << "Couldn't open file: connect.txt\n";
        }
        inconnect >> nel;
        inconnect >> nen;
        connect = (int**)calloc(nel, sizeof(int*));
        for (int i = 0; i < nel; ++i) {
            connect[i] = (int*)calloc(nen, sizeof(int));
            for (int j = 0; j < nen; ++j) {
                inconnect >> connect[i][j];
            }
        }
        // Essential Boundary Conditions
        std::ifstream inbc("bc.txt");
        if (!inbc.good()){
            std::cerr << "Couldn't open file: bc.txt\n";
        }
        inbc >> bc_size;
        bc_num = (int**)calloc(bc_size, sizeof(int*));
        bc_val = (double*)calloc(bc_size, sizeof(double));
        for (int i = 0; i < bc_size; ++i) {
            bc_num[i] = (int*)calloc(2, sizeof(int));
            inbc >> bc_num[i][0] >> bc_num[i][1] >> bc_val[i];
        }
        // Natural Boundary Conditions, either traction or pressure
        std::ifstream inload("load.txt");
        if (!inload.good()){
            std::cerr << "Couldn't open file: load.txt\n";
        }
        int load_type;
        inload >> load_size >> load_type;
        load_type -= 2;
        load_num = (int**)calloc(load_size, sizeof(int*));
        load_val = (double**)calloc(load_size, sizeof(double*));
        for (int i = 0; i < load_size; ++i) {
            load_num[i] = (int*)calloc(2, sizeof(int));
            load_val[i] = (double*)calloc(load_type, sizeof(double));
            for (int j = 0; j < load_type; ++j) {
                if (j < 2) inload >> load_num[i][j];
                else inload >> load_val[i][j-2];
            }
        }
}
