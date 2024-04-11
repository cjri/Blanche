
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>

using namespace std;

#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

struct run_params {
    int verb;
    int seed; //For random number generator
    string seqs_file;
    int dim; //Number of dimensions to try to optimise - default 2
    int fix; //Try to fix errors in sequencing.  Allow for N values in the sequence data
};

struct sparseseq {
    vector<int> locus;
    vector<char> allele;
};

struct gcoord {
    vector<double> co;
};
