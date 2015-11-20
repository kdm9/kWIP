/*
 * ============================================================================
 *
 *       Filename:  eigencalc.cc
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

MatrixXd
load_lsmat(const char *filename)
{
    ifstream ifp(filename);
    string line;
    size_t N = 0;
    size_t i, j;
    MatrixXd mat;


    while (getline(ifp, line)){
        istringstream iss(line);
        string label;
        if (N == 0) {
            // Deal with the first row.
            while (iss >> label) {
                N++;
            }
            mat.resize(N, N);
            i = 0;
        } else if (i < N) {
            double val = 0.;
            iss >> label;
            for (j = 0; j < N; j++) {
                iss >> val;
                mat(i, j) = val;
            }
            i++;
        }
    }
    return mat;
}

int
main (int argc, char *argv[])
{
    int N = 4;
    if (argc < 2) {
        cout << "USAGE: " << argv[0] << " <matrix>\n";
        return 1;
    }
    MatrixXd lm = load_lsmat(argv[1]);
    VectorXd ev = lm.eigenvalues().real();
    for (size_t i = 0; i < ev.rows(); i++) {
        double eigenv = ev[i];
        if (eigenv < 0.0) {
            cout << "NEGATIVE EIGENVALUE -- ";
        }
        cout << eigenv << endl;
    }

    return 0;
}

