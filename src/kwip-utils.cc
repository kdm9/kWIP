/*
 * ============================================================================
 *
 *       Filename:  kwip-utils.cc
 *    Description:  Utilities for kwip
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */


#include "kwip-utils.hh"

namespace kwip
{

const std::string kwip_version = KWIP_VERSION;

void
print_version()
{
    using namespace std;
    cerr << "kwip version " << kwip_version << endl;
}

void
print_lsmat(MatrixXd &mat, std::ostream &outstream)
{

}

void
load_lsmat(MatrixXd &mat, const std::string &filename)
{
    using namespace std;

    ifstream ifp(filename);
    string line;
    size_t N = 0;
    size_t i = 0, j = 0;

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
}

} // end namespace kwip
