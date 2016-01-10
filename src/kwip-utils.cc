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
#include <Eigen/Eigenvalues>


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
print_lsmat(MatrixXd &mat, std::ostream &outstream,
            std::vector<std::string> &labels)
{
    for (const auto &label: labels) {
        outstream << "\t" << label;
    }
    outstream << "\n";

    for (size_t i = 0; i < labels.size(); i++) {
        outstream << labels[i];
        for (size_t j = 0; j < labels.size(); j++) {
            outstream << "\t" << mat(i, j);
        }
        outstream << "\n";
    }
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

void
normalise_matrix(MatrixXd &norm, MatrixXd &input)
{
    size_t size = input.rows();
    norm.resize(size, size);
    norm.fill(0.0);
    auto diag = input.diagonal();

    // Normalise the diagonal of the matrix to 1 with an L2 norm
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            norm(i, j) = input(i, j) / sqrt(diag(i) * diag(j));
        }
    }

}

void
kernel_to_distance(MatrixXd &dist, MatrixXd &kernel, bool normalise)
{
    // todo check kernel is square

    size_t size = kernel.rows();
    MatrixXd norm(size, size);
    dist.resize(size, size);
    dist.fill(0.0);

    if (normalise) {
        normalise_matrix(norm, kernel);
    } else {
        norm = kernel;
    }

    // Convert the normalised kernel matrix to distance matrix
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            float d = norm(i, i) + norm(j, j) - 2 * norm(i, j);
            if (d > 0.0) {
                dist(i, j) = sqrt(d);
            } else {
                dist(i, j) = 0.;
            }
        }
    }
}

bool
matrix_is_pos_semidef(MatrixXd &mat)
{
    VectorXd eigenvalues = mat.eigenvalues().real();
    return eigenvalues.minCoeff() > -1e-5;
}

} // end namespace kwip
