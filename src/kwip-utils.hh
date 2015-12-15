/*
 * ============================================================================
 *
 *       Filename:  kwip-utils.hh
 *    Description:  Utilities for kwip
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */


#ifndef KWIP_UTILS_HH
#define KWIP_UTILS_HH

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>

using Eigen::MatrixXd;

#include <kwip-config.hh>

namespace kwip
{

void print_version();

void load_lsmat(MatrixXd &mat, const std::string &filename);
void print_lsmat(MatrixXd &mat, std::ostream &outstream);


template <typename eltype>
eltype
vec_min(std::vector<eltype> &vec)
{
    return *(std::min_element(vec.begin(), vec.end()));
}



} // end namespace kwip

#endif /* KWIP_UTILS_HH */
