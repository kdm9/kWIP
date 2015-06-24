/*
 * ============================================================================
 *
 *       Filename:  kmerclust.cc
 *    Description:  Utilities for kmerclust
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include "kwip.hh"

namespace kmerclust
{

const std::string kmerclust_version = KMERCLUST_VERSION;

void
print_version()
{
    using namespace std;
    cerr << "kmerclust version " << kmerclust_version << endl;
}

} // end namespace kmerclust
