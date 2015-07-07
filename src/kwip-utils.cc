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

} // end namespace kwip
