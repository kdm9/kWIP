/*
 * Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef D2FREQ_HH
#define D2FREQ_HH


#include "kernel.hh"

namespace kwip
{
namespace oldmetrics
{

class KernelD2freq : public Kernel
{
public:
    float kernel               (khmer::CountingHash        &a,
                                khmer::CountingHash        &b);

    const std::string       blurb =
            "D2Freq Kernel\n"
            "\n"
            "This kernel calculates the D2 Frequency kernel, the inner\n"
            "product of bin frequencies (bin value / total number of kmers).\n"
            "\n"
            "D2freq = sum (A_f * B_f)\n"
            "Where:\n"
            "A_f, B_f = X / sum X, X=A,B\n";
};

}} // end namespace kwip::oldmetrics

#endif /* D2FREQ_HH */
