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

#ifndef COUNTMIN_HH
#define COUNTMIN_HH


#include <array>
#include <vector>

namespace kwip
{

template<typename val_tp, typename hash_tp=size_t>
class CountMinSketch
{
protected:
    val_tp **_table;
    std::vector<hash_tp> _tablesizes;
    size_t _n_tables;

public:
    CountMinSketch             (std::vector<hash_tp> tablesizes);
    CountMinSketch             ();
    ~CountMinSketch            ();

    val_tp get                 (hash_tp hash);
    void set                   (hash_tp hash,
                                val_tp val);
    void increment             (hash_tp hash,
                                val_tp by);
    void decrement             (hash_tp hash,
                                val_tp by);
};

} // end namespace kwip


#endif /* COUNTMIN_HH */
