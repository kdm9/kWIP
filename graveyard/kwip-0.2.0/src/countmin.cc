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

#include "countmin.hh"

namespace kwip
{

template<typename val_tp, typename hash_tp>
CountMinSketch<val_tp, hash_tp>::
CountMinSketch(std::vector<hash_tp> tablesizes)
{
    _tablesizes = tablesizes;
    _n_tables = tablesizes.size();

    _table = new val_tp[_n_tables];
    for (size_t i = 0; i < _n_tables; i++) {
        _table[i] = new val_tp[tablesizes[i]];
    }
}

template<typename val_tp, typename hash_tp>
CountMinSketch<val_tp, hash_tp>::
CountMinSketch()
{
    _tablesizes();
    _n_tables = 0;
    _table = NULL;
}

template<typename val_tp, typename hash_tp>
CountMinSketch<val_tp, hash_tp>::
~CountMinSketch()
{
    if (_table != NULL) {
        for (size_t i = 0; i < _n_tables; i++) {
            delete[] _table[i];
        }
        delete[] _table;
    }
}

template<typename val_tp, typename hash_tp>
val_tp
CountMinSketch<val_tp, hash_tp>::
get(hash_tp hash)
{
    size_t tab;
    std::vector<val_tp> tab_vals;
    val_tp minimum;

    for (tab = 0; tab < _n_tables; tab++) {
        hash_tp bin = hash % _tablesizes[tab];
        val_tp this_val = _table[tab][bin];
        if (tab == 0 || minimum > this_val) {
            minimum = this_val;
        }
    }
}

template<typename val_tp, typename hash_tp>
void
CountMinSketch<val_tp, hash_tp>::
increment(hash_tp hash, val_tp by)
{
    size_t tab;
    std::vector<val_tp> tab_vals;

    for (tab = 0; tab < _n_tables; tab++) {
        hash_tp bin = hash % _tablesizes[tab];
        __sync_fetch_and_add(&_table[tab][bin], by);
    }
}

template<typename val_tp, typename hash_tp>
void
CountMinSketch<val_tp, hash_tp>::
decrement(hash_tp hash, val_tp by)
{
    size_t tab;
    std::vector<val_tp> tab_vals;

    for (tab = 0; tab < _n_tables; tab++) {
        hash_tp bin = hash % _tablesizes[tab];
        __sync_fetch_and_sub(&_table[tab][bin], by);
    }
}

template<typename val_tp, typename hash_tp>
void
CountMinSketch<val_tp, hash_tp>::
set(hash_tp hash, val_tp val)
{
    size_t tab;
    std::vector<val_tp> tab_vals;

    for (tab = 0; tab < _n_tables; tab++) {
        hash_tp bin = hash % _tablesizes[tab];
        _table[tab][bin] = val;
    }
}


} // end namespace kwip
