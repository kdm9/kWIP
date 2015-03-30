/*
 * Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
 * Original Source Copyright 2014, Alexander Ponomarev (lamerman)
 *
 * This source is a re-write of https://github.com/lamerman/cpp-lru-cache,
 * adding the reference counting code and has been re-licensed to GPL3. The
 * 3-clause BSD licensed orignal source is available via the above link.
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

#ifndef LRUCACHE_HH
#define LRUCACHE_HH

#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>
#include <iostream>

namespace refcounted_lru_cache
{

template<typename key_tp, typename value_tp>
class lru_cache
{
public:
    typedef typename std::pair<key_tp, int>
        key_rc_pair_tp;
    typedef typename std::list<key_rc_pair_tp>
        key_rc_list_tp;
    typedef typename key_rc_list_tp::iterator
        key_rc_list_iter_tp;
    typedef typename std::pair<value_tp, key_rc_list_iter_tp>
        val_iter_pair_tp;
    typedef typename std::unordered_map<key_tp, val_iter_pair_tp>
        item_map_tp;

private:
    key_rc_list_tp _key_rc_list;
    item_map_tp _item_map;
    size_t _capacity;

public:

    lru_cache(size_t cap) : _capacity(cap) { };

    void put(const key_tp& key, const value_tp& value)
    {
        auto it = _item_map.find(key);
        if (it != _item_map.end()) {
            return;
        }
        // New things at the front, stale things at the back
        _key_rc_list.push_front(key_rc_pair_tp(key, 0));
        _item_map[key] = val_iter_pair_tp(value, _key_rc_list.begin());

        auto rc_itr = _key_rc_list.rbegin();
        while (_item_map.size() > _capacity && \
               rc_itr != _key_rc_list.rend()) {
            // check ref counter
            if (rc_itr->second == 0) {
                _item_map.erase(rc_itr->first);
                _key_rc_list.erase(std::next(rc_itr).base());
            }
            rc_itr++;
        }
    }

    const value_tp& get(const key_tp& key)
    {
        // Obtain the value for key, and increment its reference counter
        auto it = _item_map.find(key);
        if (it == _item_map.end()) {
            throw std::range_error("There is no such key in cache");
        } else {
            _key_rc_list.splice(_key_rc_list.end(), _key_rc_list,
                                it->second.second);
            int *rc = &it->second.second->second;
            __sync_fetch_and_add(rc, 1);
            return it->second.first;
        }
    }

    void unget(const key_tp& key)
    {
        // Mark the reference as finished, decrementing the ref counter
        auto it = _item_map.find(key);
        if (it == _item_map.end()) {
            throw std::range_error("There is no such key in cache");
        }
        int *rc = &it->second.second->second;
        __sync_fetch_and_sub(rc, 1);
    }


    bool exists(const key_tp& key) const
    {
        return _item_map.find(key) != _item_map.end();
    }

    size_t size() const
    {
        return _item_map.size();
    }

    void clear()
    {
        _item_map.clear();
        _key_rc_list.clear();
    }

};

} // namespace refcounted_lru_cache
#endif /* LRUCACHE_HH */
