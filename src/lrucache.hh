/*
 * File:   lrucache.hpp
 * Author: Alexander Ponomarev
 *
 * Created on June 20, 2013, 5:09 PM
 */

#ifndef _LRUCACHE_HPP_INCLUDED_
#define	_LRUCACHE_HPP_INCLUDED_

#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>

#include <iostream>

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
        //std::cerr << "added " << key << std::endl;
        //this->print("start of put");
    	if (it != _item_map.end()) {
            return;
    	}
    	/* New things at the front, stale things at the back */
    	_key_rc_list.push_front(key_rc_pair_tp(key, 0));
    	_item_map[key] = val_iter_pair_tp(value, _key_rc_list.begin());

        auto rc_itr = _key_rc_list.rbegin();
        //this->print("put, before trim");
        while (_item_map.size() > _capacity && \
               rc_itr != _key_rc_list.rend()) {
            // check ref counter
            //std::cerr << rc_itr->first << " " << rc_itr->second << std::endl;
            if (rc_itr->second == 0) {
                //delete key_rc.second;
                _item_map.erase(rc_itr->first);
                _key_rc_list.erase(std::next(rc_itr).base());
            }
            rc_itr++;
        }
        //this->print("put, after trim");
    }

    void print(const std::string &msg)
    {
        //std::cerr << msg << std::endl;
        for (auto it = _key_rc_list.begin(); it != _key_rc_list.end(); it++) {
            //std::cerr << "\tk:" << it->first << ", rc:" << it->second << std::endl;
        }
    }

    const value_tp& get(const key_tp& key) {
        // Obtain the value for key, and increment its reference counter
    	auto it = _item_map.find(key);
        //std::cerr << "getting " << key << std::endl;
        //this->print("start of get");
    	if (it == _item_map.end()) {
    		throw std::range_error("There is no such key in cache");
    	} else {
    		_key_rc_list.splice(_key_rc_list.end(), _key_rc_list,
                                it->second.second);
            int *rc = &it->second.second->second;
            __sync_fetch_and_add(rc, 1);
            //this->print("end of get");
    		return it->second.first;
    	}
    }

    void unget(const key_tp& key) {
        // Mark the reference as finished, decrementing the ref counter
    	auto it = _item_map.find(key);
        // this->print("start of unget");
        //std::cerr << "un-getting " << key << std::endl;
    	if (it == _item_map.end()) {
    		throw std::range_error("There is no such key in cache");
    	}
        int *rc = &it->second.second->second;
        __sync_fetch_and_sub(rc, 1);
        //this->print("end of un get");
    }


    bool exists(const key_tp& key) const {
    	return _item_map.find(key) != _item_map.end();
    }

    size_t size() const {
    	return _item_map.size();
    }

    void clear()
    {
        //this->print("before clear");
        _item_map.clear();
        _key_rc_list.clear();
        //this->print("after clear");
    }

};

#endif	/* _LRUCACHE_HPP_INCLUDED_ */
