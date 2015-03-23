#include <iostream>
#include "lrucache.hh"

int
main (int argc, char *argv[])
{
    lru_cache<int, int *> cache(3);
    cache.put(1, new int(1));
    cache.put(2, new int(1));
    cache.put(3, new int(32));
    cache.get(1);
    cache.get(2);
    cache.get(2);
    cache.get(2);
    cache.put(4, new int(4));
    cache.unget(1);
    cache.put(5, new int(5));
    try {
        cache.get(1);
    } catch (std::range_error &e) {
        std::cerr << "oops, we removed it" << std::endl;
        cache.put(1, new int(12));
    }
    cache.clear();
    return EXIT_SUCCESS;
} /* ----------  end of function main  ---------- */
