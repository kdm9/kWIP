#include "kwip_utils.h"

#include <math.h>

size_t kwip_parse_size(const char *sizestr)
{
    size_t size = 0;
    char *end;
    size = strtoll(sizestr, &end, 0);
    if (*end == 'k' || *end == 'K') size <<= 10;
    else if (*end == 'm' || *end == 'M') size <<= 20;
    else if (*end == 'g' || *end == 'G') size <<= 30;
    else if (*end == 'e' || *end == 'E') {
        size_t exp = strtoll(end+1, NULL, 10);
        if (exp > 0) {
            size *= (size_t)pow(10., (double)exp);
        }
    }

    return size;
}
