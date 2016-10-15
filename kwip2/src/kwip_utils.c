#include "kwip_utils.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>

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

void kwip_print_banner(FILE *stream)
{
    const char *boldred = "", *boldblu = "", *normal = "";
    if (isatty(fileno(stream))) {
        boldred = "\033[1;32m";
        boldblu = "\033[1;34m";
        normal = "\033[0m";
    }
    //fprintf(stream, "%skWIP: the k-mer Weighted Inner Product ", boldred);
    fprintf(stream, "%skWIP2: kWIP without inner products ", boldred);
    fprintf(stream, "%s(" KWIP_VERSION ")\n\n%s", boldblu, normal);
}
void kwip_print_version(FILE *stream)
{
    fprintf(stream, "kWIP version " KWIP_VERSION "\n");
}

int kwip_mkdirp(const char *dir)
{
    int res = mkdir(dir, S_IRWXU | S_IRWXG);
    if (res == 0) return res;
    else if (errno == EEXIST) return 0;
    else return -1;
}
