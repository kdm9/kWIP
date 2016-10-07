#include "kwip_config.h"

#include <stdio.h>
#include <string.h>

#include <greatest.h>

#include "test_array.c"
#include "test_kmercount.c"


GREATEST_MAIN_DEFS();

int main(int argc, char *argv[])
{
    GREATEST_MAIN_BEGIN();
    RUN_SUITE(array);
    RUN_SUITE(counting);
    GREATEST_MAIN_END();
    return 0;
}
