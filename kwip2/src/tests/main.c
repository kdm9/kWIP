#include "kwip_config.h"

#include <stdio.h>
#include <string.h>

#include <greatest.h>

#include "test_utils.c"
#include "test_array.c"
#include "test_kmercount.c"
#include "test_kerncalc.c"


GREATEST_MAIN_DEFS();

int main(int argc, char *argv[])
{
    GREATEST_MAIN_BEGIN();
    RUN_SUITE(suite_utils);
    RUN_SUITE(array);
    RUN_SUITE(counting);
    RUN_SUITE(kerncalc);
    GREATEST_MAIN_END();
    return 0;
}
