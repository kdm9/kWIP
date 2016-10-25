#include "kwip_config.h"

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>


#include "test_utils.c"
#include "test_array.c"
#include "test_kmercount.c"
#include "test_distcalc.c"

int main(int argc, char *argv[])
{
    int ret = 0;
    ret |= cmocka_run_group_tests_name("utils", suite_utils, NULL, NULL);
    ret |= cmocka_run_group_tests_name("array", suite_array, NULL, NULL);
    ret |= cmocka_run_group_tests_name("kmercount", suite_kmercount, NULL, NULL);
    ret |= cmocka_run_group_tests_name("distcalc", suite_distcalc, NULL, NULL);
    return ret;
}
