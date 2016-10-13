#include "kwip_config.h"

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

// Patch in V1.0 cmocka defs for travis ci
#ifndef cmocka_run_group_tests_name
#define cmocka_run_group_tests_name(n, f, s, t) run_tests(t)
#define cmocka_unit_test unit_test
#define CMUnitTest UnitTest
#endif


#include "test_utils.c"
#include "test_array.c"
#include "test_kmercount.c"
#include "test_kerncalc.c"

int main(int argc, char *argv[])
{
    int ret = 0;
    ret |= cmocka_run_group_tests_name("utils", suite_utils, NULL, NULL);
    ret |= cmocka_run_group_tests_name("array", suite_array, NULL, NULL);
    ret |= cmocka_run_group_tests_name("kmercount", suite_kmercount, NULL, NULL);
    ret |= cmocka_run_group_tests_name("kernelcalc", suite_kernelcalc, NULL, NULL);
    return ret;
}
