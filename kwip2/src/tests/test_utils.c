#include "kwip_utils.h"
#include <math.h>


void test_parse_size(void **state)
{
    const char *cases[] = {
        "1e6", "1m", "1M", "1000000",
        "4M", "2k", "2G"
    };
    size_t expectations[] = {
        1000000, 1l<<20, 1l<<20, 1000000,
        4l<<20, 2l<<10, 2l<<30
    };
    const size_t ntest = sizeof(cases) / sizeof(*cases);

    for (size_t i = 0; i < ntest; i++) {
        size_t got = kwip_parse_size(cases[i]);
        assert_int_equal(got, expectations[i]);
    }
}

/*******************************************************************************
*                                    Suite                                    *
*******************************************************************************/

static const struct CMUnitTest suite_utils[] = {
    cmocka_unit_test(test_parse_size),
};
