#include "kwip_kernelcalc.h"

/*******************************************************************************
*                     Declarations for private functions                      *
*******************************************************************************/

char *kerncalc_filename_to_samplename(const char *filename);


TEST test_kerncalc_fn2sn(void)
{
    const char *cases[][2] = {
        {"test.kct", "test"},
        {"test.h5", "test"},
        {"a.sample.name", "a.sample.name"},
        {"1.234.h5", "1.234"}, 
        {"sample_1.234.h5", "sample_1.234"}, 
    };
    const size_t ntest = sizeof(cases) / sizeof(*cases);

    for (size_t i = 0; i < ntest; i++) {
        const char *fname = cases[i][0];
        const char *expect = cases[i][1];
        char *got = kerncalc_filename_to_samplename(fname);
        int res = strcmp(got, expect);
        ASSERT_EQ(res, 0);
        free(got);
    }

    PASS();
}

/*******************************************************************************
*                                    Suite                                    *
*******************************************************************************/


SUITE(kerncalc)
{
    RUN_TEST(test_kerncalc_fn2sn);
}
