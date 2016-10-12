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

TEST test_kernmatrix_ij_condensed(void)
{
    const int8_t matrix[6][6] = {
        {0 ,1 ,3 ,6 ,10,15},
        {1 ,2 ,4 ,7 ,11,16},
        {3 ,4 ,5 ,8 ,12,17},
        {6 ,7 ,8 ,9 ,13,18},
        {10,11,12,13,14,19},
        {15,16,17,18,19,20}
    };
    const size_t N = 6;

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            size_t n = kernmatrix_ij_to_condensed(i, j);
            ASSERT_EQ(n, matrix[i][j]);
        }
    }

    for (size_t i = 0; i < 100; i++) {
        for (size_t j = 0; j < 100; j++) {
            if (i < j) continue; // Lower triangular matrix
            size_t n = kernmatrix_ij_to_condensed(i, j);
            size_t i_, j_;
            int ret = 0;
            ret = kernmatrix_condensed_to_ij(&i_, &j_, n);
            ASSERT_EQ(ret, 0);
            ASSERT_EQ(i_, i);
            ASSERT_EQ(j_, j);
        }
    }

    PASS();
}

/*******************************************************************************
*                                    Suite                                    *
*******************************************************************************/


SUITE(kerncalc)
{
    RUN_TEST(test_kerncalc_fn2sn);
    RUN_TEST(test_kernmatrix_ij_condensed);
}
