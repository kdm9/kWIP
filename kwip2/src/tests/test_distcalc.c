#include "kwip_distcalc.h"

/*******************************************************************************
*                     Declarations for private functions                      *
*******************************************************************************/

char *distcalc_filename_to_samplename(const char *filename);


void test_distcalc_fn2sn(void **ctx)
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
        char *got = distcalc_filename_to_samplename(fname);
        assert_string_equal(got, expect);
        free(got);
    }
}

void test_kernmatrix_ij_condensed(void **ctx)
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
            assert_int_equal(n, matrix[i][j]);
        }
    }

    for (size_t i = 0; i < 100; i++) {
        for (size_t j = 0; j < 100; j++) {
            if (i < j) continue; // Lower triangular matrix
            size_t n = kernmatrix_ij_to_condensed(i, j);
            size_t i_, j_;
            int ret = 0;
            ret = kernmatrix_condensed_to_ij(&i_, &j_, n);
            assert_int_equal(ret, 0);
            assert_int_equal(i_, i);
            assert_int_equal(j_, j);
        }
    }
}

/*******************************************************************************
*                                    Suite                                    *
*******************************************************************************/

static const struct CMUnitTest suite_distcalc[] = {
    cmocka_unit_test(test_distcalc_fn2sn),
    cmocka_unit_test(test_kernmatrix_ij_condensed),
};