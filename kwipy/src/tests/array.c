#include "kc_array.h"

TEST array_round_trip(void)
{
    int res = 0;
    uint16_t *arr = calloc(CHUNKSIZE, sizeof(uint16_t));
    uint16_t *arr_out = NULL;
    size_t size_o = 0;

    for (size_t i = 0; i < CHUNKSIZE; i++) {
        arr[i] = i % (1<<16);
    }

    res = array_save("test.h5", "test", arr, CHUNKSIZE, H5T_NATIVE_UINT16);
    ASSERT_EQ(res, 0);

    res = array_read("test.h5", "test", (void *)&arr_out, &size_o,
                     H5T_NATIVE_UINT16);
    ASSERT_EQ(res, 0);
    ASSERT_EQ(size_o, CHUNKSIZE);

    for (size_t i = 0; i < size_o; i++) {
        ASSERT_EQ(arr[i], arr_out[i]);
    }
    remove("test.h5");
    PASS();
}

SUITE(array)
{
    RUN_TEST(array_round_trip);
}
