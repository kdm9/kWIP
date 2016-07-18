#include "kc_array.h"

TEST array_round_trip(void)
{
    int res = 0;
    uint16_t *arr = calloc(CHUNKSIZE, sizeof(uint16_t));
    uint16_t *arr_out = NULL;
    size_t size_out = 0;

    for (size_t i = 0; i < CHUNKSIZE; i++) {
        arr[i] = i % (1<<16);
    }

    res = array_save("test.h5", "test", arr, CHUNKSIZE, H5T_NATIVE_UINT16);
    ASSERT_EQ(res, 0);

    res = array_read("test.h5", "test", (void *)&arr_out, &size_out,
                     H5T_NATIVE_UINT16);
    ASSERT_EQ(res, 0);
    ASSERT_EQ(size_out, CHUNKSIZE);

    for (size_t i = 0; i < size_out; i++) {
        ASSERT_EQ(arr[i], arr_out[i]);
    }
    remove("test.h5");
    PASS();
}

TEST array_iter(void)
{
    int res = 0;
    const size_t num_chunks = 3;
    const size_t num_items = num_chunks * CHUNKSIZE;
    uint16_t *arr = calloc(CHUNKSIZE * num_chunks, sizeof(uint16_t));
    uint16_t *arr_out = NULL;
    size_t size_out = 0;
    array_blockiter_t itr;

    for (size_t i = 0; i < num_items; i++) {
        arr[i] = i % (1<<16);
    }

    res = array_save("test.h5", "test", arr, num_items, H5T_NATIVE_UINT16);
    ASSERT_EQ(res, 0);

    res = array_blockiter_init(&itr, "test.h5", "test");
    ASSERT_EQ(res, 0);
    ASSERT_EQ(itr.num_chunks, num_chunks);

    size_t i = 0;
    while (!array_blockiter_done(&itr)) {
        res = array_blockiter_next(&itr, (void *)&arr_out, &size_out, CHUNKSIZE);
        ASSERT_EQ(res, 0);
        /* we have an even number of chunks, so size will always be CHUNKSIZE */
        ASSERT(arr_out != NULL);
        ASSERT_EQ(size_out, CHUNKSIZE);
        for (size_t j = 0; i < num_items && j < CHUNKSIZE; j++, i++) {
            ASSERT_EQ(arr_out[j], i % (1<<16));
        }
    }
    ASSERT_EQ(i, num_items);

    remove("test.h5");
    PASS();
}


SUITE(array)
{
    RUN_TEST(array_round_trip);
    RUN_TEST(array_iter);
}
