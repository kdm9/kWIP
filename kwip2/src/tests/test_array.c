#include "kwip_array.h"

void test_array_round_trip(void **suite)
{
    int res = 0;
    const size_t size = KWIP_CHUNKSIZE + 100;
    uint8_t *arr = calloc(size, sizeof(uint8_t));
    uint8_t *arr_out = NULL;
    size_t size_out = 0;

    for (size_t i = 0; i < size; i++) {
        arr[i] = i % (1<<8);
    }

    res = array_save("test.h5", "test", arr, size, H5T_NATIVE_UINT8);
    assert_int_equal(res, 0);

    res = array_read("test.h5", "test", (void *)&arr_out, &size_out,
                     H5T_NATIVE_UINT8);
    assert_int_equal(res, 0);
    assert_int_equal(size_out, size);

    for (size_t i = 0; i < size_out; i++) {
        assert_int_equal(arr[i], arr_out[i]);
    }
    remove("test.h5");
    free(arr);
    free(arr_out);
}

void test_array_iter(void **suite)
{
    int res = 0;
    const size_t num_chunks = 3;
    const size_t num_items = num_chunks * KWIP_CHUNKSIZE + 10000;
    uint8_t *arr = calloc(num_items, sizeof(uint8_t));
    uint8_t *arr_out = NULL;
    size_t size_out = 0;
    array_blockiter_t itr;

    // Check that the chunksize is a multiple of 1MiB
    assert_int_equal(KWIP_CHUNKSIZE % (1<<20), 0);

    for (size_t i = 0; i < num_items; i++) {
        arr[i] = i % (1<<8);
    }

    res = array_save("test.h5", "test", arr, num_items, H5T_NATIVE_UINT8);
    assert_int_equal(res, 0);

    res = array_blockiter_init(&itr, "test.h5", "test");
    assert_int_equal(res, 0);
    assert_int_equal(itr.num_chunks, num_chunks + 1);

    size_t i = 0;
    size_t chunks = 0;
    while (!array_blockiter_done(&itr)) {
        res = array_blockiter_next(&itr, (void *)&arr_out, &size_out, KWIP_CHUNKSIZE);
        assert_int_equal(res, 0);
        /* we have an even number of chunks, so size will always be KWIP_CHUNKSIZE */
        assert_non_null(arr_out);
        if (chunks++ < num_chunks) {
            assert_int_equal(size_out, KWIP_CHUNKSIZE);
        } else {
            assert_int_equal(size_out, 10000);
        }
        for (size_t j = 0; j < size_out; j++, i++) {
            assert_int_equal(arr_out[j], i % (1<<8));
        }
    }
    assert_int_equal(i, num_items);

    remove("test.h5");
    free(arr);
    free(arr_out);
}

static const struct CMUnitTest suite_array[] = {
    cmocka_unit_test(test_array_iter),
    cmocka_unit_test(test_array_round_trip),
};
