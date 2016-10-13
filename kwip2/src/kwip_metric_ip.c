#include "kwip_metric_ip.h"

int
metric_ip_kernel(double *outp, const char *file1, const char *file2, void *extra)
{
    if (outp == NULL || file1 == NULL || file2 == NULL) return -1;

    int res;
    array_blockiter_t Aitr, Bitr;
    kc_eltype_t *A = NULL, *B = NULL;
    size_t Alen = 0, Blen = 0;


    res = array_blockiter_init(&Aitr, file1, "counts");
    if (res != 0) return res;
    res = array_blockiter_init(&Bitr, file2, "counts");
    if (res != 0) return res;

    double kernel = 0.;
    while (!array_blockiter_done(&Aitr) && !array_blockiter_done(&Bitr)) {
        res = array_blockiter_next(&Aitr, (void*)&A, &Alen, KWIP_CHUNKSIZE);
        if (res != 0) return res;
        res = array_blockiter_next(&Bitr, (void*)&B, &Blen, KWIP_CHUNKSIZE);
        if (res != 0) return res;

        if (Alen != Blen) return -1;

        for (size_t i = 0; i < Alen; i++) {
            double a = A[i], b = B[i];
            kernel += a * b;
        }
    }
    if (array_blockiter_done(&Aitr) && array_blockiter_done(&Bitr)) {
        // Check that both iterators have finished
        *outp = kernel;
        return 0;
    } else {
        *outp = 0.0;
        return -1;
    }
}

int metric_ip_normalise(double *kernelvalues, void *extra)
{
    return 0;
}
