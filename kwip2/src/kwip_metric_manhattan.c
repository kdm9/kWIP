#include "kwip_metric_manhattan.h"
#include <stdlib.h>

int
metric_manhattan_kernel(double *outp, const char *file1, const char *file2, void *extra)
{
    if (outp == NULL || file1 == NULL || file2 == NULL) return -1;
    (void) extra;

    int res;
    array_blockiter_t Aitr, Bitr;
    kc_eltype_t *A = NULL, *B = NULL;
    size_t Alen = 0, Blen = 0;


    res = array_blockiter_init(&Aitr, file1, "counts");
    if (res != 0) return res;
    res = array_blockiter_init(&Bitr, file2, "counts");
    if (res != 0) return res;

    double dist = 0.;
    double anorm = 0, bnorm = 0;
    size_t offset = 0;
    while (!array_blockiter_done(&Aitr) && !array_blockiter_done(&Bitr)) {
        res = array_blockiter_next(&Aitr, (void*)&A, &Alen);
        if (res != 0) return res;
        res = array_blockiter_next(&Bitr, (void*)&B, &Blen);
        if (res != 0) return res;

        if (Alen != Blen) return -1;

        // Block-level variables to avoid floating point round error
        double _dist = 0.;
        double _anorm = 0., _bnorm = 0.;

        offset += Alen;

        #pragma omp simd reduction(+: _dist, _anorm, _bnorm)
        for (size_t i = 0; i < Alen; i++) {
            int32_t a = A[i], b = B[i];
            if (a == 0 && b == 0) continue;
            _dist += abs(a - b);
            _anorm += a;
            _bnorm += b;
        }
        dist += _dist;
        anorm += _anorm;
        bnorm += _bnorm;
    }
    double normdist = dist / (anorm + bnorm);
    if (array_blockiter_done(&Aitr) && array_blockiter_done(&Bitr)) {
        // Check that both iterators have finished
        *outp = normdist;
        return 0;
    } else {
        *outp = 0.0;
        return -1;
    }
}

