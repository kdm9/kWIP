#include "kwip_metric_ip.h"
#include <stdlib.h>
#include <math.h>

int
metric_ip_kernel(double *outp, const char *file1, const char *file2, void *extra)
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
        res = array_blockiter_next(&Aitr, (void*)&A, &Alen, KWIP_CHUNKSIZE);
        if (res != 0) return res;
        res = array_blockiter_next(&Bitr, (void*)&B, &Blen, KWIP_CHUNKSIZE);
        if (res != 0) return res;

        if (Alen != Blen) return -1;

        // Block-level variables to avoid floating point round error
        double _dist = 0.;
        double _anorm = 0., _bnorm = 0.;

        offset += Alen;
        for (size_t i = 0; i < Alen; i++) {
            int32_t a = A[i], b = B[i];
            if (a == 0 && b == 0) continue;
            _dist += pow(abs(a - b), 2.);
            _anorm += pow(a, 2);
            _bnorm += pow(b, 2);
        }
        dist += _dist;
        anorm += _anorm;
        bnorm += _bnorm;
    }
    double normdist = dist / sqrt(anorm * bnorm);
    if (array_blockiter_done(&Aitr) && array_blockiter_done(&Bitr)) {
        // Check that both iterators have finished
        *outp = normdist;
        return 0;
    } else {
        *outp = 0.0;
        return -1;
    }
}
