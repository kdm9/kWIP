#ifndef KERNELCALC_H_B3ZZENI2
#define KERNELCALC_H_B3ZZENI2

#include "kwip_config.h"

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "kwip_array.h"
#include "kwip_kmercount.h"


typedef int (*kwip_kern_fn_t)(double *outp, const char *file1,
                              const char *file2, void *extra);

typedef struct {
    size_t num_samples;
    size_t num_samples_alloced;
    char **files;
    char **samplenames;
    char *checkpoint_dir;

    // All below here are filled in after finalisation
   
    bool have_finalised;
    size_t num_compares;
    // vectors of length num_compares (n * (n+1)/2), condensed pairwise matrix
    double *kernelvalues;
    bool *havekernel;

    kwip_kern_fn_t kernfunc;
} kwip_kerncalc_t;

typedef int (*kwip_kerncalc_finalise_fn_t)(kwip_kerncalc_t *ctx);


// set num_samples to 0 if unknown
int kerncalc_init(kwip_kerncalc_t *ctx);
int kerncalc_set_threads(kwip_kerncalc_t *ctx, int num_threads);
int kerncalc_add_sample(kwip_kerncalc_t *ctx, const char *filename, const char *samplename);
// call finalise after adding samples, before calling kerncalc_pairwise
int kerncalc_finalise(kwip_kerncalc_t *ctx, kwip_kerncalc_finalise_fn_t prepfunc);

int kerncalc_compute_kernel(kwip_kerncalc_t *ctx, size_t i, size_t j, void *extra);


// Checkpointing serialisation
int kerncalc_save(kwip_kerncalc_t *ctx, const char *filename);
int kerncalc_load(kwip_kerncalc_t *ctx, const char *filename);


void kerncalc_destroy(kwip_kerncalc_t *ctx);


static inline size_t
kernmatrix_ij_to_condensed(size_t i, size_t j)
{
    if (i < j) {
        size_t tmp = i;
        i = j;
        j = tmp;
    }

    return (i * (i + 1) / 2) + j;
}

static inline int
kernmatrix_condensed_to_ij(size_t *i, size_t *j, size_t n)
{
    if (i == NULL || j == NULL) return -1;

    // This solves (i * (i+1)/2) = n - j for i
    size_t i_ = (size_t)(sqrt(8. * n + 1.) - 1.)/2;
    size_t j_ = n - (i_ * (i_ + 1) / 2);
    *i = i_;
    *j = j_;

    return 0;
}

#endif /* end of include guard: KERNELCALC_H_B3ZZENI2 */
