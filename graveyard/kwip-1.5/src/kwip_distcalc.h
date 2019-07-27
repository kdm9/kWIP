#ifndef KWIP_DISTCALC_H_WVHJB4BU
#define KWIP_DISTCALC_H_WVHJB4BU

#include "kwip_config.h"
#include "kwip_utils.h"

#include <stdio.h>
#include <stdbool.h>
#include <math.h>


typedef int (*kwip_dist_fn_t)(double *outp, const char *file1,
                              const char *file2, void *extra);
typedef int (*kwip_norm_fn_t)(double *outp, const char *file, 
                              void *extra);

typedef struct {
    size_t num_samples;
    size_t num_samples_alloced;
    char **files;
    char **samplenames;
    char *checkpoint_dir;

    // All below here are filled in after finalisation
    bool have_finalised;
    void *extra;  // Additional info, e.g. weights

    // vectors of length num_compares (n * (n-1)/2), condensed pairwise matrix
    size_t num_compares;
    double *normvalues;
    double *havenorm;
    double *distvalues;
    bool *havedist;

    kwip_norm_fn_t normfunc;
    kwip_dist_fn_t distfunc;
} kwip_distcalc_t;


typedef int (*kwip_distcalc_finalise_fn_t)(kwip_distcalc_t *ctx, void *extra);


// set num_samples to 0 if unknown
int distcalc_init(kwip_distcalc_t *ctx);

int distcalc_set_metric(kwip_distcalc_t *ctx, kwip_dist_fn_t distfunc, kwip_norm_fn_t normfunc);

int distcalc_add_sample(kwip_distcalc_t *ctx, const char *filename, const char *samplename);

int distcalc_set_checkpoint_dir(kwip_distcalc_t *ctx, const char *dir);

int distcalc_finalise(kwip_distcalc_t *ctx, kwip_distcalc_finalise_fn_t prepfunc, void *prepfunc_extra);


// Computes the idx'th dist value
int distcalc_compute_dist(kwip_distcalc_t *ctx, size_t idx);


// Checkpointing serialisation
int distcalc_save(kwip_distcalc_t *ctx);
//int distcalc_load(kwip_distcalc_t *ctx, const char *filename);


void distcalc_destroy(kwip_distcalc_t *ctx);


static inline size_t
distmatrix_ij_to_condensed(size_t i, size_t j)
{
    if (i == j) return -1;
    if (i < j) {
        size_t tmp = i;
        i = j;
        j = tmp;
    }

    return (i * (i - 1) / 2) + j;
}

static inline int
distmatrix_condensed_to_ij(size_t *i, size_t *j, size_t n)
{
    if (i == NULL || j == NULL) return -1;

    // This solves (i * (i-1)/2) = n - j for i
    size_t i_ = (size_t)(sqrt(8. * n + 1.) + 1.)/2;
    size_t j_ = n - (i_ * (i_ - 1) / 2);
    *i = i_;
    *j = j_;

    return 0;
}

#endif /* end of include guard: KWIP_DISTCALC_H_WVHJB4BU */
