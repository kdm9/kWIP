#ifndef KERNELCALC_H_B3ZZENI2
#define KERNELCALC_H_B3ZZENI2

#include "kwip_config.h"

#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>

#include "kwip_array.h"
#include "kwip_kmercount.h"


typedef int (*kwip_kern_fn_t)(double *outp, kc_eltype_t *count1,
                              kc_eltype_t *count2, void *extra,
                              size_t length);

typedef struct {
    size_t num_samples;
    size_t num_samples_alloced;
    char **files;
    char **samplenames;
    char *checkpoint_dir;

    // All below here are filled in after finalisation
   
    bool have_finalised;
    // vector of length n * (n+1)/2, condensed pairwise matrix
    double *kernelvalues;


    kwip_kern_fn_t kernfunc;

    void *extra; // weights, etc
    size_t extra_size; // in bytes
} kwip_kerncalc_t;

typedef int (*kwip_kerncalc_finalise_fn_t)(kwip_kerncalc_t *ctx);


// set num_samples to 0 if unknown
int kerncalc_init(kwip_kerncalc_t *ctx);
int kerncalc_set_threads(kwip_kerncalc_t *ctx, int num_threads);
int kerncalc_add_sample(kwip_kerncalc_t *ctx, const char *filename, const char *samplename);
// call finalise after adding samples, before calling kerncalc_pairwise
int kerncalc_finalise(kwip_kerncalc_t *ctx, kwip_kerncalc_finalise_fn_t prepfunc);

int kerncalc_pair(kwip_kerncalc_t *ctx); // called internally


// Checkpointing serialisation
int kerncalc_save(kwip_kerncalc_t *ctx, const char *filename);
int kerncalc_load(kwip_kerncalc_t *ctx, const char *filename);


void kerncalc_destroy(kwip_kerncalc_t *ctx);

#endif /* end of include guard: KERNELCALC_H_B3ZZENI2 */
