#ifndef KERNELCALC_H_B3ZZENI2
#define KERNELCALC_H_B3ZZENI2

#include <stdio.h>
#include <stdbool.h>

#include <gsl/gsl_math.h>
#include "kwip_array.h"
#include "kwip_kmercount.h"


typedef int (*kernelfunc_t)(double *outp, kc_eltype_t *count1,
                            kc_eltype_t *count2, void *extra, size_t extra_size,
                            size_t offset, size_t count);

typedef struct {
    size_t num_samples;
    size_t num_samples_alloced;
    char **files;
    char **samplenames;
    char *checkpoint_dir;
    int threads;
    int verbosity;
    FILE *logfp;
    gsl_matrix *kmat;
    bool have_finalised;
    kernelfunc_t kernfunc;
    void *extra; // weights, etc
    size_t extra_size; // in bytes
} kerncalc_t;

typedef int (*kerncalc_prepfunc_t)(kerncalc_t *ctx);

// set num_samples to 0 if unknown
int kerncalc_init(kerncalc_t *ctx, size_t num_samples);
int kerncalc_set_threads(kerncalc_t *ctx, int num_threads);
int kerncalc_add_sample(kerncalc_t *ctx, const char *filename, const char *samplename);
// call finalise after adding samples, before calling kerncalc_pairwise
int kerncalc_finalise(kerncalc_t *ctx, kerncalc_prepfunc_t prepfunc);

int kerncalc_pairwise(kerncalc_t *ctx);
int kerncalc_pairwise_subset(kerncalc_t *ctx, size_t block, size_t num_blocks);

// Checkpointing serialisation
int kerncalc_save(kerncalc_t *ctx, const char *filename);
int kerncalc_load(kerncalc_t *ctx, const char *filename);
int kerncalc_pairwise_checkpoint(kerncalc_t *ctx, const char *tempdir);
int kerncalc_pairwise_subset_checkpoint(kerncalc_t *ctx, size_t block, size_t num_blocks, const char *tempdir);

int kerncalc_kern2dist(kerncalc_t *ctx); // called internally

int kerncalc_printdist(kerncalc_t *ctx, FILE *fp);
int kerncalc_printkernel(kerncalc_t *ctx, FILE *fp);

void kerncalc_destroy(kerncalc_t *ctx);

#endif /* end of include guard: KERNELCALC_H_B3ZZENI2 */
