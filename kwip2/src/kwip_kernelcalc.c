#include "kwip_kernelcalc.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_eigen.h>


int
kerncalc_init(kwip_kerncalc_t *ctx)
{
    if (ctx == NULL) return -1;
    memset(ctx, 0, sizeof(*ctx));
    return 0;
}

void
kerncalc_destroy(kwip_kerncalc_t *ctx)
{
    if (ctx == NULL) return;
    for (size_t i = 0; i < ctx->num_samples; i++) {
        kwip_free(ctx->files[i]);
        kwip_free(ctx->samplenames[i]);
    }
    kwip_free(ctx->files);
    kwip_free(ctx->samplenames);
    kwip_free(ctx->checkpoint_dir);
    kwip_free(ctx->kernelvalues);
}


char *
kerncalc_filename_to_samplename(const char *filename) {
    if (filename == NULL) return NULL;
    const char *start = strrchr(filename, '/');
    if (start == NULL) {
        start = filename;
    } else {
        start += 1;
    }

    size_t name_len = strlen(start);
    // Strip known extensions, leave others
    const char *extensions[] = {".kct", ".h5"};
    const size_t n_ext = sizeof(extensions) / sizeof(*extensions);
    for (size_t i = 0; i < n_ext; i++) {
        const char *found = NULL;
        found = strstr(start, extensions[i]);
        if (found != NULL) {
            name_len = found - start;
            break;
        }
    }
    return strndup(start, name_len);
}

int
kerncalc_add_sample(kwip_kerncalc_t *ctx, const char *filename, const char *samplename)
{
    if (ctx == NULL || filename == NULL) return -1;

    if (ctx->have_finalised) {
        return -1; // error if calculation in progress
    }

    // reallocate if out of space
    if (ctx->num_samples + 1 > ctx->num_samples_alloced) {
        size_t nalloc = ctx->num_samples_alloced;
        nalloc = nalloc > 0 ? nalloc * 2 : 8;
        ctx->files = realloc(ctx->files, nalloc * sizeof(*ctx->files));
        if (ctx->files == NULL) return -1;
        ctx->samplenames = realloc(ctx->samplenames, nalloc * sizeof(*ctx->samplenames));
        if (ctx->samplenames == NULL) return -1;
        ctx->num_samples_alloced = nalloc;
    }

    ctx->files[ctx->num_samples] = strdup(filename);
    assert(ctx->files[ctx->num_samples] != NULL);
    if (samplename != NULL) {
        ctx->samplenames[ctx->num_samples] = strdup(samplename);
    } else {
        ctx->samplenames[ctx->num_samples] = kerncalc_filename_to_samplename(filename);
    }
    if (ctx->samplenames[ctx->num_samples] == NULL) return -1;
        
    ctx->num_samples += 1;
    return 0;
}


// call finalise after adding samples, before calling kerncalc_pairwise
int
kerncalc_finalise(kwip_kerncalc_t *ctx, kwip_kerncalc_finalise_fn_t prepfunc)
{
    assert(ctx != NULL);

    uint64_t n_samp = ctx->num_samples;
    uint64_t n_compares = n_samp * (n_samp + 1) / 2; // Binomial coeff., including diagonal
    ctx->kernelvalues = calloc(n_compares, sizeof(*ctx->kernelvalues));
    if (prepfunc != NULL) {
        int ret = prepfunc(ctx);
        if (ret != 0) {
            return ret;
        }
    }
    ctx->have_finalised = true;
    return 0;
}

/*******************************************************************************
*                             Kernel calculation                              *
*******************************************************************************/

#if 0
static int
calculate_kernel(kwip_kerncalc_t *ctx, size_t row, size_t col)
{
    assert(ctx != NULL);
    assert(row < ctx->num_samples);
    assert(col < ctx->num_samples);
    int res;
    const char *Afile = ctx->files[row], *Bfile = ctx->files[col];
    array_blockiter_t Aitr, Bitr;
    kc_eltype_t *A, *B;
    size_t Alen = 0, Blen = 0;


    res = array_blockiter_init(&Aitr, Afile, "counts");
    if (res != 0) return res;
    res = array_blockiter_init(&Bitr, Bfile, "counts");
    if (res != 0) return res;

    double kernel = 0.;
    size_t offset = 0;
    while (!array_blockiter_done(&Aitr) &&
           !array_blockiter_done(&Bitr)) {
        double block_kern = 0;

        res = array_blockiter_next(&Aitr, (void*)&A, &Alen, KWIP_CHUNKSIZE);
        if (res != 0) return res;
        res = array_blockiter_next(&Bitr, (void*)&B, &Blen, KWIP_CHUNKSIZE);
        if (res != 0) return res;

        if (Alen != Blen) return -1;

        res = ctx->kernfunc(&block_kern, A, B, ctx->extra, ctx->extra_size, offset, Alen);
        if (res != 0) return res;
        offset += Alen;
        kernel += block_kern;
    }

    gsl_matrix_set(ctx->kmat, row, col, kernel);
    gsl_matrix_set(ctx->kmat, col, row, kernel); // symmetric!
    return 0;
}

int
kerncalc_pairwise(kwip_kerncalc_t *ctx)
{
    assert(ctx != NULL);
    int res = 0;

    for (size_t i = 0; i < ctx->num_samples; i++) {
        for (size_t j = 0; j < ctx->num_samples; j++) {
            if (j < i) continue; // upper triangle only!

            res = calculate_kernel(ctx, i, j);
            if (res != 0) return res;
        }
    }
    return 0;
}

int
kerncalc_pairwise_subset(kwip_kerncalc_t *ctx, size_t block, size_t num_blocks)
{
    return 0;
}

int kerncalc_kern2dist(kwip_kerncalc_t *ctx); // called internally

int kerncalc_printdist(kwip_kerncalc_t *ctx, FILE *fp);
int kerncalc_printkernel(kwip_kerncalc_t *ctx, FILE *fp);

void kerncalc_destroy(kwip_kerncalc_t *ctx);
#endif

