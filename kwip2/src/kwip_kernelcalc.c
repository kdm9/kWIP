#include "kwip_kernelcalc.h"

#include <string.h>
#include <assert.h>

#include <gsl/gsl_eigen.h>


int
kerncalc_init(kerncalc_t *ctx, size_t num_samples)
{
    assert(ctx != NULL);

    if (num_samples > 0) {
        ctx->num_samples_alloced = num_samples;
        ctx->files = calloc(num_samples, sizeof(*ctx->files));
        ctx->samplenames = calloc(num_samples, sizeof(*ctx->samplenames));
    }
    ctx->have_finalised = false;
    return 0;
}

int kerncalc_set_threads(kerncalc_t *ctx, int num_threads);

int
kerncalc_add_sample(kerncalc_t *ctx, const char *filename, const char *samplename)
{
    assert(ctx != NULL);
    assert(filename != NULL);

    if (ctx->have_finalised) {
        return -1; // error if calculation in progress
    }

    // reallocate if out of space
    if (ctx->num_samples + 1 > ctx->num_samples_alloced) {
        size_t nalloc = ctx->num_samples_alloced;
        nalloc = nalloc > 0 ? nalloc * 2 : 2;
        ctx->files = realloc(ctx->files, nalloc * sizeof(*ctx->files));
        assert(ctx->files);
        ctx->samplenames = realloc(ctx->samplenames, nalloc * sizeof(*ctx->samplenames));
        assert(ctx->files);
        ctx->num_samples_alloced = nalloc;
    }

    ctx->files[ctx->num_samples] = strdup(filename);
    assert(ctx->files[ctx->num_samples] != NULL);
    if (samplename != NULL) {
        ctx->samplenames[ctx->num_samples] = strdup(samplename);
        assert(ctx->samplenames[ctx->num_samples] != NULL);
    }
    ctx->num_samples += 1;
    return 0;
}


// call finalise after adding samples, before calling kerncalc_pairwise
int
kerncalc_finalise(kerncalc_t *ctx, kerncalc_prepfunc_t prepfunc)
{
    assert(ctx != NULL);

    ctx->kmat = gsl_matrix_calloc(ctx->num_samples, ctx->num_samples);
    if (prepfunc != NULL) {
        if (prepfunc(ctx) != 0) {
            return -1;
        }
    }
    ctx->have_finalised = true;
    return 0;
}

/*******************************************************************************
*                             Kernel calculation                              *
*******************************************************************************/

static int
calculate_kernel(kerncalc_t *ctx, size_t row, size_t col)
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
kerncalc_pairwise(kerncalc_t *ctx)
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
kerncalc_pairwise_subset(kerncalc_t *ctx, size_t block, size_t num_blocks)
{
    return 0;
}

int kerncalc_kern2dist(kerncalc_t *ctx); // called internally

int kerncalc_printdist(kerncalc_t *ctx, FILE *fp);
int kerncalc_printkernel(kerncalc_t *ctx, FILE *fp);

void kerncalc_destroy(kerncalc_t *ctx);

