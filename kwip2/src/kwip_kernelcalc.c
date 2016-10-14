#include "kwip_kernelcalc.h"

#include "kwip_array.h"
#include "kwip_kmercount.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

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
    kwip_free(ctx->havekernel);
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
    if (ctx->files[ctx->num_samples] == NULL) return -1;
    if (samplename != NULL) {
        ctx->samplenames[ctx->num_samples] = strdup(samplename);
    } else {
        ctx->samplenames[ctx->num_samples] = kerncalc_filename_to_samplename(filename);
    }
    if (ctx->samplenames[ctx->num_samples] == NULL) return -1;
        
    ctx->num_samples += 1;
    return 0;
}

int kerncalc_set_kernelfunction(kwip_kerncalc_t *ctx, kwip_kern_fn_t kernfunc)
{
    if (ctx == NULL) return -1;

    ctx->kernfunc = kernfunc;
    return 0;
}

int kerncalc_set_checkpoint_dir(kwip_kerncalc_t *ctx, const char *dir)
{
    if (ctx == NULL) return -1;

    ctx->checkpoint_dir = strdup(dir);
    return 0;
}

// call finalise after adding samples, before calling kerncalc_pairwise
int kerncalc_finalise(kwip_kerncalc_t *ctx, kwip_kerncalc_finalise_fn_t prepfunc,
                      void *prepfunc_extra)
{
    if (ctx == NULL) return -1;

    if (ctx->kernfunc == NULL) return -1;

    uint64_t n_samp = ctx->num_samples;
    if (n_samp < 2) return -1;

    uint64_t n_compares = n_samp * (n_samp + 1) / 2; // Binomial coeff., including diagonal
    ctx->kernelvalues = calloc(n_compares, sizeof(*ctx->kernelvalues));
    if (ctx->kernelvalues == NULL) return -1;
    ctx->havekernel = calloc(n_compares, sizeof(*ctx->havekernel));
    if (ctx->havekernel == NULL) return -1;
    ctx->num_compares = n_compares;
    if (prepfunc != NULL) {
        int ret = prepfunc(ctx, prepfunc_extra);
        if (ret != 0) {
            return ret;
        }
    }
    ctx->have_finalised = true;
    return 0;
}


int kerncalc_compute_kernel(kwip_kerncalc_t *ctx, size_t idx)
{
    if (ctx == NULL || ! ctx->have_finalised) return -1;

    int res;
    size_t row, col;

    // Get row and col into distance matrix from the comparison index.
    res = kernmatrix_condensed_to_ij(&row, &col, idx);
    if (res != 0) return res;

    if (row > ctx->num_samples || col > ctx->num_samples) return -1;

    const char *Afile = ctx->files[row];
    const char *Bfile = ctx->files[col];
    res = ctx->kernfunc(&ctx->kernelvalues[idx], Afile, Bfile, ctx->extra);
    if (res != 0) return res;
    ctx->havekernel[idx] = true;

    return 0;
}

int kerncalc_save(kwip_kerncalc_t *ctx)
{
    if (ctx == NULL || ctx->checkpoint_dir == NULL) return -1;

    char filename[4096] = "";
    int res = 0;

    res = snprintf(filename, 4096, "%s/kernellog.tab", ctx->checkpoint_dir);
    if (res >= 4096 || res < 0) return -1;

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) return -1;

    fprintf(fp, "# Kernel values generated with kWIP version %s\n", KWIP_VERSION);
    for (size_t i = 0; i < ctx->num_compares; i++) {
        size_t row, col;
        // Get row and col into distance matrix from the comparison index.
        res = kernmatrix_condensed_to_ij(&row, &col, i);
        if (res != 0) return -1;

        if (ctx->havekernel[i]) {
            fprintf(fp, "%s\t%s\t%zu\t%0.17g\n", ctx->samplenames[row],
                    ctx->samplenames[col], i, ctx->kernelvalues[i]);
        }
    }
    fclose(fp);
    return 0;
}
