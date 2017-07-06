#include "kwip_distcalc.h"

#include "kwip_array.h"
#include "kwip_kmercount.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

int
distcalc_init(kwip_distcalc_t *ctx)
{
    if (ctx == NULL) return -1;
    memset(ctx, 0, sizeof(*ctx));
    return 0;
}

void
distcalc_destroy(kwip_distcalc_t *ctx)
{
    if (ctx == NULL) return;
    for (size_t i = 0; i < ctx->num_samples; i++) {
        kwip_free(ctx->files[i]);
        kwip_free(ctx->samplenames[i]);
    }
    kwip_free(ctx->files);
    kwip_free(ctx->samplenames);
    kwip_free(ctx->checkpoint_dir);
    kwip_free(ctx->distvalues);
    kwip_free(ctx->havedist);
}


char *
distcalc_filename_to_samplename(const char *filename) {
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
distcalc_add_sample(kwip_distcalc_t *ctx, const char *filename, const char *samplename)
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
        ctx->samplenames[ctx->num_samples] = distcalc_filename_to_samplename(filename);
    }
    if (ctx->samplenames[ctx->num_samples] == NULL) return -1;
    ctx->num_samples += 1;
    return 0;
}

int distcalc_set_metric(kwip_distcalc_t *ctx,
                        kwip_dist_fn_t distfunc,
                        kwip_norm_fn_t normfunc)
{
    if (ctx == NULL) return -1;

    ctx->distfunc = distfunc;
    ctx->normfunc = normfunc;
    return 0;
}

int distcalc_set_checkpoint_dir(kwip_distcalc_t *ctx, const char *dir)
{
    if (ctx == NULL) return -1;

    ctx->checkpoint_dir = strdup(dir);
    return 0;
}

// call finalise after adding samples, before calling distcalc_pairwise
int distcalc_finalise(kwip_distcalc_t *ctx, kwip_distcalc_finalise_fn_t prepfunc,
                      void *prepfunc_extra)
{
    if (ctx == NULL) return -1;

    if (ctx->distfunc == NULL) return -1;

    uint64_t n_samp = ctx->num_samples;
    if (n_samp < 2) return -1;

    uint64_t n_compares = n_samp * (n_samp - 1) / 2; // length of lower triang. matrix
    ctx->distvalues = calloc(n_compares, sizeof(*ctx->distvalues));
    if (ctx->distvalues == NULL) return -1;
    ctx->havedist = calloc(n_compares, sizeof(*ctx->havedist));
    if (ctx->havedist == NULL) return -1;
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


int distcalc_compute_dist(kwip_distcalc_t *ctx, size_t idx)
{
    if (ctx == NULL || ! ctx->have_finalised) return -1;

    int res;
    size_t row, col;

    // Get row and col into distance matrix from the comparison index.
    res = distmatrix_condensed_to_ij(&row, &col, idx);
    if (res != 0) return res;

    if (row > ctx->num_samples || col > ctx->num_samples) return -1;

    const char *Afile = ctx->files[row];
    const char *Bfile = ctx->files[col];
    res = ctx->distfunc(&ctx->distvalues[idx], Afile, Bfile, ctx->extra);
    if (res != 0) return res;
    ctx->havedist[idx] = true;

    return 0;
}

int distcalc_save(kwip_distcalc_t *ctx)
{
    if (ctx == NULL || ctx->checkpoint_dir == NULL) return -1;

    char filename[4096] = "";
    int res = 0;

    res = snprintf(filename, 4096, "%s/distlog.tab", ctx->checkpoint_dir);
    if (res >= 4096 || res < 0) return -1;

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) return -1;

    fprintf(fp, "# dist values generated with kWIP version %s\n", KWIP_VERSION);
    for (size_t i = 0; i < ctx->num_compares; i++) {
        size_t row, col;
        // Get row and col into distance matrix from the comparison index.
        res = distmatrix_condensed_to_ij(&row, &col, i);
        if (res != 0) return -1;

        if (ctx->havedist[i]) {
            fprintf(fp, "%s\t%s\t%zu\t%0.17g\n", ctx->samplenames[row],
                    ctx->samplenames[col], i, ctx->distvalues[i]);
        }
    }
    fclose(fp);
    return 0;
}
