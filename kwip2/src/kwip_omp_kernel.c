#include "kwip_omp_kernel.h"

int
kerncalc_pairwise_omp(kwip_kerncalc_t *ctx, clg_logger_t *log, int nthreads)
{
    if (ctx == NULL || ! ctx->have_finalised) return -1;
    int res = 0;

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(ctx, res)
    for (size_t idx = 0; idx < ctx->num_compares; idx++) {
        int tl_res = kerncalc_compute_kernel(ctx, idx);

        if (tl_res != 0) {
            #pragma omp critical
            {
                res = tl_res;
                clg_log_fmt_error(log, "ERROR: comparison %zu returned %d\n", idx, tl_res);
            }
            continue;
        }

        if (res != 0) continue; // In error state

        #pragma omp critical
        {
            clg_log_fmt_progress(log, "\t- %zu done\n", idx);
        }
    }
    return res;
}
