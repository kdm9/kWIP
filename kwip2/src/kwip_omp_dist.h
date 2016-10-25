#ifndef KWIP_OMP_DIST_H_TPSMF8GN
#define KWIP_OMP_DIST_H_TPSMF8GN

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include <clogged/clogged.h>

#include "kwip_distcalc.h"

int distcalc_pairwise_omp(kwip_distcalc_t *ctx, clg_logger_t *log, int nthreads);

#endif /* end of include guard: KWIP_OMP_DIST_H_TPSMF8GN */
