#ifndef KWIP_OMP_KERNEL_H_S1EBYNQF
#define KWIP_OMP_KERNEL_H_S1EBYNQF

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "clogged/clogged.h"

#include "kwip_kernelcalc.h"


int kerncalc_pairwise_omp(kwip_kerncalc_t *ctx, clg_logger_t *log, int nthreads);

#endif /* end of include guard: KWIP_OMP_KERNEL_H_S1EBYNQF */
