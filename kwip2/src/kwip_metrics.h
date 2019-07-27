#ifndef KWIP_METRICS_H_MMGE2ACT
#define KWIP_METRICS_H_MMGE2ACT

#include "kwip_config.h"
#include "kwip_array.h"


int metric_l1_dist(double *outp, const char *file1, const char *file2, void *extra);
int metric_l1_norm(double *outp, const char *file1, void *extra);

int metric_l2_dist(double *outp, const char *file1, const char *file2, void *extra);
int metric_l2_norm(double *outp, const char *file1, void *extra);


#endif /* end of include guard: KWIP_METRICS_H_MMGE2ACT */
