#ifndef KWIP_UTILS_H_GNLDS27X
#define KWIP_UTILS_H_GNLDS27X

#include "kwip_config.h"
#include <stdlib.h>
#include <stdio.h>

size_t kwip_parse_size(const char *sizestr);

void kwip_print_banner(FILE *stream);
void kwip_print_version(FILE *stream);

int kwip_mkdirp(const char *dir);

#define kwip_free(ptr) do {if ((ptr) != NULL) { free(ptr); ptr = NULL; }} while(0)


#endif /* end of include guard: KWIP_UTILS_H_GNLDS27X */
