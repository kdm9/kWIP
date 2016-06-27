#include "greatest.h"

#include "tests/array.c"


GREATEST_MAIN_DEFS();

int main(int argc, char *argv[])
{
    GREATEST_MAIN_BEGIN();
    RUN_SUITE(array);
    GREATEST_MAIN_END();
    return 0;
}
