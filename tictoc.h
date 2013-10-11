#pragma once
#include <time.h>

inline clock_t tic() { return clock(); }
inline double toc(clock_t t) { return (clock() - t) / ((double) CLOCKS_PER_SEC) ; }

