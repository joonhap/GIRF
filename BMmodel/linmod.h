#ifndef LINMOD_H
#define LINMOD_H

#include <vector>
#include <omp.h>

extern std::vector<std::vector<std::vector<double> > > Lmat; // cholesky factor (lower triangular matrix) of the covariance matrices for different lookahead interval lengths

#endif
