#ifndef LINMOD_H
#define LINMOD_H

#include <vector>
#include <omp.h>

extern std::vector<std::vector<std::vector<double> > > Lmat; // cholesky factor (lower triangular matrix) of the covariance matrices for different lookahead interval lengths
extern bool diagCov; // whether the guide function uses a diagonal approximation to the covariance matrix of the Brownian motion

#endif
