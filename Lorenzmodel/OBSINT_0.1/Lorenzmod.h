#ifndef LORENZMOD_H
#define LORENZMOD_H

#include <omp.h>

#define OBS_INT 0.1 // observation time interval

void simulate_rprocess(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t, std::mt19937& rng);
void simulate_deterministic_process(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t);

#endif
