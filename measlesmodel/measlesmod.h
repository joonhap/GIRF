#ifndef MEASLESMOD_H
#define MEASLESMOD_H

#include <omp.h>

#define NCITY 40
#define NYEAR 21

extern double dist[NCITY][NCITY];
extern int birth_data[NYEAR][NCITY];
extern int birth_data_first_year;
extern int birth_data_last_year;
extern int pop_data[NYEAR][NCITY];
extern int pop_data_first_year;
extern int pop_data_last_year;
extern int cityindex[NCITY];
extern double dist_mean;
extern double pop_mean;
extern double school_prop;

void simulate_rprocess(double t_start, double t_end, const vector<int>& prev, vector<int>& next, const vector<double>& theta, double delta_t, std::mt19937& rng);
void simulate_deterministic_process(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t);

#pragma omp threadprivate(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)

#endif
