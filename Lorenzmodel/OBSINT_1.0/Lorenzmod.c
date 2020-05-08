#include <cfloat>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include "../../stf_parallel.h"
#include "Lorenzmod.h"

using std::vector;
using std::cout;
using std::endl;

// this version (measlesmod_momentmatching.c) computes the guide function by adding the forecast variability to the variance of the measurement process.
// The forecast variability is estimated by simulating forward the rprocess multiple times and computing the sample statistics (e.g., inter-quartile range or variance).

bool annealed = true; // forecast likelihood estimates are raised to fractional powers depending on the forecast time length

double npdf(double value, double mean, double sd, bool logarithm) { // normal pdf
  double l = -(value-mean)*(value-mean)/(2*sd*sd)-log(2*M_PI)/2-log(sd); // log pdf
  return (logarithm ? l : exp(l));
}


// set initial state based on the initial value parameters
template <typename state_T, typename obs_T>
  void stf<state_T, obs_T>::set_initstates(vector<vector<vector<state_T> > >& initstates, const vector<vector<vector<state_T> > >& init_template, const vector<vector<vector<double> > >& thetaSwarm) {
  // the initial value parameter vector defines the proportion of each compartment at t=0
  // as a precaution, check if the length of the first dimension of init_tempate and thetaSwarm match.
  if (init_template.size() != thetaSwarm.size()) {
    cout << "Error in initializing particle states: the length of init_template and that of thetaSwarm do not match." << endl;
    exit(1);
  }
  initstates = init_template;
}


// // // deterministic skeleton // // //
void simulate_deterministic_process(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t) {
  // process simulator for either deterministic skeleton
  // dx_i/dt = (x_{i+1} - x_{i-2}) * x_{i-1} - x_i + F

  const double F = theta[0]; // external forcing
  const int dim = prev.size();

  vector<double> pv = prev; // spaceholder for the previous state during incremental updates
  vector<double> nx = pv; // spaceholder for the next sate during incremental updates

  const double tol=1e-5;

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end
    }

    nx[0] = pv[0] + t_inc * ( (pv[1] - pv[dim-2]) * pv[dim-1] - pv[0] + F);
    nx[1] = pv[1] + t_inc * ( (pv[2] - pv[dim-1]) * pv[0] - pv[1] + F);
    for (int k=2; k<(dim-1); k++) {
      nx[k] = pv[k] + t_inc * ( (pv[k+1] - pv[k-2]) * pv[k-1] - pv[k] + F);
    }
    nx[dim-1] = pv[dim-1] + t_inc * ( (pv[0] - pv[dim-3]) * pv[dim-2] - pv[dim-1] + F);
    pv = nx;
  }
  next = nx;
}


// // // rprocess // // //
void simulate_rprocess(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t, std::mt19937& rng) {
  // process simulator for rprocess
  // dx_i = [ (x_{i+1} - x_{i-2}) * x_{i-1} - x_i + F ] * dt + sigma_p * dB_t

  const double F = theta[0]; // external forcing
  const double sig_p = exp(theta[1]); // process noise standard deviation
  const int dim = prev.size();

  std::normal_distribution<double> normal(0.0, 1.0);

  vector<double> pv = prev; // spaceholder for the previous state during incremental updates
  vector<double> nx = pv; // spaceholder for the next sate during incremental updates

  const double tol=1e-5;

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end
    }

    nx[0] = pv[0] + t_inc * ( (pv[1] - pv[dim-2]) * pv[dim-1] - pv[0] + F) + sqrt(t_inc) * sig_p * normal(rng);
    nx[1] = pv[1] + t_inc * ( (pv[2] - pv[dim-1]) * pv[0] - pv[1] + F) + sqrt(t_inc) * sig_p * normal(rng);
    for (int k=2; k<(dim-1); k++) {
      nx[k] = pv[k] + t_inc * ( (pv[k+1] - pv[k-2]) * pv[k-1] - pv[k] + F) + sqrt(t_inc) * sig_p * normal(rng);
    }
    nx[dim-1] = pv[dim-1] + t_inc * ( (pv[0] - pv[dim-3]) * pv[dim-2] - pv[dim-1] + F) + sqrt(t_inc) * sig_p * normal(rng);
    pv = nx;
  }
  next = nx;
}


template<>
void stf<double, double>::rprocess(int n, const vector<double>& prev, vector<double>& next, const vector<double>& theta, std::mt19937& rng) {
  // n: observation time index (related to time t via the equation t = OBS_INT * n)

  simulate_rprocess(n*OBS_INT, (n+1)*OBS_INT, prev, next, theta, 0.01, rng);
}


// // // rmeasure // // //
template<>
void stf<double, double>::rmeasure(int n, const std::vector<double>& state, std::vector<double>& meas, const std::vector<double>& theta, std::mt19937& rng) {
  const double sig_m = exp(theta[2]); // measurement noise standard deviation
  const int dim = state.size();
  
  std::normal_distribution<double> normal(0.0, 1.0);

  for (int k = 0; k < dim; k++) {
    meas[k] = state[k] + sig_m * normal(rng);
  }
}


// // // dmeasure // // //
template<>
double stf<double, double>::dmeasure(int n, const std::vector<double>& state, const std::vector<double>& meas, const std::vector<double>& theta, bool logarithm) {
  const double sig_m = exp(theta[2]); // measurement noise standard deviation
  const int dim = state.size();

  double dm = 0.0; // log density of measurement
  
  for (int k = 0; k < dim; k++) {
    dm += npdf(meas[k], state[k], sig_m, true);
  }

  return (logarithm ? dm : exp(dm));
}



// // // intermediary processes for imfilter // // //
// // // imrprocess // // // 
template<>
void stf<double, double>::imrprocess(int t, int s, int im, const vector<double>& prev, vector<double>& next, const vector<double>& theta, std::mt19937& rng) {
  // the simulator of process distribution of x_{t + (s+1)/im} given x_{t + s/im} using a default time increment of 0.01.

  simulate_rprocess((t+(s+0.0)/im)*OBS_INT, (t+(s+1.0)/im)*OBS_INT, prev, next, theta, 0.01, rng);
}



// // // imdmeasure // // //
template<>
vector<double> stf<double, double>::imdmeasure(int n, int s, vector<int> lookaheads, int im, const vector<double>& state, const vector<vector<double> >& meas, const vector<double>& theta, vector<vector<double> >& forecast_var, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{n+lookahead} given x_{(n+s/im)*OBS_INT}
  // note that the input s will take values between 1 and im (does not take value zero). Check the code stf_parallel.tpp
  // we use projection of state to n+lookahead to compute pseudo-predictive likelihood
  // when lookahead = 1 and s = im, this function should be the same as the dmeasure function.
  // the length of lookaheads and the length of meas should equal.

  const double sig_m = exp(theta[2]); // measurement noise standard deviation
  const int dim = state.size();

  // measurement density (log)
  vector<double> dm(lookaheads.size(), 0.0);
  vector<double> exp_dm(lookaheads.size(), 0.0);

  int lasize = lookaheads.size(); // the total number of lookaheads

  bool estimate_forecast_var = (s==1); // whether we update the estimated forecast variability at the current time point. Should always be re-estimated at s=1.
  int s_last_fve = 1; // last intermediate step at which the forecast variability was estimated. Should correspond with how estimate_forecast_var was set.
  if (estimate_forecast_var) {
    int nforecast = 40; // number of random forecasts to make
    vector<vector<vector<double> > > forecasts(nforecast, vector<vector<double> >(lasize)); // array structure: [forecast_id][lookahead_entry_id][space_dim]
    std::mt19937 rng(4255427);
    double delta_t = 0.01;
    for (int fc = 0; fc < nforecast; fc++) {
      simulate_rprocess((n+(s+0.0)/im)*OBS_INT, (n+lookaheads[0])*OBS_INT, state, forecasts[fc][0], theta, delta_t, rng);
      for (int pos = 1; pos < lasize; pos++) {
	simulate_rprocess((n+lookaheads[pos-1])*OBS_INT, (n+lookaheads[pos])*OBS_INT, forecasts[fc][pos-1], forecasts[fc][pos], theta, delta_t, rng);
      }
    }

    for (int pos = 0; pos < lasize; pos++) {
      forecast_var[pos].resize(dim);
      for (int k = 0; k < dim; k++) {
	vector<double> forecasts_rar(nforecast); // rearranged forecasts
	for (int fc = 0; fc < nforecast; fc++) {
	  forecasts_rar[fc] = forecasts[fc][pos][k];
	}
	quicksort(forecasts_rar, 0, nforecast-1);
	double iqr = forecasts_rar[int(nforecast*.75)] - forecasts_rar[int(nforecast*.25)]; // inter-quartile range
	double var_to_iqr2 = .55 * (1+2/sqrt(nforecast)); // distribution-specific ratio between the variance and the squared iqr (for scale families) (~0.55 for normal). Take a conservative number to account for the variability in the estimates of forecast variability itself. This translates into inflating the number by a factor (1+2*sqrt(nforecast)).
	forecast_var[pos][k] = iqr*iqr*var_to_iqr2;
      }
    }
  }

  vector<double> forc_base = state; // initial state for making a deterministic forecast
  vector<double> det_forc; // deterministic forecast

  for (int pos = 0; pos < lasize; pos++) {
    double delta_t = 0.01;
    simulate_deterministic_process((pos == 0 ? (n+(s+0.0)/im)*OBS_INT : (n+lookaheads[pos-1])*OBS_INT), (n+lookaheads[pos])*OBS_INT, forc_base, det_forc, theta, delta_t);
    forc_base = det_forc;
    double scale_fv = (lookaheads[pos] - (s+0.0)/im) / std::max((lookaheads[pos] - (s_last_fve+0.0)/im), 1e-6); // scaling factor for forecast variance
    for (int k = 0; k < dim; k++) {
      dm[pos] += npdf(meas[pos][k], det_forc[k], std::sqrt(sig_m*sig_m + scale_fv * forecast_var[pos][k]), true);
    }
  }
  if (!logarithm) {
    for (int pos = 0; pos < lasize; pos++)
      exp_dm[pos] = exp(dm[pos]);
  }
  return (logarithm ? dm : exp_dm);
}



template<>
double stf<double, double>::imdmeasure(int n, int s, int im, const vector<double>& state, const vector<double>& meas, const vector<double>& theta, vector<double>& forecast_var_1step, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+1} given x_{t+s/im} (ONLY to the immediate next time point)
  // when s = im, this function should be the same as the dmeasure function.

  // set default values
  vector<int> lookaheads(1, 1);
  vector<vector<double> > observations(1, meas);
  vector<vector<double> > forecast_var(1, forecast_var_1step);

  vector<double> dm = imdmeasure(n, s, lookaheads, im, state, observations, theta, forecast_var, logarithm);

  return (logarithm ? dm[0] : exp(dm[0]));
}


template class stf<double, double>;
