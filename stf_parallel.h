#ifndef STF_PARALLEL_H
#define STF_PARALLEL_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <iterator>
#include <cmath>
#include <numeric>
#include <random>
#include <functional>

#define SEED_INC_PER_IF_ITER 98663 // the mersenne twister rng seed increases by this amount at each iteration of IF2

extern bool annealed; // whether the forecast likelihood estimates should be raised to fractional powers depending on the forecast length when evaluating the guide function. Can be defined in the code defining the model or in the code running experiments.

template <typename state_T, typename obs_T>
class stf {
 public:
  int nvar; // no. of variables (state dimension)
  int nvar_obs; // no. of variables (measurement dimension)
  int ntimes; // time length

  // obs stores ntimes-by-nvar elements in a vectorized form;
  // where row index changes faster than column index 
  std::vector<std::vector<obs_T> > obs; // observation
  std::vector<std::vector<state_T> > hidden_state; // hidden true state
  std::vector<double> theta; // true parameter
  std::vector<int> pos_ivp; // position of initial value parameters. All initial value parameters should come at the end of the theta vector.

  // simulator
  void simulate(std::vector<state_T>, int initial_seed);

  //// methods ////
  // constructor
  stf(int nv, int nv_obs, int nt) : nvar(nv), nvar_obs(nv_obs), ntimes(nt) { 
    // obs stores ntimes-by-nvar elements in a matrix format
    obs.resize(ntimes);
    for (int t=0; t<ntimes; t++) {
      obs[t].resize(nvar_obs);
    }
    // hidden_state stores (ntimes+1)-by-nvar elements in a matrix format
    hidden_state.resize(ntimes+1);
    for (int t=0; t<ntimes+1; t++) {
      hidden_state[t].resize(nvar);
    }
  }
  
  // destructor
  ~stf() { }
  
  // plain particle filter
  void ppfilter(const std::vector<std::vector<state_T> >& initstate, int J, std::vector<std::vector<double> > stateTrajectory, int initial_seed);

  // run intermediary filter
  void imfilter(const std::vector<std::vector<std::vector<state_T> > >& initstates, int im, int R, int J, std::vector<std::vector<double> >& lestimate, std::vector<std::vector<double> >& state_mean, std::vector<std::vector<state_T> >& q10, std::vector<std::vector<state_T> >& q50, std::vector<std::vector<state_T> >& q90, double (*f1)(double), std::vector<std::vector<double> >& f1_mean, std::vector<std::vector<std::vector<double> > >& ESS, std::vector<double> theta_point, bool weighted_interaction, int max_lookahead, std::function<bool (int, int)> island_interact_condition, int initial_seed); // filter without IF2
  void imfilter(const std::vector<std::vector<std::vector<state_T> > >& initstates, int im, int R, int J, std::vector<std::vector<double> >& lestimate, std::vector<std::vector<double> >& state_mean, std::vector<std::vector<state_T> >& q10, std::vector<std::vector<state_T> >& q50, std::vector<std::vector<state_T> >& q90, double (*f1)(double), std::vector<std::vector<double> >& f1_mean, std::vector<std::vector<std::vector<double> > >& ESS, const std::vector<std::vector<std::vector<double> > >& initThetaSwarm, std::vector<double> init_perturb_size, int M, std::function<double (int)> coolschedule, std::vector<std::vector<std::vector<double> > >& thetaTrajectory, std::vector<std::vector<std::vector<double> > >& finalThetaSwarm, bool weighted_interaction, int max_lookahead, std::function<bool (int, int)> island_interact_condition, int initial_seed); // IF2

  
  //// model ////
  
  double dprocess(int t, const std::vector<state_T>& prev, const std::vector<state_T>& next, const std::vector<double>& theta, bool logarithm);
  void rprocess(int t, const std::vector<state_T>& prev, std::vector<state_T>& next, const std::vector<double>& theta, std::mt19937& rng); 
  void rprocess(int t, const std::vector<state_T>& prev, std::vector<state_T>& next, const std::vector<double>& theta, double delta_t, std::mt19937& rng); // continuous time version
  double dmeasure(int t, const std::vector<state_T>& state, const std::vector<obs_T>& meas, const std::vector<double>& theta, bool logarithm);
  void rmeasure(int t, const std::vector<state_T>& state, std::vector<obs_T>& meas, const std::vector<double>& theta, std::mt19937& rng);

  //// models in imfilter (in addition to those already defined) ////
  void imrprocess(int t, int s, int im, const std::vector<state_T>& prev, std::vector<state_T>& next, const std::vector<double>& theta, std::mt19937& rng);
  void imrprocess(int t, int s, int im, const std::vector<state_T>& prev, std::vector<state_T>& next, const std::vector<double>& theta, double delta_t, std::mt19937& rng); // continuous time version
  std::vector<double> imdmeasure(int t, int s, std::vector<int> lookaheads, int im, const std::vector<state_T>& state, const std::vector<std::vector<obs_T> >& meas, const std::vector<double>& theta, bool logarithm); // imdmeasure with lookahead steps
  double imdmeasure(int t, int s, int im, const std::vector<state_T>& state, const std::vector<obs_T>& meas, const std::vector<double>& theta, bool logarithm); // default lookahead = 1
  std::vector<double> imdmeasure(int t, int s, std::vector<int> lookaheads, int im, const std::vector<state_T>& state, const std::vector<std::vector<obs_T> >& meas, const std::vector<double>& theta, std::vector<std::vector<double> >& forecast_var, bool logarithm); // imdmeasure with lookahead steps, using stored forecast variability estimates
  double imdmeasure(int t, int s, int im, const std::vector<state_T>& state, const std::vector<obs_T>& meas, const std::vector<double>& theta, std::vector<double>& forecast_var, bool logarithm); // default lookahead = 1, using stored forecast variability estimates

  //// set initial particle states using initial value parameter
  void set_initstates(std::vector<std::vector<std::vector<state_T> > >& initstates, const std::vector<std::vector<std::vector<state_T> > >& init_template, const std::vector<std::vector<std::vector<double> > >& thetaSwarm);

};



#include "stf_parallel.tpp" // code file containing key template class methods

#endif
