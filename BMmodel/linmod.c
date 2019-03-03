// defines dprocess, rprocess, dmeasure, rmeasure, mdprocess, mrprocess
#include "../stf_parallel.h"
#include "linmod.h"

using std::vector;
using std::cout;
using std::endl;

bool annealed = false; // the forecast likelihoods are NOT raised to fractional powers depending on the forecast time length

double ncdf(double value) { // normal cdf
  return 0.5 * std::erfc(-value * M_SQRT1_2);
}

double npdf(double value, double sd, bool logarithm) { // normal pdf
  return logarithm ?
    log(M_SQRT1_2 / (std::sqrt(M_PI) * sd)) - value * value / (2 * sd * sd) :
    M_SQRT1_2 / (std::sqrt(M_PI) * sd) * std::exp(-value * value / (2 * sd * sd));
}

void ltsolve(const vector<vector<double> >& L, vector<double>& x, const vector<double>& y) {
  // solve the linear equation Lx=y for lower triangular matrix L using forward substitution.
  int n = L.size();
  x = y;
  for (int i=0; i<n; i++) {
    for (int j=0; j<i; j++) 
      x[i] -= L[i][j] * x[j];
    x[i] /= L[i][i];
  }
}

// set initial state based on the initial value parameters
template <typename state_T, typename obs_T>
  void stf<state_T, obs_T>::set_initstates(vector<vector<vector<state_T> > >& initstates, const vector<vector<vector<state_T> > >& init_template, const vector<vector<vector<double> > >& thetaSwarm) {
  // for this example, we do not estimate the starting value.
  
  initstates = init_template;
}


// parameters: exp(theta[0]): s (process sd), exp(theta[1]): d (measurement sd), logit(theta[2]): alpha (off-diagonal correlation of multidimensional jumps)
// MODEL DESCRIPTION
// X_{t} = X_{t-1} + epsilon, where epsilon has covariance structure of compound symmetry (all off-diagonal entries are logit(theta[2]) times the diagonal entries)
// This process can be thought of as having a jump expressed as the sum of a common jump on all dimensions and an independent jump for each dimension.
// Y_{k,t} = X_{k,t} + N(0,1)
// imrprocess (continuous time rprocess) is defined such that the covariance of the jump is N times smaller.

// // // dprocess // // //
template<>
double stf<double, double>::dprocess(int t, const vector<double>& prev, const vector<double>& next, const vector<double>& theta, bool logarithm) {

  double dp = 0.0; // log density of process

  // TODO: compute the correct dprocess here (with correlated noise with compound symmetry)
  for (int k = 0; k < nvar; k++)
    dp += npdf(next[k] - prev[k], exp(theta[0]), true);
  
  return (logarithm ? dp : exp(dp));
}


// // // rprocess // // //
template<>
void stf<double, double>::rprocess(int t, const vector<double>& prev, vector<double>& next, const vector<double>& theta, std::mt19937& rng) {

  std::normal_distribution<double> normal(0.0, 1.0);

  next = prev;

  double jumpcommon = sqrt(1.0/(1.0+exp(-theta[2]))) * normal(rng);
  
  for (int k = 0; k < nvar; k++) {
    next[k] += exp(theta[0]) * (jumpcommon + sqrt(1.0-1.0/(1.0+exp(-theta[2]))) * normal(rng));
  } 
}


// // // dmeasure // // //
template<>
double stf<double, double>::dmeasure(int t, const vector<double>& state, const vector<double>& meas, const vector<double>& theta, bool logarithm) {

  double dm = 0.0; // log density of measurement

  for (int k=0; k<nvar; k++)
    dm += npdf(meas[k] - state[k], exp(theta[1]), true);

  return (logarithm ? dm : exp(dm));
}


// // // rmeasure // // //
template<>
void stf<double, double>::rmeasure(int t, const vector<double>& state, vector<double>& meas, const vector<double>& theta, std::mt19937& rng) {

  std::normal_distribution<double> normal(0.0, 1.0);

  for (int k = 0; k < nvar; k++)
    meas[k] = state[k] + exp(theta[1]) * normal(rng);
}




// // // intermediary processes for imfilter // // //
// // // imrprocess // // // 
template<>
void stf<double, double>::imrprocess(int t, int s, int im, const vector<double>& prev, vector<double>& next, const vector<double>& theta, std::mt19937& rng) {
  // the simulator of auxiliary process distribution of x_{t + (s+1)/im} given x_{t + s/im}

  std::normal_distribution<double> normal(0.0, 1.0);

  next = prev;

  double jumpcommon = sqrt(1.0/(1.0+exp(-theta[2]))) * normal(rng);

  for(int k=0; k<nvar; k++) 
    next[k] += ((jumpcommon + sqrt(1.0-1.0/(1.0+exp(-theta[2]))) * normal(rng)) * exp(theta[0]) / sqrt(im));

}


// // // imdmeasure // // //
template<>
vector<double> stf<double, double>::imdmeasure(int t, int s, vector<int> lookaheads, int im, const vector<double>& state, const vector<vector<double> >& meas, const vector<double>& theta, vector<vector<double> >& forecast_var, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+lookahead} given x_{t + s/im}
  // when lookahead = 1 and s = im, this function should be the same as the dmeasure function.
  // the length of lookaheads and the length of measurements should equal.

  vector<double> dm(lookaheads.size(), 0.0); // log density of measurement

  if (diagCov) { // when the guide function uses a diagonal approximation to the covariance matrix
    for (size_t pos = 0; pos < lookaheads.size(); pos++) {
      for (int k=0; k<nvar; k++)
	dm[pos] += npdf(meas[pos][k] - state[k], sqrt(exp(2*theta[1])+(lookaheads[pos]-(s+0.0)/im)*exp(2*theta[0])), true);
    }
  }
  if (!diagCov) { // when the guide function usese the exact covariance matrix
    for (size_t pos = 0; pos < lookaheads.size(); pos++) {
      int lookahead_steps = lookaheads[pos]*im-s; // the number of lookahead steps 
      vector<vector<double> > L = Lmat[lookahead_steps];
      vector<double> meas_minus_state = meas[pos];
      for (int k=0; k<nvar; k++) meas_minus_state[k] -= state[k];
      vector<double> z; // z = L^{-1}(meas-state)
      ltsolve(L,z,meas_minus_state);
      for (int k=0; k<nvar; k++)
	dm[pos] += (npdf(z[k], 1.0, true) - log(L[k][k])); // multivariate normal density
    }
  }
  
  vector<double> exp_dm(lookaheads.size());
  for (size_t pos=0; pos<lookaheads.size(); pos++)
    exp_dm[pos] = exp(dm[pos]);
  
  return (logarithm ? dm : exp_dm);
}


template<>
double stf<double, double>::imdmeasure(int t, int s, int im, const vector<double>& state, const vector<double>& meas, const vector<double>& theta, vector<double>& forecast_var_1step, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+1} given x_{t + s/im}
  // when s = im, this function should be the same as the dmeasure function.

  // set default values
  vector<int> lookaheads(1, 1);
  vector<vector<double> > observations(1, meas);
  vector<vector<double> > forecast_var(1, forecast_var_1step);

  vector<double> dm = imdmeasure(t, s, lookaheads, im, state, observations, theta, forecast_var, true);

  return (logarithm ? dm[0] : exp(dm[0]));
}

template class stf<double, double>;
