#include <cfloat>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include "../stf_parallel.h"
#include "measlesmod.h"

using std::vector;
using std::cout;
using std::endl;

// this version (measlesmod_variational.c) uses variational-inference-like method to compute the guide function.
// Specifically, this version estimates the marginal likelihoods of data corresponding to one city using SMC using only the data for that one city.
// The guide function is taken to be the product of marginal likelihood estimates.

double dnpmf(int value, double mean, double sd) { // discrete normal pmf
  if (value > 0.5) {
    if (value < mean)
      return 0.5 * (std::erfc(-(value + 0.5 - mean) / sd * M_SQRT1_2) - std::erfc(-(value - 0.5 - mean) / sd * M_SQRT1_2));
    else
      return -0.5 * (std::erfc((value + 0.5 - mean) / sd * M_SQRT1_2) - std::erfc((value - 0.5 - mean) / sd * M_SQRT1_2));
  }
  else
    return 0.5 * std::erfc(-(0 + 0.5 - mean) / sd * M_SQRT1_2);
}

double dnpmf(double value, double mean, double sd) { // discrete normal pmf where value can take real values
  if (value > 0.5) {
    if (value < mean)
      return 0.5 * (std::erfc(-(value + 0.5 - mean) / sd * M_SQRT1_2) - std::erfc(-(value - 0.5 - mean) / sd * M_SQRT1_2));
    else
      return -0.5 * (std::erfc((value + 0.5 - mean) / sd * M_SQRT1_2) - std::erfc((value - 0.5 - mean) / sd * M_SQRT1_2));
  }
  else
    return 0.5 * std::erfc(-(0 + 0.5 - mean) / sd * M_SQRT1_2);
}

double npdf(double value, double mean, double sd, bool logarithm) { // normal pdf
  double l = -(value-mean)*(value-mean)/(2*sd*sd)-log(2*M_PI)/2-log(sd); // log pdf
  return (logarithm ? l : exp(l));
}

double negbinpmf(int value, double mean, double var, bool logarithm) { // negative binomial distribution with given mean and variance
  if (value < 0) { return (logarithm ? log(0.0) : 0.0); }
  if (mean < 1e-6) { double pmf = (value==0 ? 1.0 : 0.0); return (logarithm ? log(pmf) : pmf); } // the case when mean = 0.0;
  if (mean >= var) { double l = -lgamma(value+1) - mean + value * log(mean); return (logarithm ? l : exp(l)); } // if mean >= var, use Poisson pmf
  double p = 1.0-mean/var; // success probability
  double r = mean*(1.0-p)/p;
  double l = lgamma(value+r)-lgamma(value+1)-lgamma(r)+r*log(1-p)+value*log(p); // log pmf
  return (logarithm ? l : exp(l));
}

double hol_eff(int d, double hol_coeff) { // compute holiday effect coefficient
  
  if (d < 7) return hol_coeff; // New year
  if (99 < d && d < 116) return hol_coeff; // Easter
  if (198 < d && d < 253) return hol_coeff; // summer
  if (299 < d && d < 309) return hol_coeff; // autumn
  if (355 < d) return hol_coeff; // Chirstmas

  return 1.0;
}


// set initial state based on the initial value parameters
template <typename state_T, typename obs_T>
  void stf<state_T, obs_T>::set_initstates(vector<vector<vector<state_T> > >& initstates, const vector<vector<vector<state_T> > >& init_template, const vector<vector<vector<double> > >& thetaSwarm) {
  // the initial value parameter vector (length 1) defines the proportion of the susceptible population at t=0
  // the length of the first dimension of init_tempate and thetaSwarm should match.
  if (init_template.size() != thetaSwarm.size()) {
    cout << "Error in initializing particle states: the length of init_template and that of thetaSwarm do not match." << endl;
    exit(1);
  }
  initstates = init_template;
}


// // // deterministic skeleton // // //
void simulate_deterministic_process(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t) {
  // process simulator for either deterministic skeleton or rprocess (as specified by the parameter 'deterministic')
  // prev, next: state vector. Length 5*k. Consists of k copies of (S,E,I,I_to_R_transition,P). P (population) is (piecewise) constant.
  // I_to_R_transition is the aggregate transitions from I to R for each observation time interval.
  // t: week of a year, delta_t: euler increment step. This rprocess operates between time t and t+1.

  const int ncity = prev.size() / 5; // number of cities being analyzed
  const double R0 = exp(theta[0]); // basic reproduction number
  const double amp = 1.0/(1.0+exp(-theta[1]))*.68; // seasonal amplitude
  const double gen = exp(theta[4]); // generation cycle (1/nu_EI+1/nu_IR, log-transformed)
  const double infec_period_prop = 1.0/(1.0+exp(-theta[5])); // infectious period proportion
  const double nu_IR = 1.0/(gen*infec_period_prop); // 1/infectious period (log-transformed)
  const double nu_EI = 1.0/(gen*(1.0-infec_period_prop)); // 1/latent period
  const double beta_s = (1.0+2*(1.0-school_prop)*amp)*R0*nu_IR; // school time transmission coefficient
  const double hol_coeff = (1-2*school_prop*amp) / (1+2*(1-school_prop)*amp); // holiday effect on transmission coefficient (<1)
  const double alpha = theta[2]; // mixing exponent
  const double mu = theta[3]; // mortality rate
  // const double sigma2 = exp(theta[6]); // variability constant (for gamma process, sigma^2. see Box 1 in He et al 2009) (unused for deterministic process)
  const double G = exp(theta[9]); // gravity model coefficient
  const double source_pow = theta[10]; // exponent to source population size in gravity model
  const double dest_pow = theta[11]; // exponent to destination population size in gravity model
  const double dist_pow = theta[12]; // exponent to distance in gravity model
  const double cohort_entry = 1.0/(1.0+exp(-theta[13])); // cohort entry fraction
  const int schoolstart = int(theta[14]+.5); // school start day (given as double in theta)
  const int entryage = int(theta[15]+.5); // school entry age (given as double in theta)
  const int startyear = int(theta[16]+.5); // start year (the year at which t=0)

  vector<double> pv(5*ncity); // for deterministic skeleton, the state needs to take real values
  for (int n=0; n<5*ncity; n++)
    pv[n] = prev[n];
  vector<double> nx = pv;

  const double tol=1e-6;
  const double oti=14.0; // observation time interval in unit of days (14.0 means biweekly data)

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end
    }

    // initialize total I to R transitions per observation time interval (if an integer is in the interval [time, time+t_inc), set the cumulative sum to zero.)
    if (int(time-tol) < int(time+t_inc-tol)) {
      for (int k=0; k < ncity; k++)
	pv[5*k+3] = 0; 
    }

    int year = int(time * oti / 365.25) + startyear; 
    int day = int(time * oti + tol - (year-startyear)*365.25);

    // change population size at the start of year
    if (year < int((time+t_inc) * oti / 365.25) + startyear) {
      // if population data is not found for the new year, exit the program
      if ((int((time+t_inc) * oti / 365.25) + startyear < pop_data_first_year) || (int((time+t_inc) * oti / 365.25) + startyear > pop_data_last_year)) {
	cout << "Population data not found for the new year (Population data needed for year " << (int((time+t_inc) * oti / 365.25) + startyear) << ".) Exit the program." << endl;
	exit(1);
      }
      for (int k=0; k<ncity; k++)
	pv[5*k+4] = pop_data[int((time+t_inc) * oti / 365.25) + startyear-pop_data_first_year][cityindex[k]];
    }
    
    // recruit incremental susceptibles
    double day_double = time*oti - int(time*oti/365.25) * 365.25; // days since the start of the year
    for (int k=0; k<ncity; k++) 
      pv[5*k] += (int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * (day_double+t_inc*oti) / 365.25) - int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * day_double / 365.25));

    // recruit susceptibles on school start day
    if((day < schoolstart) & (schoolstart <= int((time+t_inc)*oti + tol - (year-startyear)*365.25))) { // make sure not to recruit susceptibles multiples time a day
      // if birth data is not found for the cohort, exit the program
      if ((year-entryage < birth_data_first_year) || (year-entryage > birth_data_last_year)) {
	cout << "Birth data not found for the current year (Birth data needed for year " << (year - entryage) << ".) Exit the program." << endl;
	exit(1);
      }

      for (int k=0; k<ncity; k++) {
	pv[5*k] += int(cohort_entry * birth_data[year-entryage-birth_data_first_year][cityindex[k]]);
      }
    }

    for(int k=0; k<ncity; k++) {

      double nu_SE;
      nu_SE = pow(double(pv[5*k+2])/double(pv[5*k+4]), alpha); // per capita exposure proportion

      vector<double> ntraveler(ncity, 0.0);
      for(int l=0; l<ncity; l++) {
	if (k==l)
	  continue;
	ntraveler[l] = G * pow(double(pv[5*k+4])/pop_mean, source_pow) * pow(double(pv[5*l+4])/pop_mean, dest_pow) / pow(dist[k][l]/dist_mean, dist_pow);
      }
      double tot_traveler = std::accumulate(ntraveler.begin(), ntraveler.end(), 0.0);
      double travel_prob = 1.0 - exp(-tot_traveler/pv[5*k+4]); // in case tot_traveler exceeds the city's population (in IF2), transform the value.
      if (tot_traveler > tol) {
	for(int l=0; l<ncity; l++) {
	  if (k==l)
	    continue;
	  nu_SE += ( travel_prob * ntraveler[l] / tot_traveler * (pow(double(pv[5*l+2])/pv[5*l+4], alpha) - pow(double(pv[5*k+2])/pv[5*k+4], alpha)) );
	}
      }
      nu_SE *= (beta_s * hol_eff(day, hol_coeff));
      
      double gamma_SE = t_inc * oti / 7.0;
      double gamma_SD = t_inc * oti / 7.0;
      double gamma_EI = t_inc * oti / 7.0;
      double gamma_ED = t_inc * oti / 7.0;
      double gamma_IR = t_inc * oti / 7.0;
      double gamma_ID = t_inc * oti / 7.0;
      
      double outrate_S = nu_SE*gamma_SE + mu*gamma_SD;
      double p_SE = nu_SE*gamma_SE / outrate_S * (1-std::exp(-outrate_S));
      double p_SD = mu*gamma_SD / outrate_S * (1-std::exp(-outrate_S));
      double dN_SE = pv[5*k] * p_SE; 
      double dN_SD = (pv[5*k]-dN_SE) * (p_SD/(1-p_SE)); 

      double outrate_E = nu_EI*gamma_EI + mu*gamma_ED;
      double p_EI = nu_EI*gamma_EI / outrate_E * (1-std::exp(-outrate_E));
      double p_ED = mu*gamma_ED / outrate_E * (1-std::exp(-outrate_E));
      double dN_EI = pv[5*k+1] * p_EI; 
      double dN_ED = (pv[5*k+1]-dN_EI) * (p_ED/(1-p_EI)); 
      
      double outrate_I = nu_IR*gamma_IR + mu*gamma_ID;
      double p_IR = nu_IR*gamma_IR / outrate_I * (1-std::exp(-outrate_I));
      double p_ID = mu*gamma_ID / outrate_I * (1-std::exp(-outrate_I));
      double dN_IR = pv[5*k+2] * p_IR; 
      double dN_ID = (pv[5*k+2]-dN_IR) * (p_ID/(1-p_IR));

      nx[5*k] = pv[5*k] -dN_SE - dN_SD;
      nx[5*k+1] = pv[5*k+1] + dN_SE - dN_EI - dN_ED;
      nx[5*k+2] = pv[5*k+2] + dN_EI - dN_IR - dN_ID;
      nx[5*k+3] = pv[5*k+3] + dN_IR;
      nx[5*k+4] = pv[5*k+4];
    }
    pv = nx;
  }
  next = nx;
}


// // // rprocess // // //

void simulate_rprocess(double t_start, double t_end, const vector<int>& prev, vector<int>& next, const vector<double>& theta, double delta_t, std::mt19937& rng) {
  // process simulator for either deterministic skeleton or rprocess (as specified by the parameter 'deterministic')
  // prev, next: state vector. Length 5*k. Consists of k copies of (S,E,I,I_to_R_transition,P). P (population) is (piecewise) constant.
  // I_to_R_transition is the aggregate transitions from I to R for each observation time interval.
  // t: week of a year, delta_t: euler increment step. This rprocess operates between time t and t+1.

  const int ncity = prev.size() / 5; // number of cities being analyzed
  const double R0 = exp(theta[0]); // basic reproduction number
  const double amp = 1.0/(1.0+exp(-theta[1]))*.68; // seasonal amplitude
  const double gen = exp(theta[4]); // generation cycle (1/nu_EI+1/nu_IR, log-transformed)
  const double infec_period_prop = 1.0/(1.0+exp(-theta[5])); // infectious period proportion
  const double nu_IR = 1.0/(gen*infec_period_prop); // 1/infectious period (log-transformed)
  const double nu_EI = 1.0/(gen*(1.0-infec_period_prop)); // 1/latent period
  const double beta_s = (1.0+2*(1.0-school_prop)*amp)*R0*nu_IR; // school time transmission coefficient
  const double hol_coeff = (1-2*school_prop*amp) / (1+2*(1-school_prop)*amp); // holiday effect on transmission coefficient (<1)
  const double alpha = theta[2]; // mixing exponent
  const double mu = theta[3]; // mortality rate
  const double sigma2 = exp(theta[6]); // variability constant (for gamma process, sigma^2. see Box 1 in He et al 2009) (log-transformed)
  const double G = exp(theta[9]); // gravity model coefficient (log-transformed)
  const double source_pow = theta[10]; // exponent to source population size in gravity model
  const double dest_pow = theta[11]; // exponent to destination population size in gravity model
  const double dist_pow = theta[12]; // exponent to distance in gravity model
  const double cohort_entry = 1.0/(1.0+exp(-theta[13])); // cohort entry fraction
  const int schoolstart = int(theta[14]+.5); // school start day (given as double in theta)
  const int entryage = int(theta[15]+.5); // school entry age (given as double in theta)
  const int startyear = int(theta[16]+.5); // start year (the year at which t=0)

  vector<int> pv = prev; // starting state of the infinitesimal interval
  vector<int> nx = pv; // end state of the infinitesimal interval

  const double tol=1e-6;
  const double oti=14.0; // observation time interval in unit of days (14.0 means biweekly data)

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end.
    }

    // initialize total I to R transitions per observation time interval (if an integer is in the interval [time, time+t_inc), set the cumulative sum to zero.)
    if (int(time-tol) < int(time+t_inc-tol)) {
      for (int k=0; k < ncity; k++)
	pv[5*k+3] = 0; 
    }

    int year = int(time * oti / 365.25) + startyear; 
    int day = int(time * oti + tol - (year-startyear)*365.25);

    // change population size at the start of year
    if (year < int((time+t_inc) * oti / 365.25) + startyear) {
      // if population data is not found for the new year, exit the program
      if ((int((time+t_inc) * oti / 365.25) + startyear < pop_data_first_year) || (int((time+t_inc) * oti / 365.25) + startyear > pop_data_last_year)) {
	cout << "Population data not found for the new year (Population data needed for year " << (int((time+t_inc) * oti / 365.25) + startyear) << ".) Exit the program." << endl;
	exit(1);
      }
      for (int k=0; k<ncity; k++)
	pv[5*k+4] = pop_data[int((time+t_inc) * oti / 365.25) + startyear-pop_data_first_year][cityindex[k]]; 
    }

    // recruit incremental susceptibles
    double day_double = time*oti - int(time*oti/365.25) * 365.25; // days since the start of the year

    for (int k=0; k<ncity; k++) {
      pv[5*k] += (int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * (day_double+t_inc*oti) / 365.25) - int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * day_double / 365.25));
    }

    // recruit susceptibles on school start day
    if((day < schoolstart) & (schoolstart <= int((time+t_inc)*oti + tol - (year-startyear)*365.25))) { // make sure not to recruit susceptibles multiples time a day
      // if birth data is not found for the cohort, exit the program
      if ((year-entryage < birth_data_first_year) || (year-entryage > birth_data_last_year)) {
	cout << "Birth data not found for the current year (Birth data needed for year " << (year - entryage) << ".) Exit the program." << endl;
	exit(1);
      }

      for (int k=0; k<ncity; k++) {
	pv[5*k] += int(cohort_entry * birth_data[year-entryage-birth_data_first_year][cityindex[k]]);
      }
    }

    for(int k=0; k<ncity; k++) {

      double nu_SE;
      nu_SE = pow(double(pv[5*k+2])/double(pv[5*k+4]), alpha); // per capita exposure proportion

      vector<double> ntraveler(ncity, 0.0);
      for(int l=0; l<ncity; l++) {
	if (k==l)
	  continue;
	ntraveler[l] = G * pow(double(pv[5*k+4])/pop_mean, source_pow) * pow(double(pv[5*l+4])/pop_mean, dest_pow) / pow(dist[k][l]/dist_mean, dist_pow);
      }
      double tot_traveler = std::accumulate(ntraveler.begin(), ntraveler.end(), 0.0);
      double travel_prob = 1.0 - exp(-tot_traveler/pv[5*k+4]); // in case tot_traveler exceeds the city's population (in IF2), transform the value.
      if (tot_traveler > tol) {
	for(int l=0; l<ncity; l++) {
	  if (k==l)
	    continue;
	  nu_SE += ( travel_prob * ntraveler[l] / tot_traveler * (pow(double(pv[5*l+2])/pv[5*l+4], alpha) - pow(double(pv[5*k+2])/pv[5*k+4], alpha)) );
	}
      }
      nu_SE *= (beta_s * hol_eff(day, hol_coeff));
      
      std::gamma_distribution<double> gamma(t_inc*oti/7.0/sigma2, sigma2);
      double gamma_SE = gamma(rng); // gamma variable for S to E transition
      double gamma_SD = gamma(rng); // gamma variable for S to D (death) transition
      double gamma_EI = gamma(rng); // gamma variable for E to I transition
      double gamma_ED = gamma(rng); // gamma variable for E to D transition
      double gamma_IR = gamma(rng); // gamma variable for I to R transition
      double gamma_ID = gamma(rng); // gamma variable for I to D transition

      double outrate_S = nu_SE*gamma_SE + mu*gamma_SD;
      double p_SE = nu_SE*gamma_SE / outrate_S * (1-std::exp(-outrate_S));
      double p_SD = mu*gamma_SD / outrate_S * (1-std::exp(-outrate_S));
      int dN_SE;
      { std::binomial_distribution<int> binom(pv[5*k], p_SE); dN_SE = binom(rng); }
      int dN_SD;
      { std::binomial_distribution<int> binom(pv[5*k]-dN_SE, p_SD/(1-p_SE)); dN_SD = binom(rng); }

      double outrate_E = nu_EI*gamma_EI + mu*gamma_ED;
      double p_EI = nu_EI*gamma_EI / outrate_E * (1-std::exp(-outrate_E));
      double p_ED = mu*gamma_ED / outrate_E * (1-std::exp(-outrate_E));
      int dN_EI;
      { std::binomial_distribution<int> binom(pv[5*k+1], p_EI); dN_EI = binom(rng); }
      int dN_ED;
      { std::binomial_distribution<int> binom(pv[5*k+1]-dN_EI, p_ED/(1-p_EI)); dN_ED = binom(rng); }

      double outrate_I = nu_IR*gamma_IR + mu*gamma_ID;
      double p_IR = nu_IR*gamma_IR / outrate_I * (1-std::exp(-outrate_I));
      double p_ID = mu*gamma_ID / outrate_I * (1-std::exp(-outrate_I));
      int dN_IR;
      { std::binomial_distribution<int> binom(pv[5*k+2], p_IR); dN_IR = binom(rng); }
      int dN_ID;
      { std::binomial_distribution<int> binom(pv[5*k+2]-dN_IR, p_ID/(1-p_IR)); dN_ID = binom(rng); }

      nx[5*k] = pv[5*k] - dN_SE - dN_SD;
      nx[5*k+1] = pv[5*k+1] + dN_SE - dN_EI - dN_ED;
      nx[5*k+2] = pv[5*k+2] + dN_EI - dN_IR - dN_ID;
      nx[5*k+3] = pv[5*k+3] + dN_IR;
      nx[5*k+4] = pv[5*k+4];
      
    }
    pv = nx;
  }

  next = nx;
}


template<>
void stf<int, int>::rprocess(int t, const vector<int>& prev, vector<int>& next, const vector<double>& theta, double delta_t, std::mt19937& rng) {
  // t: week of a year, delta_t: euler increment step. This rprocess operates between time t and t+1.

  simulate_rprocess(t+0.0, t+1.0, prev, next, theta, delta_t, rng);
}


template<>
void stf<int, int>::rprocess(int t, const vector<int>& prev, vector<int>& next, const vector<double>& theta, std::mt19937& rng) {
  // t: week of a year, delta_t: euler increment step. This rprocess operates between time t and t+1.
  // This version uses a default time increment of 0.05.

  simulate_rprocess(t+0.0, t+1.0, prev, next, theta, .05, rng);
}



// // // rmeasure // // //
template<>
void stf<int, int>::rmeasure(int t, const std::vector<int>& state, std::vector<int>& meas, const std::vector<double>& theta, std::mt19937& rng) {
  const int ncity = state.size() / 5; // number of cities being analyzed
  const double rep = 1.0/(1.0+exp(-theta[7])); // reporting probability
  const double repOD = exp(theta[8]); // reporting overdispersion

  std::normal_distribution<double> normal(0.0, 1.0);

  for (int k = 0; k < ncity; k++) {
    meas[k] = int(state[5*k+3]*rep + normal(rng) * std::sqrt(state[5*k+3]*rep*(1-rep) + repOD*repOD*state[5*k+3]*state[5*k+3]*rep*rep + 1) + 0.5);
    if (meas[k] < 0) meas[k] = 0;
  }
}



// // // dmeasure // // //
template<>
double stf<int, int>::dmeasure(int t, const std::vector<int>& state, const std::vector<int>& meas, const std::vector<double>& theta, bool logarithm) {
  const int ncity = state.size() / 5; // number of cities being analyzed
  const double rep = 1.0/(1.0+exp(-theta[7])); // reporting probability
  const double repOD = exp(theta[8]); // reporting overdispersion

  double dm = 0.0; // log density of measurement

  for (int k = 0; k < ncity; k++) {
    dm += log(dnpmf(meas[k], state[5*k+3] * rep, std::sqrt(state[5*k+3]*rep*(1-rep) + repOD*repOD*state[5*k+3]*state[5*k+3]*rep*rep + 1)));
  }

  return (logarithm ? dm : exp(dm));
}



// // // intermediary processes for imfilter // // //
// // // imrprocess // // // 
template<>
void stf<int, int>::imrprocess(int t, int s, int im, const vector<int>& prev, vector<int>& next, const vector<double>& theta, double delta_t, std::mt19937& rng) {
  // the simulator of process distribution of x_{t + (s+1)/im} given x_{t + s/im}

  simulate_rprocess(t+(s+0.0)/im, t+(s+1.0)/im, prev, next, theta, delta_t, rng);
}


template<>
void stf<int, int>::imrprocess(int t, int s, int im, const vector<int>& prev, vector<int>& next, const vector<double>& theta, std::mt19937& rng) {
  // the simulator of process distribution of x_{t + (s+1)/im} given x_{t + s/im}
  // This version uses a default time increment of 0.05.

  simulate_rprocess(t+(s+0.0)/im, t+(s+1.0)/im, prev, next, theta, 0.05, rng);
}



// // // imdmeasure // // //
double marginal_dmeasure(int t, int k, const std::vector<int>& state, const std::vector<int>& meas, const std::vector<double>& theta) {
  // k: the space component for which the marginal likelihood of observation is evaluated
  const double rep = 1.0/(1.0+exp(-theta[7])); // reporting probability
  const double repOD = exp(theta[8]); // reporting overdispersion

  double dm; // log density of measurement

  dm = dnpmf(meas[k], state[5*k+3] * rep, std::sqrt(state[5*k+3]*rep*(1-rep) + repOD*repOD*state[5*k+3]*state[5*k+3]*rep*rep + 1));

  return dm;
}


template<>
vector<double> stf<int, int>::imdmeasure(int t, int s, vector<int> lookaheads, int im, const vector<int>& state, const vector<vector<int> >& meas, const vector<double>& theta, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+lookahead} given x_{t+s/im}
  // when lookahead = 1 and s = im, this function should be the same as the dmeasure function.
  // the length of lookaheads and the length of measurements should equal.

  int Jprime = 100; // number of particles used to estimate marginal likelihood of each component of y

  const int ncity = state.size() / 5; // number of cities being analyzed

  int rngseed = 834521; // use fixed seed to reduce variability in the ratio of two adjacent guide function evaluations
  std::mt19937 g_rng(rngseed); // rng used for guide function evaluation
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  vector<double> ll(lookaheads.size(), 0.0); // log likelihood estimates of observation components

  // the projections to (t+1) will be shared by all spatial coordinates. Store them to save computations.
  vector<vector<int> > g_xP_t_plus_1(Jprime, vector<int>(5*ncity, 0));
  for (int jp = 0; jp < Jprime; jp++)
    simulate_rprocess(t+(s+0.0)/im, t+lookaheads[0]+0.0, state, g_xP_t_plus_1[jp], theta, 0.25, g_rng);

  for (int k = 0; k < ncity; k++) {
    // set up particles
    vector<vector<int> > g_xF(Jprime, vector<int>(5*ncity, 0));
    vector<vector<int> > g_xP(Jprime, vector<int>(5*ncity, 0));

    for (size_t pos = 0; pos < lookaheads.size(); pos++) {
      vector<double> weights(Jprime);
      if (pos == 0) {
	g_xP = g_xP_t_plus_1;
      }
      for (int jp = 0; jp < Jprime; jp++) {
	// propagate particles
	if (pos > 0) {
	  simulate_rprocess(t+(s+0.0)/im, t+lookaheads[pos]+0.0, g_xF[jp], g_xP[jp], theta, 0.25, g_rng);
	}
	// evaluate particles
	weights[jp] = marginal_dmeasure(t+lookaheads[pos], k, g_xP[jp], meas[pos], theta);
      }
      double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
      ll[pos] += log(sum_weights/Jprime);
      if (pos < lookaheads.size() - 1) { // for the last lookahead, no need for resampling
	if (sum_weights > 0.0) { // resample
	  for (int jp = 0; jp < Jprime; jp++) {
	    weights[jp] /= sum_weights;
	  }
	  double p1 = uniform(g_rng) / Jprime;
	  double p2 = 0.0;
	  int j1 = -1;
	  for(int j=0; j<Jprime; j++) {
	    while(p2 < p1 && j1<Jprime-1)
	      p2 += weights[++j1];
	    g_xF[j] = g_xP[j1];
	    p1 += 1.0/Jprime;
	  }
	}
	else {
	  g_xF = g_xP;
	}
      }
    }
  }
  return ll;
}



template<>
double stf<int, int>::imdmeasure(int t, int s, int im, const vector<int>& state, const vector<int>& meas, const vector<double>& theta, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+1} given x_{t+s/im}
  // when s = im, this function should be the same as the dmeasure function.

  // set default values
  vector<int> lookaheads(1, 1);
  vector<vector<int> > observations(1, meas);

  vector<double> dm = imdmeasure(t, s, lookaheads, im, state, observations, theta, logarithm);

  return (logarithm ? dm[0] : exp(dm[0]));
}



template class stf<int, int>;
