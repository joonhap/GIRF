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

// this version (measlesmod_momentmatching.c) computes the guide function by adding the forecast variability to the variance of the measurement process.
// The forecast variability is estimated by simulating forward the rprocess multiple times and computing the sample statistics (e.g., inter-quartile range or variance).

bool annealed = true; // the forecast likelihoods are raised to fractional powers depending on the forecast time length

// use city-specific measles case reporting ratio estimated by dividing the total number of reported cases in year 1948--64 by the total number of births in year 1944-1960.
double repdata[] =   {0.52, 0.56, 0.69, 0.62, 0.64, 0.59, 0.58, 0.62, 0.58, 0.69, 0.69, 0.61, 0.64, 0.59, 0.61, 0.68, 0.61, 0.51, 0.51, 0.58, 0.61, 0.59, 0.68, 0.59, 0.70, 0.64, 0.57, 0.56, 0.50, 0.44, 0.67, 0.59, 0.77, 0.61, 0.51, 0.54, 0.68, 0.55, 0.69, 0.55};


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
  // the initial value parameter vector defines the proportion of each compartment at t=0
  // as a precaution, check if the length of the first dimension of init_tempate and thetaSwarm match.
  if (init_template.size() != thetaSwarm.size()) {
    cout << "Error in initializing particle states: the length of init_template and that of thetaSwarm do not match." << endl;
    exit(1);
  }
  initstates = init_template;

  int R = initstates.size();
  int J = initstates[0].size();
  int ncity = initstates[0][0].size()/5;
  for (int i=0; i<R; i++)
    for (int j=0; j<J; j++) 
      for (int k=0; k<ncity; k++) {
	initstates[i][j][5*k+0] = initstates[i][j][5*k+4] / (1.0+exp(-thetaSwarm[i][j][pos_ivp[3*k+0]]));
	initstates[i][j][5*k+1] = initstates[i][j][5*k+4] / (1.0+exp(-thetaSwarm[i][j][pos_ivp[3*k+1]]));
	initstates[i][j][5*k+2] = initstates[i][j][5*k+4] / (1.0+exp(-thetaSwarm[i][j][pos_ivp[3*k+2]]));
      }
}


// // // deterministic skeleton // // //
void simulate_deterministic_process(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t) {
  // process simulator for deterministic skeleton
  // prev, next: state vector. Length 5*k. Consists of k copies of (S,E,I,weekly_I_to_R_transition,P). P (population) is (piecewise) constant.
  // weekly_I_to_R_transition is the weekly aggregate transition from I to R.
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

  const double tol=1e-5;
  const double dpu=365.25/26.0; // days per time unit interval

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end
    }

    // initialize biweekly total I to R transitions (if an integer is in the interval [time, time+t_inc), set the biweekly sum to zero.)
    if (int(time-tol) < int(time+t_inc-tol)) {
      for (int k=0; k < ncity; k++)
	pv[5*k+3] = 0; 
    }

    int year = int(time * dpu / 365.25) + startyear; 
    int day = int(time * dpu + tol - (year-startyear)*365.25);

    // change population size at the start of year
    if (year < int((time+t_inc) * dpu / 365.25) + startyear) {
      // if population data is not found for the new year, exit the program
      if ((int((time+t_inc) * dpu / 365.25) + startyear < pop_data_first_year) || (int((time+t_inc) * dpu / 365.25) + startyear > pop_data_last_year)) {
	cout << "Population data not found for the new year (Population data needed for year " << (int((time+t_inc) * dpu / 365.25) + startyear) << ".) Exit the program." << endl;
	exit(1);
      }
      for (int k=0; k<ncity; k++)
	pv[5*k+4] = pop_data[int((time+t_inc) * dpu / 365.25) + startyear-pop_data_first_year][cityindex[k]];
    }
    
    // recruit incremental susceptibles
    double day_double = time*dpu - int(time*dpu/365.25) * 365.25;
    for (int k=0; k<ncity; k++) 
      pv[5*k] += (int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * (day_double+t_inc*dpu) / 365.25) - int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * day_double / 365.25));

    // recruit susceptibles on school start day
    if((day < schoolstart) & (schoolstart <= int((time+t_inc)*dpu + tol - (year-startyear)*365.25))) { // make sure not to recruit susceptibles multiples time a day
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
      
      double gamma_SE = t_inc;
      double gamma_SD = t_inc;
      double gamma_EI = t_inc;
      double gamma_ED = t_inc;
      double gamma_IR = t_inc;
      double gamma_ID = t_inc;
      
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
  // process simulator for rprocess
  // prev, next: state vector. Length 5*k. Consists of k copies of (S,E,I,weekly_I_to_R_transition,P). P (population) is (piecewise) constant.
  // weekly_I_to_R_transition is the weekly aggregate transition from I to R.
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

  const double tol=1e-5;
  const double dpu=365.25/26.0; // days in a biweek

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end.
    }

    // initialize weekly total I to R transitions (if an integer is in the interval [time, time+t_inc), set the weekly sum to zero.)
    if (int(time-tol) < int(time+t_inc-tol)) {
      for (int k=0; k < ncity; k++)
	pv[5*k+3] = 0; 
    }

    int year = int(time * dpu / 365.25) + startyear; 
    int day = int(time * dpu + tol - (year-startyear)*365.25);

    // change population size at the start of year
    if (year < int((time+t_inc) * dpu / 365.25) + startyear) {
      // if population data is not found for the new year, exit the program
      if ((int((time+t_inc) * dpu / 365.25) + startyear < pop_data_first_year) || (int((time+t_inc) * dpu / 365.25) + startyear > pop_data_last_year)) {
	cout << "Population data not found for the new year (Population data needed for year " << (int((time+t_inc) * dpu / 365.25) + startyear) << ".) Exit the program." << endl;
	exit(1);
      }
      for (int k=0; k<ncity; k++)
	pv[5*k+4] = pop_data[int((time+t_inc) * dpu / 365.25) + startyear-pop_data_first_year][cityindex[k]]; 
    }

    // recruit incremental susceptibles
    double day_double = time*dpu - int(time*dpu/365.25) * 365.25;

    for (int k=0; k<ncity; k++) {
      pv[5*k] += (int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * (day_double+t_inc*dpu) / 365.25) - int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * day_double / 365.25));
    }

    // recruit susceptibles on school start day
    if((day < schoolstart) & (schoolstart <= int((time+t_inc)*dpu + tol - (year-startyear)*365.25))) { // make sure not to recruit susceptibles multiples time a day
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
      
      std::gamma_distribution<double> gamma(t_inc/sigma2, sigma2);
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
  vector<double> rep_all; // reporting probability for all cities in the data file
  rep_all.assign(repdata, repdata+NCITY);
  vector<double> rep(ncity); // reporting probabilities for the cities being analyzed
  for (int k = 0; k < ncity; k++)
    rep[k] = rep_all[cityindex[k]];
  const double repOD = exp(theta[8]); // reporting overdispersion

  std::normal_distribution<double> normal(0.0, 1.0);

  for (int k = 0; k < ncity; k++) {
    meas[k] = int(state[5*k+3]*rep[cityindex[k]] + normal(rng) * std::sqrt(state[5*k+3]*rep[cityindex[k]]*(1-rep[cityindex[k]]) + repOD*repOD*state[5*k+3]*state[5*k+3]*rep[cityindex[k]]*rep[cityindex[k]] + 1) + 0.5);
    if (meas[k] < 0) meas[k] = 0;
  }
}



// // // dmeasure // // //
template<>
double stf<int, int>::dmeasure(int t, const std::vector<int>& state, const std::vector<int>& meas, const std::vector<double>& theta, bool logarithm) {
  const int ncity = state.size() / 5; // number of cities being analyzed
  vector<double> rep_all; // reporting probability for all cities in the data file
  rep_all.assign(repdata, repdata+NCITY);
  vector<double> rep(ncity); // reporting probabilities for the cities being analyzed
  for (int k = 0; k < ncity; k++)
    rep[k] = rep_all[cityindex[k]];
  const double repOD = exp(theta[8]); // reporting overdispersion

  double dm = 0.0; // log density of measurement

  for (int k = 0; k < ncity; k++) {
    dm += log(dnpmf(meas[k], state[5*k+3] * rep[k], std::sqrt(state[5*k+3]*rep[k]*(1-rep[k]) + repOD*repOD*state[5*k+3]*state[5*k+3]*rep[k]*rep[k] + 1)));
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
template<>
vector<double> stf<int, int>::imdmeasure(int t, int s, vector<int> lookaheads, int im, const vector<int>& state, const vector<vector<int> >& meas, const vector<double>& theta, vector<vector<double> >& forecast_var, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+lookahead} given x_{t+s/im}
  // note that the input s will take values between 1 and im (does not take value zero). Check stf_parallel.tpp
  // we use projection of state to t+lookahead to compute pseudo-predictive likelihood
  // when lookahead = 1 and s = im, this function should be the same as the dmeasure function.
  // the length of lookaheads and the lengthe measurements should equal.

  const int ncity = state.size()/5;
  vector<double> rep_all; // reporting probability for all cities in the data file
  rep_all.assign(repdata, repdata+NCITY);
  vector<double> rep(ncity); // reporting probabilities for the cities being analyzed
  for (int k = 0; k < ncity; k++)
    rep[k] = rep_all[cityindex[k]];
  const double repOD = exp(theta[8]); // reporting overdispersion

  // measurement density (log)
  vector<double> dm(lookaheads.size(), 0.0);
  vector<double> exp_dm(lookaheads.size(), 0.0);

  int lasize = lookaheads.size(); // the total number of lookaheads

  bool estimate_forecast_var = (s==1); // whether we update the estimated forecast variability at the current time point. Should always be re-estimated at s=1.
  int s_last_fve = 1; // last intermediate step at which the forecast variability was estimated. Should correspond with how estimate_forecast_var was set.
  if (estimate_forecast_var) {
    int nforecast = 40; // number of random forecasts to make
    vector<vector<vector<int> > > forecasts(nforecast, vector<vector<int> >(lasize)); // array structure: [forecast_id][lookahead_entry_id][space_dim]
    std::mt19937 rng(4255427);
    double delta_t = 0.25;
    for (int fc = 0; fc < nforecast; fc++) {
      simulate_rprocess(t+(s+0.0)/im, double(t+lookaheads[0]), state, forecasts[fc][0], theta, delta_t, rng);
      for (int pos = 1; pos < lasize; pos++) {
	simulate_rprocess(double(t+lookaheads[pos-1]), double(t+lookaheads[pos]), forecasts[fc][pos-1], forecasts[fc][pos], theta, delta_t, rng);
      }
    }

    for (int pos = 0; pos < lasize; pos++) {
      forecast_var[pos].resize(ncity);
      for (int k = 0; k < ncity; k++) {
	vector<double> repCase_forecasts(nforecast); // forecasts of reported biweekly measles cases for this time point and space dimension
	for (int fc = 0; fc < nforecast; fc++) {
	  repCase_forecasts[fc] = forecasts[fc][pos][5*k+3] * rep[k];
	}
	quicksort(repCase_forecasts, 0, nforecast-1);
	double iqr = repCase_forecasts[int(nforecast*.75)] - repCase_forecasts[int(nforecast*.25)]; // inter-quartile range
	double var_to_iqr2 = .55 * (1+2/sqrt(nforecast)); // distribution-specific ratio between the variance and the squared iqr (for scale families) (~0.55 for normal). Take a conservative number that takes into account the variability in the estimates of forecast variability itself. This translates into inflating the number by a factor (1+2*sqrt(nforecast)).
	forecast_var[pos][k] = iqr*iqr*var_to_iqr2;
      }
    }
  }

  vector<double> forc_base(state.size()); // initial state for making deterministic forcast
  for (size_t k = 0; k < state.size(); k++) {
    forc_base[k] = double(state[k]);
  }
  vector<double> det_forc; // deterministic forecast

  for (int pos = 0; pos < lasize; pos++) {
    double delta_t = 0.25;
    simulate_deterministic_process((pos == 0 ? t+(s+0.0)/im : double(t+lookaheads[pos-1])), double(t+lookaheads[pos]), forc_base, det_forc, theta, delta_t);
    forc_base = det_forc;
    double scale_fv = (lookaheads[pos] - (s+0.0)/im) / std::max((lookaheads[pos] - (s_last_fve+0.0)/im), 1e-6); // scaling factor for forecast variance
    for (int k = 0; k < ncity; k++) {
      double meas_variance = det_forc[5*k+3]*rep[k]*(1-rep[k]) + repOD*repOD*det_forc[5*k+3]*det_forc[5*k+3]*rep[k]*rep[k] + 1; // measurement variance
      dm[pos] += log(dnpmf(meas[pos][k], det_forc[5*k+3] * rep[k], std::sqrt(meas_variance + scale_fv * forecast_var[pos][k])));
    }
    exp_dm[pos] = exp(dm[pos]);
  }
  return (logarithm ? dm : exp_dm);
}




template<>
double stf<int, int>::imdmeasure(int t, int s, int im, const vector<int>& state, const vector<int>& meas, const vector<double>& theta, vector<double>& forecast_var_1step, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+1} given x_{t+s/im} (ONLY to the immediate next time point)
  // when s = im, this function should be the same as the dmeasure function.

  // set default values
  vector<int> lookaheads(1, 1);
  vector<vector<int> > observations(1, meas);
  vector<vector<double> > forecast_var(1, forecast_var_1step);

  vector<double> dm = imdmeasure(t, s, lookaheads, im, state, observations, theta, forecast_var, true);

  return (logarithm ? dm[0] : exp(dm[0]));
}



template class stf<int, int>;
