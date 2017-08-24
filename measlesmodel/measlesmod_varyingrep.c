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
  

  /*
  int ncity = initstates[0][0].size()/5;
  int startyear = int(thetaSwarm[0][0][16]+.5); // start year (the year at which t=0)
  int entryage = int(thetaSwarm[0][0][15]+.5); // school entry age (given as double in theta)
  vector<double> birthrate(ncity, 0.0); // average birth rate over the time period (with a lag of entryage)
  vector<double> avgPop(ncity, 0.0); // average population over the time period
  for (int k=0; k<ncity; k++) {
    for (int y=0; y<ntimes/52; y++) {
      birthrate[k] += birth_data[startyear-entryage-birth_data_first_year+y][cityindex[k]];
      avgPop[k] += pop_data[startyear-pop_data_first_year+y][cityindex[k]]; 
    }
    birthrate[k] /= ((ntimes/52)*365.25/7.0);
    avgPop[k] /= (ntimes/52);
  }

  vector<double> rep(20); // reporting probability
  double repdata[] = {.75, .60, .13, .46, .63, .78, .31, .65, .70, .60, .60, .58, .61, .63, .67, .65, .55, .49, .56, .49};
  rep.assign(repdata, repdata+20);

  double hol_prop = .261; // the proportion of school holidays in a year
  for (size_t i=0; i<initstates.size(); i++)
    for (size_t j=0; j<initstates[0].size(); j++) {
      double R0 = exp(thetaSwarm[i][j][0]); // basic reproduction number
      double amp = 1.0/(1.0+exp(-thetaSwarm[i][j][1]))*.68; // seasonal amplitude
      double gen = exp(thetaSwarm[i][j][4]); // generation cycle length (1/nu_IR+1/nu_EI, log-transformed)
      double infec_period_prop = 1.0/(1.0+exp(-thetaSwarm[i][j][5])); // infectious period proportion
      double nu_IR = 1.0/(gen*infec_period_prop); // 1/infectious period (log-transformed)
      double nu_EI = 1.0/(gen*(1.0-infec_period_prop)); // 1/latent period
      double beta_s = (1.0+2*(1.0-school_prop)*amp)*R0*nu_IR; // school time transmission coefficient
      double hol_coeff = (1-2*school_prop*amp) / (1+2*(1-school_prop)*amp); // holiday effect on transmission coefficient (<1)
      double alpha = thetaSwarm[i][j][2]; // mixing exponent
      double mu = thetaSwarm[i][j][3]; // mortality rate
      double beta_m = beta_s * ( (1.0-school_prop)*hol_coeff + school_prop ); // annual average beta
      for (size_t k=0; k<initstates[0][0].size()/5; k++) {
	// initstates[i][j][5*k] = int(initstates[i][j][5*k+4] * 1.0/(1.0+exp(-thetaSwarm[i][j][pos_ivp[0]]))); // old version
	double stationary_I = birthrate[k]*nu_EI/(nu_EI+mu)/(nu_IR+mu) - avgPop[k]*mu/beta_m; // the size of I in stationary distribution for the computation of stationary size of S (under the assumption that alpha = 1.0. The alpha value has little effect on the stationary size of I)
	initstates[i][j][5*k] = (nu_EI+mu)*(nu_IR+mu)/nu_EI*std::pow(stationary_I, 1.0-alpha)*std::pow(avgPop[k], alpha)/beta_m; // set the starting value of S at the stationary value given the current set of parameters.
	vector<double> a12(2); // the (1,2) entry of rate matrix A which satisfies d/dt(E I)^t = A * (E I)^t (for the first and the second weeks)
	a12[0] = initstates[i][j][5*k] * beta_s * hol_eff(0, hol_coeff) / avgPop[k];
	a12[1] = initstates[i][j][5*k] * beta_s * hol_eff(7, hol_coeff) / avgPop[k];
	vector<double> disc(2); // discriminant of the characteristic polynomial
	disc[0] = sqrt((nu_EI - nu_IR)*(nu_EI - nu_IR) + 4 * nu_EI * a12[0]); // for A(0)
	disc[1] = sqrt((nu_EI - nu_IR)*(nu_EI - nu_IR) + 4 * nu_EI * a12[1]); // for A(1)
	vector<vector<double> > lm(2, vector<double>(2)); // eigenvalues of matrix A
	lm[0][0] = .5 * (-nu_EI - nu_IR - 2*mu + disc[0]); // for A(0)
	lm[0][1] = .5 * (-nu_EI - nu_IR - 2*mu - disc[0]);
	lm[1][0] = .5 * (-nu_EI - nu_IR - 2*mu + disc[1]); // for A(1)
	lm[1][1] = .5 * (-nu_EI - nu_IR - 2*mu - disc[1]);

	vector<vector<vector<double> > > C(2, vector<vector<double> >(2, vector<double>(2, 0.0))); 
	// the matrix C[0] is such that (N_IR(t=1)-N_IR(0)) = nu_IR * (0 1) * C[0] * (E(t=0) I(0))^t and
	// C[1] such that (N_IR(2) - N_IR(1)) = nu_IR * (0 1) * C[1] * (E(0) I(0))^t
	// C[0] = (INT_0^1 exp(A(0) s) ds), C[1] = (INT_0^1 exp(A(1) s) ds) * exp(A(0))
	// eigendecomposition A = B L B^{-1}
	vector<vector<vector<double> > > B(2, vector<vector<double> >(2, vector<double>(2)));
	vector<vector<vector<double> > > Bi(2, vector<vector<double> >(2, vector<double>(2))); // B-inverse
	B[0][0][0] = -a12[0];
	B[0][0][1] = -a12[0];
	B[0][1][0] = -(nu_EI+mu) - lm[0][0];
	B[0][1][1] = -(nu_EI+mu) - lm[0][1];
	B[1][0][0] = -a12[1];
	B[1][0][1] = -a12[1];
	B[1][1][0] = -(nu_EI+mu) - lm[1][0];
	B[1][1][1] = -(nu_EI+mu) - lm[1][1];
	vector<double> detB(2); // determinant of B
	detB[0] = B[0][0][0]*B[0][1][1] - B[0][0][1]*B[0][1][0];
	detB[1] = B[1][0][0]*B[1][1][1] - B[1][0][1]*B[1][1][0];
	Bi[0][0][0] = B[0][1][1]/detB[0];
	Bi[0][0][1] = -B[0][0][1]/detB[0];
	Bi[0][1][0] = -B[0][1][0]/detB[0];
	Bi[0][1][1] = B[0][0][0]/detB[0];
	Bi[1][0][0] = B[1][1][1]/detB[1];
	Bi[1][0][1] = -B[1][0][1]/detB[1];
	Bi[1][1][0] = -B[1][1][0]/detB[1];
	Bi[1][1][1] = B[1][0][0]/detB[1];
	// compute C matrices
	vector<double> d1(2); // eigenvalues of C(t=0)
	d1[0] = (exp(lm[0][0])-1)/lm[0][0];
	d1[1] = (exp(lm[0][1])-1)/lm[0][1];
	C[0][0][0] = B[0][0][0]*d1[0]*Bi[0][0][0] + B[0][0][1]*d1[1]*Bi[0][1][0];
	C[0][0][1] = B[0][0][0]*d1[0]*Bi[0][0][1] + B[0][0][1]*d1[1]*Bi[0][1][1];
	C[0][1][0] = B[0][1][0]*d1[0]*Bi[0][0][0] + B[0][1][1]*d1[1]*Bi[0][1][0];
	C[0][1][1] = B[0][1][0]*d1[0]*Bi[0][0][1] + B[0][1][1]*d1[1]*Bi[0][1][1];
	vector<double> d2(2); // eigenvalues of exp(A(0))
	vector<vector<double> > expA0(2, vector<double>(2)); // exp(A(0))	
	d2[0] = exp(lm[0][0]);
	d2[1] = exp(lm[0][1]);
	expA0[0][0] = B[0][0][0]*d2[0]*Bi[0][0][0] + B[0][0][1]*d2[1]*Bi[0][1][0];
	expA0[0][1] = B[0][0][0]*d2[0]*Bi[0][0][1] + B[0][0][1]*d2[1]*Bi[0][1][1];
	expA0[1][0] = B[0][1][0]*d2[0]*Bi[0][0][0] + B[0][1][1]*d2[1]*Bi[0][1][0];
	expA0[1][1] = B[0][1][0]*d2[0]*Bi[0][0][1] + B[0][1][1]*d2[1]*Bi[0][1][1];
	vector<double> d3(2); // eigenvalues of INT_0^1 exp(A(1) s) ds
	vector<vector<double> > intExpA1(2, vector<double>(2)); // INT_0^1 exp(A(1) s) ds
	d3[0] = (exp(lm[1][0])-1)/lm[1][0];
	d3[1] = (exp(lm[1][1])-1)/lm[1][1];
	intExpA1[0][0] = B[1][0][0]*d3[0]*Bi[1][0][0] + B[1][0][1]*d3[1]*Bi[1][1][0];
	intExpA1[0][1] = B[1][0][0]*d3[0]*Bi[1][0][1] + B[1][0][1]*d3[1]*Bi[1][1][1];
	intExpA1[1][0] = B[1][1][0]*d3[0]*Bi[1][0][0] + B[1][1][1]*d3[1]*Bi[1][1][0];
	intExpA1[1][1] = B[1][1][0]*d3[0]*Bi[1][0][1] + B[1][1][1]*d3[1]*Bi[1][1][1];
	// C(1) = (INT_0^1 exp(A(1) s) ds) * exp(A(0))
	C[1][0][0] = intExpA1[0][0]*expA0[0][0] + intExpA1[0][1]*expA0[1][0];
	C[1][0][1] = intExpA1[0][0]*expA0[0][1] + intExpA1[0][1]*expA0[1][1];
	C[1][1][0] = intExpA1[1][0]*expA0[0][0] + intExpA1[1][1]*expA0[1][0];
	C[1][1][1] = intExpA1[1][0]*expA0[0][1] + intExpA1[1][1]*expA0[1][1];

	// solve equations  nu_IR * (0 1) * C(0) * (E(0) I(0))^t = N_IR(1) - N_IR(0) and
	// nu_IR * (0 1) * C(1) * (E(0) I(0))^t = N_IR(2) - N_IR(1)
	// call F = nu_IR * (C(0)_21, C(0)_22 // C(1)_21, C(1)_22) and then F * (E(0) I(0))^t = (N_IR(1) - N_IR(0), N_IR(2) - N_IR(1))^t
	vector<vector<double> > F(2, vector<double>(2));
	F[0][0] = nu_IR * C[0][1][0];
	F[0][1] = nu_IR * C[0][1][1];
	F[1][0] = nu_IR * C[1][1][0];
	F[1][1] = nu_IR * C[1][1][1];
	double detF = F[0][0]*F[1][1] - F[0][1]*F[1][0];
	initstates[i][j][5*k+1] = (F[1][1] * obs[0][k] - F[0][1] * obs[1][k]) / rep[cityindex[k]] / detF;
	initstates[i][j][5*k+2] = (-F[1][0] * obs[0][k] + F[0][0] * obs[1][k]) / rep[cityindex[k]] / detF;
	initstates[i][j][5*k+1] *= (initstates[i][j][5*k+1] < 0 ? 0 : 1);
	initstates[i][j][5*k+2] *= (initstates[i][j][5*k+2] < 0 ? 0 : 1);
      }
    }
  */
}


// // // deterministic skeleton // // //
void simulate_deterministic_process(double t_start, double t_end, const vector<double>& prev, vector<double>& next, const vector<double>& theta, double delta_t) {
  // process simulator for either deterministic skeleton or rprocess (as specified by the parameter 'deterministic')
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

  const double tol=1e-6;

  for(double time = t_start; time < t_end-tol; time += delta_t) {
    double t_inc = delta_t; // time increment for this step
    if (time > t_end-delta_t+tol) {
      t_inc = t_end - time; // if the last time increment is less than delta_t, make it such that the final time equals t_end
    }

    // initialize weekly total I to R transitions (if an integer is in the interval [time, time+t_inc), set the weekly sum to zero.)
    if (int(time-tol) < int(time+t_inc-tol)) {
      for (int k=0; k < ncity; k++)
	pv[5*k+3] = 0; 
    }

    int year = int(time * 7.0 / 365.25) + startyear; 
    int day = int(time * 7.0 + tol - (year-startyear)*365.25);

    // change population size at the start of year
    if (year < int((time+t_inc) * 7.0 / 365.25) + startyear) {
      // if population data is not found for the new year, exit the program
      if ((int((time+t_inc) * 7.0 / 365.25) + startyear < pop_data_first_year) || (int((time+t_inc) * 7.0 / 365.25) + startyear > pop_data_last_year)) {
	cout << "Population data not found for the new year (Population data needed for year " << (int((time+t_inc) * 7.0 / 365.25) + startyear) << ".) Exit the program." << endl;
	exit(1);
      }
      for (int k=0; k<ncity; k++)
	pv[5*k+4] = pop_data[int((time+t_inc) * 7.0 / 365.25) + startyear-pop_data_first_year][cityindex[k]];
    }
    
    // recruit incremental susceptibles
    double day_double = time*7.0 - int(time*7.0/365.25) * 365.25;
    for (int k=0; k<ncity; k++) 
      pv[5*k] += (int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * (day_double+t_inc*7.0) / 365.25) - int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * day_double / 365.25));

    // recruit susceptibles on school start day
    if((day < schoolstart) & (schoolstart <= int((time+t_inc)*7.0 + tol - (year-startyear)*365.25))) { // make sure not to recruit susceptibles multiples time a day
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
  // process simulator for either deterministic skeleton or rprocess (as specified by the parameter 'deterministic')
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

  const double tol=1e-6;

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

    int year = int(time * 7.0 / 365.25) + startyear; 
    int day = int(time * 7.0 + tol - (year-startyear)*365.25);

    // change population size at the start of year
    if (year < int((time+t_inc) * 7.0 / 365.25) + startyear) {
      // if population data is not found for the new year, exit the program
      if ((int((time+t_inc) * 7.0 / 365.25) + startyear < pop_data_first_year) || (int((time+t_inc) * 7.0 / 365.25) + startyear > pop_data_last_year)) {
	cout << "Population data not found for the new year (Population data needed for year " << (int((time+t_inc) * 7.0 / 365.25) + startyear) << ".) Exit the program." << endl;
	exit(1);
      }
      for (int k=0; k<ncity; k++)
	pv[5*k+4] = pop_data[int((time+t_inc) * 7.0 / 365.25) + startyear-pop_data_first_year][cityindex[k]]; 
    }

    // recruit incremental susceptibles
    double day_double = time*7.0 - int(time*7.0/365.25) * 365.25;
    for (int k=0; k<ncity; k++) {
      pv[5*k] += (int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * (day_double+t_inc*7.0) / 365.25) - int((1.0-cohort_entry) * birth_data[year-entryage-birth_data_first_year][cityindex[k]] * day_double / 365.25));
    }

    // recruit susceptibles on school start day
    if((day < schoolstart) & (schoolstart <= int((time+t_inc)*7.0 + tol - (year-startyear)*365.25))) { // make sure not to recruit susceptibles multiples time a day
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
  vector<double> rep; // reporting probability
  double repdata[] = {.75, .60, .13, .46, .63, .78, .31, .65, .70, .60, .60, .58, .61, .63, .67, .65, .55, .49, .56, .49};
  rep.assign(repdata, repdata+20);
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
  vector<double> rep; // reporting probability
  double repdata[] = {.75, .60, .13, .46, .63, .78, .31, .65, .70, .60, .60, .58, .61, .63, .67, .65, .55, .49, .56, .49};
  rep.assign(repdata, repdata+20);
  const double repOD = exp(theta[8]); // reporting overdispersion

  double dm = 0.0; // log density of measurement

  for (int k = 0; k < ncity; k++) {
    dm += log(dnpmf(meas[k], state[5*k+3] * rep[cityindex[k]], std::sqrt(state[5*k+3]*rep[cityindex[k]]*(1-rep[cityindex[k]]) + repOD*repOD*state[5*k+3]*state[5*k+3]*rep[cityindex[k]]*rep[cityindex[k]] + 1)));
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
vector<double> stf<int, int>::imdmeasure(int t, int s, vector<int> lookaheads, int im, const vector<int>& state, const vector<vector<int> >& meas, const vector<double>& theta, bool logarithm) {
  // the density of auxiliary measurement distribution of y_{t+lookahead} given x_{t+s/im}
  // we use projection of state to t+1 to compute pseudo-predictive likelihood
  // when lookahead = 1 and s = im, this function should be the same as the dmeasure function.
  // the length of lookaheads and the lengthe measurements should equal.

  const int ncity = state.size() / 5; // number of cities being analyzed 
  const double gen = exp(theta[4]); // generation cycle (1/nu_EI+1/nu_IR, log-transformed)
  const double infec_period_prop = 1.0/(1.0+exp(-theta[5])); // infectious period proportion
  const double nu_IR = 1.0/(gen*infec_period_prop); // 1/infectious period (log-transformed)
  const double nu_EI = 1.0/(gen*(1.0-infec_period_prop)); // 1/latent period
  const double sigma2 = exp(theta[6]); // variability constant (for gamma process, sigma^2. see Box 1 in He et al 2009) (log-transformed)
  vector<double> rep; // reporting probability
  double repdata[] = {.75, .60, .13, .46, .63, .78, .31, .65, .70, .60, .60, .58, .61, .63, .67, .65, .55, .49, .56, .49};
  rep.assign(repdata, repdata+20);
  const double repOD = exp(theta[8]); // reporting overdispersion

  const double R0 = exp(theta[0]); // basic reproduction number
  const double amp = 1.0/(1.0+exp(-theta[1]))*.68; // seasonal amplitude
  const double beta_s = (1.0+2*(1.0-school_prop)*amp)*R0*nu_IR; // school time transmission coefficient
  const double hol_coeff = (1-2*school_prop*amp) / (1+2*(1-school_prop)*amp); // holiday effect on transmission coefficient (<1)
  const double alpha = theta[2]; // mixing exponent
  const double mu = theta[3]; // mortality rate
  const int startyear = int(theta[16]+.5); // start year (the year at which t=0)
  const double tol = 1e-6; // tolerance for machine precision

  // measurement density (log)
  vector<double> dm(lookaheads.size(), 0.0);
  vector<double> exp_dm(lookaheads.size(), 0.0);

  // make projections to quarter points (numbers of the form n/4)
  int lasize = lookaheads.size(); // the total number of lookaheads
  int nproj_max = lookaheads[lasize-1]; // the maximum lookahead
  int q_first = int(s/double(im)*4-tol)+1; // the first quarter point to make projection to
  int nproj_q = nproj_max*4-q_first+1; // number of projections to make
  vector<vector<double> > proj_q(nproj_q+1, vector<double>(5*ncity)); // state simulated to quarter time points (the zeroth entry is the current state)
  for (int k=0; k<5*ncity; k++)
    proj_q[0][k] = double(state[k]);
  simulate_deterministic_process(t+(s+0.0)/im, t+q_first/4.0, proj_q[0], proj_q[1], theta, 0.05);
  for (int qpoint = 1; qpoint < nproj_q; qpoint++) 
    simulate_deterministic_process(t+(q_first+qpoint-1)/4.0, t+(q_first+qpoint)/4.0, proj_q[qpoint], proj_q[qpoint+1], theta, 0.05);

  // compute the variance of N_IR(t_{n+b})-N_IR(t_{n+b-1}) conditioned on the present (t_{n,s})
  // The matrix A(t) in the system of equations d/dt(E(t) I(t))^t = A(t) (E(t) I(t))^t. see supplementary material of the paper.
  // A(t) = (-(nu_EI+m), beta*S*(I/P)^alpha/I ; nu_EI, -(nu_IR+m))
  // compute 1. the entry a_12 in the rate of change matrix A(t), 2. the eigenvalues of A(t), lm.
  vector<vector<double> > a12(nproj_max, vector<double>(ncity)); // the (1,2)-th entry of A matrix
  vector<vector<vector<double> > > lm(nproj_max, vector<vector<double> >(ncity, vector<double>(2)));
  for (int a = 0; a < nproj_max; a++) { // compute A matrix for the interval [t_{n+a}, t_{n+a+1}]
    int year = int((t+a+1) * 7.0 / 365.25) + startyear; 
    int day = int((t+a+1) * 7.0 + tol - (year-startyear)*365.25);
    double beta = beta_s * hol_eff(day, hol_coeff);
    for (int k=0; k<ncity; k++) {
      int start_q = (4*a-q_first+1 < 0 ? 0 : 4*a-q_first+1); // starting position for this week in proj_q
      int end_q = 4*(a+1)-q_first+1; // end position for this week in proj_q
      double meanS = .5 * (proj_q[start_q][5*k] + proj_q[end_q][5*k]);
      double meanOoP = .5 * ((proj_q[start_q][5*k+2]<tol ? 1.0/proj_q[start_q][5*k+4] : pow(proj_q[start_q][5*k+2]/proj_q[start_q][5*k+4], alpha) / proj_q[start_q][5*k+2]) + (proj_q[end_q][5*k+2]<tol ? 1.0/proj_q[end_q][5*k+4] : pow(proj_q[end_q][5*k+2]/proj_q[end_q][5*k+4], alpha) / proj_q[end_q][5*k+2]));
      a12[a][k] = meanS * beta * meanOoP;
      double disc = sqrt((nu_EI - nu_IR)*(nu_EI - nu_IR) + 4 * nu_EI * a12[a][k]); // discriminant of the characteristic polynomial
      lm[a][k][0] = .5 * (-nu_EI - nu_IR - 2*mu + disc);
      lm[a][k][1] = .5 * (-nu_EI - nu_IR - 2*mu - disc);
    }
  }
  
  vector<vector<double> > proj(lasize, vector<double>(5*ncity)); // state simulated to integer lookahead points
  for (int pos=0; pos<lasize; pos++) 
    proj[pos] = proj_q[lookaheads[pos]*4-q_first+1];

  // compute approximate predictive likelihood
  for (int pos=0; pos<lasize; pos++) {
    int b = lookaheads[pos];
    vector<double> var(ncity, 0.0); // variance of projection
    vector<vector<vector<double> > > C(ncity, vector<vector<double> >(2, vector<double>(2, 0.0))); 

    for (int k=0; k<ncity; k++) {
      // compute variance of projection
      for (int qpoint = lookaheads[pos]*4-q_first+1; qpoint >= 1; qpoint--) { // start adding conditional variance from the end
	// compute conditional variances : E( Var( E(N_IR(t_{n+b})-N_IR(t_{n+b-1}) | t_{n+(qpoint+q_first-1)/4} ) | t_{n+(qpoint+q_first-2)/4} ) )
	// t_{n+(qpoint+q_first-1)/4} = t_{n+a+a'/4}, a, a' integer, a'<= 4
	int a = (q_first+qpoint-1)/4; // the interval [t_{n+(q_first+qpoint-1)/4.0}, t_{n+(q_first+qpoint)/4.0}] is in the interval [t_{n+a}, t_{n+a+1}].
	a = (a>=b ? b-1 : a); // at the last step (a=b), the C matrix is zero, and no need to compute it. 
	// the matrix C is such that (N_IR(t_{n+b})-N_IR(max(t_{n+b-1},t_{n+a+a'/4})) = (0 1)*C*(E(t_{n+a+a'/4}) I(t_{n+a+a'/4}))^t
	// C = (INT_0^u exp(A_{t_{n+b-1}} s) ds)*exp(A_{t_{n+b-2}})*...*exp(A_{t_{n+a}} (1-a'/4)) where u = t_{n+b}-max(t_{n+a+a'/4},t_{n+b-1})
	// A[a] = B L B^{-1} : diagonalization
	vector<vector<double> > B(2, vector<double>(2));
	vector<vector<double> > Bi(2, vector<double>(2)); // B inverse
	B[0][0] = -nu_IR-mu-lm[a][k][0];
	B[0][1] = -nu_IR-mu-lm[a][k][1];
	B[1][0] = -nu_EI;
	B[1][1] = -nu_EI;
	double detB = B[0][0]*B[1][1] - B[0][1]*B[1][0];
	detB = (detB<tol && detB>-tol ? tol*(detB>=0 ? 1 : -1) : detB); // regularize detB to be at least tol in absolute value
	Bi[0][0] = B[1][1]/detB;
	Bi[0][1] = -B[0][1]/detB;
	Bi[1][0] = -B[1][0]/detB;
	Bi[1][1] = B[0][0]/detB;
	// compute C
	if ((q_first+qpoint-1 >= 4*(b-1)) & (q_first+qpoint-1 < 4*b)) { 
	  // C[k] = B*(INT_0^u exp(Lm s) ds)*Bi where u = t+b-(t+a+a'/4) and Lm = diag(lm1 lm2)
	  double u = (4*b-q_first-qpoint+1)/4.0;
	  // (INT_0^u exp(Lm s) ds) = diag(d1 d2)
	  double d1 = (exp(lm[a][k][0]*u)-1)/lm[a][k][0];
	  double d2 = (exp(lm[a][k][1]*u)-1)/lm[a][k][1];
	  C[k][0][0] = B[0][0]*d1*Bi[0][0] + B[0][1]*d2*Bi[1][0];
	  C[k][0][1] = B[0][0]*d1*Bi[0][1] + B[0][1]*d2*Bi[1][1];
	  C[k][1][0] = B[1][0]*d1*Bi[0][0] + B[1][1]*d2*Bi[1][0];
	  C[k][1][1] = B[1][0]*d1*Bi[0][1] + B[1][1]*d2*Bi[1][1];
	}
	if (q_first+qpoint-1 < 4*(b-1)) {
	  // C[k] = C[k] * exp(A[a] * .25)
	  vector<vector<double> > Ccopy = C[k];
	  vector<vector<double> > expA(2, vector<double>(2));
	  double d1 = exp(lm[a][k][0] * .25);
	  double d2 = exp(lm[a][k][1] * .25);
	  expA[0][0] = B[0][0]*d1*Bi[0][0] + B[0][1]*d2*Bi[1][0];
	  expA[0][1] = B[0][0]*d1*Bi[0][1] + B[0][1]*d2*Bi[1][1];
	  expA[1][0] = B[1][0]*d1*Bi[0][0] + B[1][1]*d2*Bi[1][0];
	  expA[1][1] = B[1][0]*d1*Bi[0][1] + B[1][1]*d2*Bi[1][1];
	  C[k][0][0] = Ccopy[0][0]*expA[0][0] + Ccopy[0][1]*expA[1][0];
	  C[k][0][1] = Ccopy[0][0]*expA[0][1] + Ccopy[0][1]*expA[1][1];
	  C[k][1][0] = Ccopy[1][0]*expA[0][0] + Ccopy[1][1]*expA[1][0];
	  C[k][1][1] = Ccopy[1][0]*expA[0][1] + Ccopy[1][1]*expA[1][1];
	}

	// compute conditional variance  Var( E(N_IR(t_{n+b})-N_IR(t_{n+b-1}) | t_{n+a+a'/4}) | t_{n+a+(a'-1)/4} )
	// if (t_{n+a+a'/4} > t_{n+b-1}), E(N_IR(t_{n+b})-N_IR(t_{n+b-1}) | t_{n+a+a'/4}) = c1 * E(t_{n+a+a'/4}) + c2 * I(t_{n+a+a'/4}) + N_IR(t_{n+a+a'/4}) - N_IR(t_{n+b-1})
	// where c1 = nu_IR*C[k][1][0] and c2 = nu_IR*C[k][1][1]
	// otherwise, the last term does not exist.
	// c1*E(t_{n+a+a'/4}) + c2*I(t_{n+a+a'/4}) + N_IR(t_{n+a+a'/4})-N_IR(t_{n+a+(a'-1)/4}) = c1*E(t_{n+a+(a'-1)/4})+c2*I(t_{n+a+(a'-1)/4}) + c1*(N_SE(t_{n+a+a'/4})-N_SE(t_{n+a+(a'-1)/4})) + (c2-c1)*(N_EI(t_{n+a+a'/4})-N_EI(t_{n+a+(a'-1)/4})) + (1-c2)*(N_IR(t_{n+a+a'/4})-N_IR(t_{n+a+(a'-1)/4}))
	// the conditional variance of the above quantity given the state at time t_{n+a+(a'-1)/4} equals  (c1^2 condVar(N_SE) + (c1-c2)^2 condVar(N_EI) + (1-c2)^2 condVar(N_IR))
	// if (t_{n+a+a'/4} > t_{n+b-1}), the term N_IR(t_{n+a+a'/4})-N_IR(t_{n+a+(a'-1)/4}) is not included.
	double time1 = t+(q_first+qpoint-2)/4.0;
	time1 = (time1 < t+s/double(im) ? t+s/double(im) : time1); // when qpoint=1, the starting point is t+s/im
	double time2 = t+(q_first+qpoint-1)/4.0;
	double tinc = time2-time1;
	int year1 = int(time1 * 7.0 / 365.25) + startyear; 	
	int year2 = int(time2 * 7.0 / 365.25) + startyear; 
	int day1 = int(time1 * 7.0 + tol - (year1-startyear)*365.25);
	int day2 = int(time2 * 7.0 + tol - (year2-startyear)*365.25);
	double meanInfectionRate = .5 * (proj_q[qpoint-1][5*k]*beta_s*hol_eff(day1, hol_coeff)*pow(proj_q[qpoint-1][5*k+2]/proj_q[qpoint-1][5*k+4],alpha) + proj_q[qpoint][5*k]*beta_s*hol_eff(day2, hol_coeff)*pow(proj_q[qpoint][5*k+2]/proj_q[qpoint][5*k+4],alpha));
	double varN_SE = tinc * meanInfectionRate * (1+sigma2*meanInfectionRate);
	double meanE = .5 * (proj_q[qpoint-1][5*k+1] + proj_q[qpoint][5*k+1]);
	double varN_EI = tinc * meanE*nu_EI * (1+sigma2*meanE*nu_EI);
	double meanI = .5 * (proj_q[qpoint-1][5*k+2] + proj_q[qpoint][5*k+2]);
	double varN_IR = tinc * meanI*nu_IR * (1+sigma2*meanI*nu_IR);
	double c1 = nu_IR*C[k][1][0];
	double c2 = nu_IR*C[k][1][1];
	bool inLastInt = (q_first+qpoint-1 > 4*(b-1)); // Is (t_{n+(q_first+qpoint-1)/4} > t_{n+b-1})?
	var[k] += (c1*c1*varN_SE + (c1-c2)*(c1-c2)*varN_EI + (-c2+inLastInt)*(-c2+inLastInt)*varN_IR);
      }
      
      // N_IRbar: projection (with mean process) of N_IR(max(t_{n,s}, t_{n+b-1})) - N_IR(t_{n+b-1}), b=lookaheads[n].
      // Y_{n+b} | X_{t_{n,s}} ~ rep * (NegBin(mean=dN_IRbar1, var=dN_IR) + dN_IR2) + sqrt{V_meas}*N(0,1),
      // decompose increments: dN_IRbar := dN_IRbar1 + dN_IR2
      // dN_IRbar1 := N_IRbar(t_{n+b}) - N_IRbar(max(t_{n+b-1}, t_{n,s})), dN_IR2 := N_IR(max(t_{n+b-1}, t_{n,s})) - N_IR(t_{n+b-1})
      double dN_IRbar = proj[pos][5*k+3];
      int dN_IR2 = (lookaheads[pos]==1 ? state[5*k+3] : 0);
      double dN_IRbar1 = dN_IRbar - dN_IR2;

      // V_meas :=rep(1-rep)dN_IRbar+rep^2repOD^2dN_IRbar^2+1.
      double var_meas = dN_IRbar*rep[cityindex[k]]*(1-rep[cityindex[k]]) + repOD*repOD*dN_IRbar*dN_IRbar*rep[cityindex[k]]*rep[cityindex[k]] + 1.0;
      // p_{Y_{n+b}|X_{t_{n,s}}} can be expressed as a convolution of a negative binomial distribution and a normal distribution,
      // roughly approximated as TERM1 + TERM2 where
      // TERM1 = sum_{i=-3:3} Prob[y_{n+b} | dN_IR = (y_{n+b}+i*sqrt(V_meas))/rep] * Prob[dN_IR = (y_{n+b}+i*sqrt(V_meas))/rep]*sqrt(V_meas)/rep
      // TERM2 = sum_{i=-3:3} Prob[y_{n+b} | dN_IR = dN_IRbar+i*sqrt(V_dN_IR)] * Prob[dN_IR = dN_IRbar+i*sqrt(V_dN_IR)]*sqrt(V_dN_IR)
      vector<int> dN_IR; // determine the points at which of dN_IR the convolution probability will be computed
      for (int i = -3; i<=3; i++) {
	dN_IR.push_back(int(.5+(meas[pos][k]+i*sqrt(var_meas))/rep[cityindex[k]]));
	dN_IR.push_back(int(.5+dN_IRbar+i*sqrt(var[k])));
      }
      quicksort(dN_IR, 0, 13);
      double pred_pmf = 0.0; // predictive pmf for y given the current state
      double p_prev = 0.0; // placeholder for previous probability term
      for (int i=0; i<14; i++) {
	if (i>0 && dN_IR[i] == dN_IR[i-1]) { continue; }
	double p_this = dnpmf(meas[pos][k], dN_IR[i]*rep[cityindex[k]], sqrt(var_meas)) * negbinpmf(dN_IR[i]-dN_IR2, dN_IRbar1, var[k], false);
	pred_pmf += p_this;
	if (i>0) {
	  pred_pmf += (dN_IR[i]-dN_IR[i-1]-1) * (p_this > p_prev ? p_prev : p_this); // add the probability computed for this dN_IR[i] for all the numbers between dN_IR[i-1]+1, ..., dN_IR[i]-1
	}
	p_prev = p_this;
      }
      dm[pos] += log(pred_pmf);
    }
    exp_dm[pos] = exp(dm[pos]);
  }

    /*      
    vector<vector<double> > proj(lasize, vector<double>(5*ncity)); // state simulated to integer lookahead points
    for (int pos=0; pos<lasize; pos++) 
      proj[pos] = proj_q[lookaheads[pos]*4-q_first+1];

      // compute approximate predictive likelihood
      dm[pos] += log(dnpmf(meas[pos][k], proj[pos][5*k+3] * rep[cityindex[k]], std::sqrt(proj[pos][5*k+3]*rep[cityindex[k]]*(1-rep[cityindex[k]]) + repOD*repOD*proj[pos][5*k+3]*proj[pos][5*k+3]*rep[cityindex[k]]*rep[cityindex[k]] + 1 + var[k])));
    }
    exp_dm[pos] = exp(dm[pos]);
    }
    */
      /*
      // for debugging purpose
      double var_dN_IR = (state[5*k+2]*nu_IR*sigma2 + 1.0)*state[5*k+2];
      for (int nn = 0; nn < lookaheads[n] ; nn++)
	var_dN_IR += ( (proj[nn][5*k+2]*nu_IR*sigma2+1.0)*proj[nn][5*k+2] );
      var_dN_IR *= ( nu_IR * (lookaheads[n]-(s+0.0)/im) / (lookaheads[n] + 1.0) );
      // variance of projected cumulative I to R transition, computed as (I * nu_IR * sigma2 + 1) * I * nu_IR * delta_t. Since the state changes while projecting, we take average of the variance estimated from the initial state and the intermediate projections (including the end point).
      */




  /*
  vector<double> state_double(5*ncity);
  for (size_t k=0; k<state.size(); k++)
    state_double[k] = double(state[k]);
  int maxla = lookaheads[lookaheads.size()-1];
  vector<vector<double> > projections(maxla, vector<double>(5*ncity)); // the state simulated through the end of time interval (t+lookahead)
  simulate_deterministic_process(t+(s+0.0)/im, t+1.0, state_double, projections[0], theta, 0.05);
  for (int c = 1; c < maxla; c++)
    simulate_deterministic_process(t+c+0.0, t+c+1.0, projections[c-1], projections[c], theta, 0.05);

  // new version
  for (size_t pos = 0; pos < lookaheads.size(); pos++) {
    for (int k = 0; k < ncity; k++) {
      // b:=lookaheads[pos], t_start:=max(t_{n+b-1},t_{n,s}), t_end:=t_{n+b}
      // variance of projected cumulative I to R transition, computed as INT (I * nu_IR * sigma2 + 1) * I * nu_IR * delta_t. Since the state changes while projecting, we take average of the variance estimated from the initial state and the intermediate projections (including the end point).
      // var(N_IR(t_end)-N_IR(t_start)) = E(Var(N_IR(t_end)-N_IR(t_start))|B_{t_end}) + Var(E(N_IR(t_start)-N_IR(t_end)|B_{t_end})) ~~ (Iavg*nu_IR*sigma2+1)*Iavg*nu_IR*(t_end-t_start) + Var(nu_IR*(t_end-t_start)*(I(t_{n,s})+N_EI(t_start)-N_EI(t_{n,s})-N_IR(t_start)+N_IR(t_{n,s})))
      // ~~ (sum of (Iavg*nu_IR*sigma2+1)*Iavg*nu_IR*interval_length for all intervals up to lookaheads[pos]) + (sum of (Eavg*nu_EI*sigma2+1)*Eavg*nu_EI*interval_length for all intervals up to lookaheads[pos]-1)
      double var_I_tstart = 0.0; // variance of I(t_start)
      for (int c = 0; c < lookaheads[pos]-1; c++) {
	double Iavg = .5*( (c==0 ? state_double[5*k+2] : projections[c-1][5*k+2]) + projections[c][5*k+2] ); //average of I(t) between the two endpoints of this interval
	double Eavg = .5*( (c==0 ? state_double[5*k+1] : projections[c-1][5*k+1]) + projections[c][5*k+1] ); //average of E(t) between the two endpoints of this interval
	var_I_tstart += (Iavg*nu_IR*sigma2+1)*Iavg*nu_IR*(c==0 ? 1.0-(s+0.0)/im : 1.0); // variance of dN_IR for this interval
	var_I_tstart += (Eavg*nu_EI*sigma2+1)*Eavg*nu_EI*(c==0 ? 1.0-(s+0.0)/im : 1.0); // variance of dN_EI for this interval except the last one
      }
      double Iavg = .5*( (lookaheads[pos]==1 ? state_double[5*k+2] : projections[lookaheads[pos]-2][5*k+2]) + projections[lookaheads[pos]-1][5*k+2] ); //average of I(t) between t_start and t_end
      double var_dN_IR = var_I_tstart * nu_IR*nu_IR * 1.0*1.0  +  (Iavg*nu_IR*sigma2+1)*Iavg*nu_IR*(lookaheads[pos]==1 ? 1.0-(s+0.0)/im : 1.0); // // Var(E(dN_IR|B_{t_start})) + E(Var(dN_IR|B_{t_start}))

      // N_IRbar: projection with mean process
      // decompose increments: dN_IRbar := dN_IRbar1 + dN_IR2
      // dN_IRbar1 := N_IRbar(t_{n+b}) - N_IRbar(max(t_{n+b-1}, t_{n,s})), dN_IR2 := N_IR(max(t_{n+b-1}, t_{n,s})) - N_IR(t_{n+b-1})
      double dN_IRbar = projections[pos][5*k+3];
      int dN_IR2 = (lookaheads[pos]==1 ? state[5*k+3] : 0);
      double dN_IRbar1 = dN_IRbar - dN_IR2;
      // N_IR(max(t_{n,s}, t_{n+b-1})) - N_IR(t_{n+b-1}), b=lookaheads[pos].      // Y_{n+b} | X_{t_{n,s}} ~ rep * (NegBin(mean=dN_IRbar1, var=dN_IR) + dN_IR2) + sqrt{V_meas}*N(0,1),
      // V_meas :=rep(1-rep)dN_IRbar+rep^2repOD^2dN_IRbar^2+1.
      double var_meas = dN_IRbar*rep[cityindex[k]]*(1-rep[cityindex[k]]) + repOD*repOD*dN_IRbar*dN_IRbar*rep[cityindex[k]]*rep[cityindex[k]] + 1.0;
      // p_{Y_{n+b}|X_{t_{n,s}}} can be expressed as a convolution of a negative binomial distribution and a normal distribution,
      // roughly approximated as TERM1 + TERM2 where
      // TERM1 = sum_{i=-3:3} Prob[y_{n+b} | dN_IR = (y_{n+b}+i*sqrt(V_meas))/rep * Prob[dN_IR = (y_{n+b}+i*sqrt(V_meas))/rep]*sqrt(V_meas)/rep
      // TERM2 = sum_{i=-3:3} Prob[y_{n+b} | dN_IR = dN_IRbar+i*sqrt(V_dN_IR)] * Prob[dN_IR = dN_IRbar+i*sqrt(V_dN_IR)]*sqrt(V_dN_IR)
      vector<int> dN_IR; // determine the points at which of dN_IR the convolution probability will be computed
      for (int i = -3; i<=3; i++) {
	dN_IR.push_back(int(.5+(meas[pos][k]+i*sqrt(var_meas))/rep[cityindex[k]]));
	dN_IR.push_back(int(.5+dN_IRbar+i*sqrt(var_dN_IR)));
      }
      quicksort(dN_IR, 0, 13);
      double pred_pmf = 0.0; // predictive pmf for y given the current state
      double p_prev = 0.0; // placeholder for previous probability term
      for (int i=0; i<14; i++) {
	if (i>0 && dN_IR[i] == dN_IR[i-1]) { continue; }
	double p_this = dnpmf(meas[pos][k], dN_IR[i]*rep[cityindex[k]], sqrt(var_meas)) * negbinpmf(dN_IR[i]-dN_IR2, dN_IRbar1, var_dN_IR, false);
	pred_pmf += p_this;
	if (i>0) {
	  pred_pmf += (dN_IR[i]-dN_IR[i-1]-1) * (p_this > p_prev ? p_prev : p_this); // add the probability computed for this dN_IR[i] for all the numbers between dN_IR[i-1]+1, ..., dN_IR[i]-1
	}
	p_prev = p_this;
      }
      dm[pos] += log(pred_pmf);
    }
    exp_dm[pos] = exp(dm[pos]);
  }
  */

  /* old version  
  for (size_t pos = 0; pos < lookaheads.size(); pos++) {
    for (int k = 0; k < ncity; k++) {
      double var_dN_IR = (state[5*k+2]*nu_IR*sigma2 + 1.0)*state[5*k+2];
      for (int nn = 0; nn < lookaheads[pos]; nn++)
	var_dN_IR += ( (projections[nn][5*k+2]*nu_IR*sigma2+1.0)*projections[nn][5*k+2] );
      var_dN_IR *= ( nu_IR * (lookaheads[pos]-(s+0.0)/im) / (lookaheads[pos] + 1.0) );
      // variance of projected cumulative I to R transition, computed as (I * nu_IR * sigma2 + 1) * I * nu_IR * delta_t. Since the state changes while projecting, we take average of the variance estimated from the initial state and the intermediate projections (including the end point).
      double pred_pmf = log(dnpmf(meas[pos][k], projections[pos][5*k+3] * rep[cityindex[k]], std::sqrt(projections[pos][5*k+3]*rep[cityindex[k]]*(1-rep[cityindex[k]]) + repOD*repOD*projections[pos][5*k+3]*projections[pos][5*k+3]*rep[cityindex[k]]*rep[cityindex[k]] + 1 + var_dN_IR)));
      dm[pos] += pred_pmf;
    }
    exp_dm[pos] = exp(dm[pos]);
  }
  */
  
  return (logarithm ? dm : exp_dm);
  
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
