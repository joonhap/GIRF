#include <fstream>
#include <algorithm>
#include <omp.h>
#include "stf_parallel.h"

using std::vector;
using std::cout;
using std::endl;

template <typename T>
void quicksort(vector<T>& v, int start, int end) {
  int i=start;
  int j=end;
  T pivot = v[(start+end)/2];

  while (i<=j) {
    while (v[i]<pivot)
      i++;
    while (v[j]>pivot)
      j--;
    if (i<=j) {
      T tmp = v[i];
      v[i] = v[j];
      v[j] = tmp;
      i++;
      j--;
    }
  }

  if (start < j) { quicksort(v, start, j); }
  if (i < end) { quicksort(v, i, end); }
}


template <typename state_T, typename obs_T>
void stf<state_T, obs_T>::simulate(std::vector<state_T> state0, int initial_seed) {
  // simulate states and observations, when the underlying states take integer values

  std::mt19937 rng(initial_seed);

  std::vector<state_T> oldstate(nvar);
  for (int k = 0; k < nvar; k++)
    oldstate[k] = state0[k];

  hidden_state[0] = state0; // save true hidden_state

  std::vector<state_T> state(nvar);

  for (int t = 0; t != ntimes; t++) {
    rprocess(t, oldstate, state, theta, rng);
    oldstate = state;

    hidden_state[t+1] = state;

    rmeasure(t+1, state, obs[t], theta, rng);
  }
}



// // // intermediary filter // // //
// iterated filtering (IF2)
template <typename state_T, typename obs_T>
  void stf<state_T, obs_T>::imfilter(const std::vector<std::vector<std::vector<state_T> > >& initstates, int im, int R, int J, std::vector<std::vector<double> >& lestimate, std::vector<std::vector<double> >& state_mean, std::vector<std::vector<state_T> >& q10, std::vector<std::vector<state_T> >& q50, std::vector<std::vector<state_T> >& q90, double (*f1)(double), std::vector<std::vector<double> >& f1_mean, std::vector<std::vector<std::vector<double> > >& ESS, const std::vector<std::vector<std::vector<double> > >& initThetaSwarm, std::vector<double> init_perturb_size, int M, std::function<double (int)> coolschedule, std::vector<std::vector<std::vector<double> > >& thetaTrajectory, std::vector<std::vector<std::vector<double> > >& finalThetaSwarm, bool weighted_interaction, int max_lookahead, std::function<bool (int, int)> island_interact_condition, int initial_seed) {
  // initstates: initial states; im: the number of subdivisions of a time interval; R: the number of repititions of intermediary filtering at each time interval; J: the number of particles used in each intermediary filtering; lestimate: log likelihood esitmate for each time interval; state_mean: state estimate; q10, q50, q90: quantiles of state variable, f1: function of state to be estimated; f1_mean: estimate of f1 of state; ESS: effective sample size at each resampling stage; initThetaSwarm: initial parameter vectors for IF2; init_perturb_size: initial parameter perturbation size for IF2; M: number of IF2 iterations; coolschedule: function specifying how fast perturbation size decreases; thetaTrajectory: the trajectory of the average of parameter estimates over IF2; finalThetaSwarm: final parameter locations (for subsequently continued IF2 iterations), weighted_interaction: if true, duplicate or discard islands according to their weights (cf. del moral et al's island particle filter.); max_lookahead: maximum number of observations to look ahead for the computation of the guide function; initial_seed: seed value for rng initialization; island_interact_condition: a function returning true/false whether the island interaction should occur at the current time point

  std::vector<std::vector<std::vector<state_T> > > prev; // 'prev' are the collection of all monte carlo samples at each time step (R*J*K)
  std::vector<std::vector<std::vector<double> > > prev_theta = initThetaSwarm; // collection of parameter vectors. (R*J*param_dim)
  std::vector<std::vector<double> > oldltilde_storage(R); // global storage of oldltilde (see below)
  vector<vector<vector<vector<double> > > > forecast_var(R, vector<vector<vector<double> > >(J)); // global storage of forecast variances

  lestimate.resize(M);
  thetaTrajectory.resize(M);

  state_mean.resize(ntimes);
  f1_mean.resize(ntimes);
  q10.resize(ntimes);
  q50.resize(ntimes);
  q90.resize(ntimes);

  ESS.resize(R);

  for (int m=0; m<M; m++) {
    std::uniform_real_distribution<double> uniform_global(0.0, 1.0);
    std::normal_distribution<double> normal_global(0.0, 1.0);
    std::mt19937 rng_global(initial_seed + SEED_INC_PER_IF_ITER*m);

    thetaTrajectory[m].resize(ntimes);

    std::vector<double> perturb_size(theta.size());
    // change the perturbation size for parameter vector
    for(std::size_t param_dim = 0; param_dim < theta.size(); param_dim++) 
      perturb_size[param_dim] = coolschedule(m) * init_perturb_size[param_dim];

    // perturb initial value parameters
    for (int i=0; i<R; i++)
      for (int j=0; j<J; j++)
	for(std::size_t ivp_dim = 0; ivp_dim < pos_ivp.size(); ivp_dim++)
	  prev_theta[i][j][pos_ivp[ivp_dim]] += perturb_size[pos_ivp[ivp_dim]] * normal_global(rng_global);

    set_initstates(prev, initstates, prev_theta);

    std::vector<double> island_logweights(R, 0.0);

#pragma omp parallel num_threads(R)
    {
      int i = omp_get_thread_num();

      std::mt19937 rng(1299721*i + SEED_INC_PER_IF_ITER*m + initial_seed + 25828);
      std::normal_distribution<double> normal_loc(0.0, 1.0);
      std::uniform_real_distribution<double> uniform_loc(0.0, 1.0);

      std::vector<std::vector<state_T> > prev_loc = prev[i]; // 'prev_loc' are the local particles
      std::vector<std::vector<double> > prev_theta_loc = prev_theta[i]; // local parameter vectors
      std::vector<double> oldltilde(J, 0.0); // local estimate of log likelihood of next y given previous (intermediary) x under auxiliary measurement distribution
	
      std::vector<std::vector<double> > ESS_loc(ntimes, std::vector<double>(im)); // local effective sample sizes
      
      double cumulative_l = 0.0; // local cumulative sum of intermediate log likelihood until island interaction

      for (int t=0; t<ntimes; t++) {
      // compute the pseudo-predictive likelihood from the last step, i.e., the denominator of weight for s=0.
	if (t > 0) {
          for (int j=0; j<J; j++) 
            oldltilde[j] -= dmeasure(t, prev_loc[j], obs[t-1], prev_theta_loc[j], true);
        }

        // lookahead steps and the corresponding observations for computing the guide function (imdmeasure)
        vector<int> lookaheads;
        vector<vector<obs_T> > observations;
        for (int lookahead = 1; lookahead <= std::min(max_lookahead, ntimes-t); lookahead++) {
          lookaheads.push_back(lookahead);
          observations.push_back(obs[t+lookahead-1]);
        }

        vector<vector<vector<double> > > forecast_var_loc(J, vector<vector<double> >(lookaheads.size())); // estimated variability in the time-forward projections (forecasts). Used for computing the guide function (imdmeasure)

        for (int s=0; s<im; s++) {
          std::vector<std::vector<state_T> > x(J, std::vector<state_T>(nvar)); // proposed particles at this intermediary step
	  std::vector<std::vector<double> > theta_loc(J, std::vector<double>(theta.size())); // perturbed parameter vectors.
	  std::vector<double> ltilde(J, 0.0); // estimate of log likelihood of next y given current (intermediary) x under auxiliary measurement distribution
          std::vector<double> logw(J, 0.0); // weights corresponding to the pairs of proposed intermediary particles and their parents.

	  for(int j=0; j<J; j++) {
	    theta_loc[j] = prev_theta_loc[j];
	    for(std::size_t param_dim = 0; param_dim < theta.size()-pos_ivp.size(); param_dim++)
	      theta_loc[j][param_dim] += perturb_size[param_dim] / sqrt(im) * normal_loc(rng);

	    imrprocess(t, s, im, prev_loc[j], x[j], theta_loc[j], rng);

	    // evaluate fitness to the next few observations
            vector<double> pseudo_predictive_loglikelihoods = imdmeasure(t, s+1, lookaheads, im, x[j], observations, theta_loc[j], forecast_var_loc[j], true);

            for (size_t pos = 0; pos < lookaheads.size(); pos++) {
	      if (annealed) { // if annealed predictive likelihood is used
                ltilde[j] += (pseudo_predictive_loglikelihoods[pos] * (1.0 - (lookaheads[pos]-(s+1.0)/im) / std::max(2, std::min(max_lookahead, t+lookaheads[pos]))));
              } else { // if non-annealed predictive likelihood is used
                ltilde[j] += pseudo_predictive_loglikelihoods[pos];
	      }
            }
	    logw[j] = ltilde[j] - oldltilde[j];
	  }

          // if at least one particle has ltilde that is not -Inf, update cumulative_l and resample.
	  // otherwise, keep oldltilde and just move on.
	  bool any_survival = false;
	  for (int j=0; j<J; j++) {
	    if (ltilde[j] > -1.0/0.0) {
	      any_survival = true;
	      break;
	    }
	  }

	  double effSS; // effective sample size
	  
          if (any_survival) {
	    // if any of the parent particles have zero weights, resample with ltilde directly.
	    bool any_parent_zero_weight = false;
	    for (int j=0; j<J; j++) 
	      if (oldltilde[j] == -1.0/0.0) {
	        any_parent_zero_weight = true;
		break;
              }
	    if (any_parent_zero_weight) logw = ltilde;

  	    double avg_logw = 0.0; // an offset to compute sum_w (in case w is less than machine precision)
            int no_finite_logw = 0; // number of finite logw
	    for (int j=0; j<J; j++)
	      if (std::isfinite(logw[j])) {
	        avg_logw += logw[j];
                no_finite_logw++;
              }
	    avg_logw /= no_finite_logw;
	  
            vector<double> scaled_w(J); // w divided by the (exponentiated) offset
	    for (int j=0; j<J; j++)
	      scaled_w[j] = exp(logw[j] - avg_logw);
            double sum_scaled_w = 0.0;
	    for (int j=0; j<J; j++)
	      sum_scaled_w += scaled_w[j];
            double logsum_w = avg_logw + log(sum_scaled_w);

            // compute effective sample size
            double ssw = 0.0;
	    for (int j=0; j<J; j++)
	      ssw += scaled_w[j]*scaled_w[j];
	    effSS = sum_scaled_w*sum_scaled_w / ssw; // effective sample size

  	    // normalize weights
	    for(int j=0; j<J; j++)
	      logw[j] -= logsum_w;

            // resampling
	    vector<vector<vector<double> > > forecast_var_loc_storage = forecast_var_loc;
 	    char resampling_mode = 's'; // 's': systematic resampling, 'm': multinomial resampling
	    if (resampling_mode == 's') {
	      double p1 = uniform_loc(rng) / J;
	      double p2 = 0.0;
	      int j1 = -1;
	      for(int j=0; j<J; j++) {
	        while(p2 < p1 && j1<J-1)
	          p2 += exp(logw[++j1]);
	        prev_loc[j] = x[j1];
	        prev_theta_loc[j] = theta_loc[j1];
	        oldltilde[j] = ltilde[j1];
		forecast_var_loc[j] = forecast_var_loc_storage[j1];
	        p1 += 1.0 / J;
	      }
	    } else if (resampling_mode == 'm') {
	      for(int j=0; j<J; j++) {
	        double p1 = uniform_loc(rng);
	        double p2 = 0.0;
	        int j1 = -1;
	        while (p2 < p1 && j1<J-1)
 	          p2 += exp(logw[++j1]);
	        prev_loc[j] = x[j1];
	        prev_theta_loc[j] = theta_loc[j1];
	        oldltilde[j] = ltilde[j1];
		forecast_var_loc[j] = forecast_var_loc_storage[j1];
	      }
            }

            // update likelihood estimate
            cumulative_l += (any_parent_zero_weight ? 1.0/0.0 : (logsum_w - log(J))); 

          } else { // when all particles have zero weight, move along with existing particles.
	    cout << " All particle have zero weight." << "  t: " << t << " s: " << s << " particle island: " << i << endl;
	    prev_loc = x;
	    prev_theta_loc = theta_loc;
	    effSS = 0.0/0.0;
	  }

	  if (m == (M-1)) { ESS_loc[t][s] = effSS; }

	  //double sdw = sqrt((1.0/effSS - 1.0/J)/(J-1));
	  //std::cout << "i : " << i << "  t : " << t << ",  s : " << s << ",   ESS : " << effSS << ",   CV (particle weight) : " << (sdw*J) << std::endl; // ESS=1/sum_sq, CV=sd/mean

	  // island interaction
	  if (island_interact_condition(t,s) || (t==(ntimes-1)&&s==(im-1))) {
	    prev[i] = prev_loc; // copy local particles into the central storage
	    prev_theta[i] = prev_theta_loc; // copy local parameter vectors into the central storage
            if (!any_survival) oldltilde = vector<double>(J, -1.0/0.0);
	    oldltilde_storage[i] = oldltilde; // copy local oldltilde into the central storage
	    forecast_var[i] = forecast_var_loc; // copy forecast variance into the central storage
  	    island_logweights[i] = (any_survival ? cumulative_l : -1.0/0.0); // copy local cumulative log likelihood estimate to the central storage

#pragma omp barrier
#pragma omp single
	    {
              double avg_island_logweights = 0.0; // an offset to compute sum_island_weights (in case island_weights are less than machine precision)
	      int no_finite_island_logweight = 0; // number of finite island logweight
	      for (int r=0; r<R; r++)
	      	if (std::isfinite(island_logweights[r])) {
	          avg_island_logweights += island_logweights[r];
                  no_finite_island_logweight++;
		}
              avg_island_logweights /= no_finite_island_logweight;

	      vector<double> scaled_island_weights(R);
	      for (int r=0; r<R; r++)
	        scaled_island_weights[r] = (island_logweights[r] == 1.0/0.0 ? 0.0 : exp(island_logweights[r] - avg_island_logweights)); // if island_logweight equals Inf, give the island the weight equal to avg_island_logweights.
              double logsum_island_weights = avg_island_logweights + log(accumulate(scaled_island_weights.begin(), scaled_island_weights.end(), 0.0));

              lestimate[m].push_back(logsum_island_weights - log(R)); // update likelihood estimate
	      
	      std::vector<int> pid(R*J); // parent particle ID
	      for (int n=0; n<R*J; n++)
		pid[n] = n;
	      
	      bool shuffle = false; // if true, shuffle all particles between islands
	  
	      if (weighted_interaction && (logsum_island_weights > -1.0/0.0)) {
                // normalize island weights
		for (int r=0; r<R; r++)
		  island_logweights[r] -= logsum_island_weights;
		// resample islands
		std::vector<int> id_before_resample(R);
		double p1 = uniform_global(rng_global) / R;
		double p2 = 0.0;
		int r1 = -1;
		for(int r=0; r<R; r++) {
		  while(p2 < p1 && r1<R-1)
		    p2 += exp(island_logweights[++r1]);
		  id_before_resample[r] = r1;
		  p1 += 1.0 / R;
		}
		for (int r=0; r<R; r++)
		  for (int j=0; j<J; j++)
		    pid[r*J+j] = id_before_resample[r]*J+j;

		//double i_ssw = 0.0;
		//for (int r=0; r<R; r++)
		//  i_ssw += exp(2*island_logweights[r]);
		//double i_sdw = sqrt((i_ssw - 1.0/R)/(R-1));
		//std::cout << "[ISLAND] t : " << t << "  s : " << s << ",   ESS : " << 1.0/i_ssw << ",   CV (particle weight) : " << (i_sdw*R) << std::endl;
	      }
	      if (shuffle) {   
		// shuffle particles between islands
		for (int n=0; n<R*J; n++) {
		  int random_position = int(uniform_global(rng_global)*(R*J-n));
		  int replacement = pid[random_position];
		  pid[random_position] = pid[R*J-n-1];
		  pid[R*J-n-1] = replacement;
		}
	      }

	      // update the central particle storage
	      std::vector<std::vector<std::vector<state_T> > > prev_copy = prev;
	      std::vector<std::vector<std::vector<double> > > prev_theta_copy = prev_theta;
	      std::vector<std::vector<double> > oldltilde_storage_copy = oldltilde_storage;
	      vector<vector<vector<vector<double> > > > forecast_var_copy = forecast_var;
	      for (int r=0; r<R; r++)
		for (int j=0; j<J; j++) {
		  int id = pid[r*J+j];
		  prev[r][j] = prev_copy[id/J][id%J];
		  prev_theta[r][j] = prev_theta_copy[id/J][id%J];
		  oldltilde_storage[r][j] = oldltilde_storage_copy[id/J][id%J];
		  forecast_var[r][j] = forecast_var_copy[id/J][id%J];
		}

              // compute average theta value
	      if (s == im-1) {
		std::vector<double> avg_theta(theta.size(), 0.0);
		for (int r=0; r<R; r++)
		  for (int j=0; j<J; j++)
	            for (size_t param_dim=0; param_dim<theta.size(); param_dim++)
	              avg_theta[param_dim] += prev_theta[r][j][param_dim];
	        for (size_t param_dim=0; param_dim<theta.size(); param_dim++)
	          avg_theta[param_dim] /= double(R*J);
		thetaTrajectory[m][t] = avg_theta;
              }

              // end the iteration if there are too many nan or -inf's.
              int no_fin = 0;
              for (int tt=0; tt<=t; tt++)
                no_fin += (std::isfinite(lestimate[m][tt]));
              if (no_fin < t*0.75 && t>20) {
                cout << "All particles were lost too many times in this iteration up to time " << t << "." << endl;
                exit(1);
              }
            } // end of pragma omp single

	    // copy into local storages
	    prev_loc = prev[i];
	    prev_theta_loc = prev_theta[i];
	    oldltilde = oldltilde_storage[i];
	    forecast_var_loc = forecast_var[i];
	    
            // reset local cumulative likelihood estimates
	    cumulative_l = 0.0;
	    
	  } // end of if (interaction condition)

	} // end of for (s) loop
	
// compute averages
	if (m == (M-1)) { // compute state average only at the last IF2 iteration
#pragma omp barrier
#pragma omp single
	  {
	    // compute particle average
	    state_mean[t] = std::vector<double>(nvar, 0.0);
	    f1_mean[t] = std::vector<double>(nvar, 0.0);
	    for (int r=0; r<R; r++)
	      for (int j=0; j<J; j++){
		for (int k=0; k<nvar; k++) {
		  state_mean[t][k] += prev[r][j][k];
		  f1_mean[t][k] += (*f1)(prev[r][j][k]);
		}
	      }
	    for (int k=0; k<nvar; k++) {
	      state_mean[t][k] /= double(R*J);
	      f1_mean[t][k] /= double(R*J);
	    }

	    // compute particle quantiles
	    vector<vector<state_T> > prev_transpose(nvar, vector<state_T>(R*J));
	    for (int r=0; r<R; r++)
	      for (int j=0; j<J; j++)
		for (int k=0; k<nvar; k++)
		  prev_transpose[k][r*J+j] = prev[r][j][k];

	    for (int k=0; k<nvar; k++) {
	      quicksort(prev_transpose[k], 0, R*J-1);
	      q10[t].resize(nvar);
	      q50[t].resize(nvar);
	      q90[t].resize(nvar);	  
	      q10[t][k] = prev_transpose[k][ceil(R*J*.10)-1];
	      q90[t][k] = prev_transpose[k][ceil(R*J*.90)-1];
	      q50[t][k] = prev_transpose[k][ceil(R*J*.50)-1];
	    }
	  } // end of omp single
	} // end of if (m==(M-1))
      } // end of for (t) loop

      ESS[i] = ESS_loc;
    } // end of omp parallel
    
  } // end of for (m) loop

  finalThetaSwarm = prev_theta;
}



// non-interated filtering
template <typename state_T, typename obs_T>
  void stf<state_T, obs_T>::imfilter(const std::vector<std::vector<std::vector<state_T> > >& initstates, int im, int R, int J, std::vector<std::vector<double> >& lestimate, std::vector<std::vector<double> >& state_mean, std::vector<std::vector<state_T> >& q10, std::vector<std::vector<state_T> >& q50, std::vector<std::vector<state_T> >& q90, double (*f1)(double), std::vector<std::vector<double> >& f1_mean, std::vector<std::vector<std::vector<double> > >& ESS, std::vector<double> theta_point, bool weighted_interaction, int max_lookahead, std::function<bool (int, int)> island_interact_condition, int initial_seed) {
  // initstates: initial states; im: the number of subdivisions of a time interval; R: the number of repititions of intermediary filtering at each time interval; J: the number of particles used in each intermediary filtering

  std::vector<std::vector<std::vector<double> > > initThetaSwarm(R, std::vector<std::vector<double> >(J, theta_point));
  std::vector<double> perturb_size(theta.size(), 0.0);
  std::vector<std::vector<std::vector<double> > > thetaTrajectory;
  std::function<double (int)> coolschedule = [] (int m) -> double { return 1.0; };
  std::vector<std::vector<std::vector<double> > > finalThetaSwarm;

  imfilter(initstates, im, R, J, lestimate, state_mean, q10, q50, q90, f1, f1_mean, ESS, initThetaSwarm, perturb_size, 1, coolschedule, thetaTrajectory, finalThetaSwarm, weighted_interaction, max_lookahead, island_interact_condition, initial_seed);
}
 




// // // plain particle filter // // //
template <typename state_T, typename obs_T>
  void stf<state_T, obs_T>::ppfilter(const std::vector<std::vector<state_T> >& initstates, int J, std::vector<std::vector<double> > stateTrajectory, int initial_seed) {
  // initstate: initial states; J: number of particles

  std::mt19937 rng(initial_seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  std::vector<std::vector<state_T> > prev = initstates;
  std::vector<std::vector<state_T> > xP(J, std::vector<state_T>(nvar));
  std::vector<std::vector<state_T> > xF(J, std::vector<state_T>(nvar));

  stateTrajectory.resize(ntimes);

  for(int t=0; t<ntimes; t++) {
    std::vector<double> w(J, 1.0);

    for(int j=0; j<J; j++) { 
      rprocess(t, prev[j], xP[j], theta, rng);
      w[j] = dmeasure(t, xP[j], obs[t], theta, false);
    }

    double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
    for(int j=0; j<J; j++)
      w[j] /= sum_w;

    double p1 = uniform(rng) / J;
    double p2 = 0.0;
    int j1 = -1;
    for(int j=0; j<J; j++) {
      while(p2 < p1 && j1<J-1)
	p2 += w[++j1];
      xF[j] = xP[j1];
      p1 += 1.0/J;
    }
    
    prev = xF;

    // For evaluation: Coefficient of Variation of weights (w[j])
    double ssw = 0.0;
    for (int j=0; j<J; j++)
      ssw += w[j] * w[j];
    double sdw = sqrt((ssw - 1/J)/(J-1));
    std::cout << "time " << t << ",   ESS : " << 1/ssw << ",   CV (particle weight) : " << (sdw*J) << std::endl;
    // ^ ESS = 1 / sum of squares, CV = standard deviation / mean
    
    // take average of xF and then save
    std::vector<double> meanXF(nvar, 0.0);
    for (int j=0; j<J; j++)
      for (int k=0; k<nvar; k++)
	meanXF[k] += xF[j][k];
    for (int k=0; k<nvar; k++)
      meanXF[k] /= J;

    stateTrajectory[t] = meanXF;

  }
}
