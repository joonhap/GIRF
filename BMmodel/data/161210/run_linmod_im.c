#include <fstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cmath>
#include "../../stf_parallel.h"
#include "../../linmod.h"

using std::vector;
using std::cout; 
using std::endl;

double square(double x) {
  return x*x;
}

void cholesky(const vector<vector<double> >& A, vector<vector<double> >& L) {
  int n = A.size();
  L = vector<vector<double> >(n, vector<double>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      double s = 0;
      for (int k = 0; k < j; k++)
	s += L[i][k] * L[j][k];
      L[i][j] = (i == j) ?
	sqrt(A[i][i] - s) :
	(1.0 / L[j][j] * (A[i][j] - s));
    }
}

vector<vector<vector<double> > > Lmat;

int main(int argc, char** args)
{

  ////////////////////////////////////////
  int nvar = std::stoi(args[1]);
  stf<double, double> sobj(nvar,nvar,50); // state_dim, measurement_dim, ntimes

  vector<double> state0(sobj.nvar);
  // initial state
  for (int k = 0; k < sobj.nvar; k++) 
    state0[k] = 0.0; // initial state

  double alpha = std::stof(args[2]);
  double mytheta[] = {log(1.0), log(1.0), log(alpha/(1.0-alpha))}; // {s, d, alpha}
  // ^ s: process sd, d: measurement sd, alpha: interaction strength
  sobj.theta.resize(3);
  sobj.theta.assign(mytheta, mytheta+3);
  
  /////////////////////////////////////////
  // set write or read
  char write_or_read = 'r'; // read ('r') or write ('w')
  
  if (write_or_read == 'w') {
    int simulation_seed = 91414 + 981427*sobj.nvar + 2358*int(8158*alpha); // random seed for data generation
    sobj.simulate(state0, simulation_seed);
    
    // write obs into a file
    std::ofstream ofile;
    std::string ofilename = "../obs_K" + std::to_string(sobj.nvar) + "alpha" + std::to_string(alpha).substr(0,4) + ".txt";
    ofile.open(ofilename);
    for (int t = 0; t < sobj.ntimes; t++) {
      for (int k = 0; k < sobj.nvar; k++)
	ofile << sobj.obs[t][k] << " ";
      ofile << endl;
    }
    ofile.close();
    
    // write states into a file
    std::ofstream sfile;
    std::string sfilename = "../state_K" + std::to_string(sobj.nvar) + "alpha" + std::to_string(alpha).substr(0,4) + ".txt";
    sfile.open(sfilename);
    for (int t = 0; t < sobj.ntimes+1 ; t++) {
      for (int k = 0; k < sobj.nvar; k++)
	sfile << sobj.hidden_state[t][k] << " ";
      sfile << endl;
    }
    sfile.close();
    cout << "simulation finished." << endl;
    
  } else if (write_or_read == 'r') {
    
    // read from a file
    double read_double = 0.0;
    std::string ifilename = "../obs_K" + std::to_string(sobj.nvar) + "alpha" + std::to_string(alpha).substr(0,4) + ".txt";
    std::ifstream ifile(ifilename);
    if (!ifile.is_open()) {
      cout << "Data file does not exist." << endl;
      exit(1);
    }
    for (int t = 0; t < sobj.ntimes; t++)
      for (int k = 0; k < sobj.nvar; k++) {
	ifile >> read_double;
	sobj.obs[t][k] = read_double;
      }
  }

  // // // // // // // //
  // space-time filter //
  // // // // // // // //
  
  ///////////////////////////////////////////
  int R = 60;
  int J = 1000;
    
  int N = sobj.nvar;
  
  std::string save_directory = "";

  // storage matrix
  vector<vector<double> > state_mean;
  vector<vector<double> > q10;
  vector<vector<double> > q50;
  vector<vector<double> > q90;
  vector<vector<double> > f1_mean;
  vector<vector<double> > lestimate;
  vector<vector<vector<double> > > ESS;
  vector<vector<vector<double> > > thetaTrajectory;

  vector<vector<vector<double> > > initstates(R, vector<vector<double> >(J, state0));
    
  ///////////////////////////////////////////
  // filenames
  bool run_IF2 = false;
  bool weighted_interaction = true; // should properly weighted island interaction be used?
  int max_lookahead = 2; // max number of lookahead steps for intermediate model
  if (argc==1) {
    cout << "no exp no provided." << endl;
    exit(1);
  }
  int exp_no = std::stoi(args[3]);
  std::string common_string = "_K" + std::to_string(sobj.nvar) + "alpha" + std::to_string(alpha).substr(0,4) + "R" + std::to_string(R) + "J" + std::to_string(J) + "N" + std::to_string(N) + "T" + std::to_string(sobj.ntimes) + "if_" + (run_IF2 ? "T" : "F") + "wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + ".txt";

  Lmat.resize(max_lookahead*N+1);
  // compute the cholesky factor of covariance matrices
  // in this case we are not going to run mif, and do just particle filter with known parameter value.
  for (int n = 0; n <= max_lookahead*N; n++) {
    vector<vector<double> > V(nvar, vector<double>(nvar)); // covariance matrix
    double diagonal = (1.0/(N+0.0)*n)*exp(2*sobj.theta[0]) + exp(2*sobj.theta[1]);
    double off_diag = (1.0/(N+0.0)*n)*exp(2*sobj.theta[0])*(1.0/(1.0+exp(-sobj.theta[2])));
    for (int k=0; k<nvar; k++)
      for (int l=0; l<nvar; l++)
	V[k][l] = (k == l) ? diagonal : off_diag;
    vector<vector<double> > L;
    cholesky(V, L);
    Lmat[n] = L;
    }

  //////////////////////////////////////////
  int M;

  std::clock_t stf_start = clock();

  if (run_IF2) {
    vector<vector<vector<double> > > initThetaSwarm(R, vector<vector<double> >(J, sobj.theta));
    ////////////////////////////////////
    bool continued = false; // indicator of whether IF2 is continued from a previous one
    int start_m; // starting IF2 iteration number (if a continued IF2, m is non-zero)

    if (!continued) {
      start_m = 0;
      for (int i=0; i<R; i++)
	for (int j=0; j<J; j++) {
	  ///////////////////////////////////////////
	  std::mt19937 rng(852316);
	  std::normal_distribution<double> normal(0.0, 1.0);
	  initThetaSwarm[i][j][0] = sobj.theta[0] + log(2.0) + normal(rng) * 0.05;
	  initThetaSwarm[i][j][1] = sobj.theta[1] + log(2.0) + normal(rng) * 0.05;
	  initThetaSwarm[i][j][2] = sobj.theta[2] + log(2.0) + normal(rng) * 0.05;
	}
    }
    else if (continued) {
      ////////////////////////////////////////
      start_m = 10; // start from iteration number start_m

      double read_double;
      std::string iswarmfilename = save_directory + "thetaswarm_m" + std::to_string(start_m) + common_string;
      std::ifstream iswarmfile(iswarmfilename);
      if (!iswarmfile.is_open()) {
	cout << "The theta swarm file does not exist." << endl;
	exit(1);
      }

      for (int i=0; i<R; i++)
	for (int j=0; j<J; j++)
	  for (std::size_t theta_dim=0; theta_dim < sobj.theta.size(); theta_dim++) {
	    iswarmfile >> read_double;
	    initThetaSwarm[i][j][theta_dim] = read_double;
	  }
    }
    
    vector<double> perturb_size(sobj.theta.size(), 0.0);
    perturb_size[0] = 0.01;
    perturb_size[1] = 0.01;
    perturb_size[2] = 0.01;

    /////////////////////////////////
    M = 1; // how many iterations of mif will take place?
    std::function<double (int)> coolschedule = [start_m] (int m) -> double { return pow(0.95, m + start_m); };
    /////////////////////////////////
    bool save_theta_swarm = true;

    vector<vector<vector<double> > > finalThetaSwarm;

    int initial_seed = 642345;
    sobj.imfilter(initstates, N, R, J, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, initThetaSwarm, perturb_size, M, coolschedule, thetaTrajectory, finalThetaSwarm, weighted_interaction, max_lookahead, initial_seed, 20); // run IF2

    if (save_theta_swarm) {
      std::ofstream oswarmfile;
      std::string oswarmfilename = save_directory + "thetaswarm_m" + std::to_string(start_m+M) + common_string;      
      oswarmfile.open(oswarmfilename);
      for (int i=0; i<R; i++)
	for (int j=0; j<J; j++) {
	  for (std::size_t theta_dim=0; theta_dim<sobj.theta.size(); theta_dim++)
	    oswarmfile << finalThetaSwarm[i][j][theta_dim] << " ";
	  oswarmfile << endl;
	}
      oswarmfile.close();
    }
    
  }

  if (!run_IF2) {
    vector<double> theta_point;
    theta_point = sobj.theta;
      
    int initial_seed = 642345 + 72521*sobj.nvar + 6442083*exp_no + 12471*int(alpha*815);
    sobj.imfilter(initstates, N, R, J, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, theta_point, weighted_interaction, max_lookahead, initial_seed, 20); // run filter
    M = 1;
  }

  std::clock_t stf_end = clock();

  double stf_time = (double) (stf_end - stf_start) / CLOCKS_PER_SEC;
  cout << "Total time elapsed (s): " << stf_time << endl;

  std::ofstream mfile;
  std::string mfilename = save_directory + "im_est_state_mean" + common_string;
    
  std::ofstream m1file;
  std::string m1filename = save_directory + "im_est_f1_mean" + common_string;
    
  std::ofstream q10file;
  std::string q10filename = save_directory + "im_est_q10" + common_string;
    
  std::ofstream q50file;
  std::string q50filename = save_directory + "im_est_q50" + common_string;

  std::ofstream q90file;
  std::string q90filename = save_directory + "im_est_q90" + common_string;

  std::ofstream tfile;
  std::string tfilename = save_directory + "im_time" + common_string;
    
  std::ofstream lfile;
  std::string lfilename = save_directory + "im_est_likelihood" + common_string;

  std::ofstream ESSfile;
  std::string ESSfilename = save_directory + "im_ESS" + common_string;
    
  std::ofstream thetafile;
  std::string thetafilename = save_directory + "im_theta" + common_string;
    
  mfile.open(mfilename); q10file.open(q10filename); q50file.open(q50filename); q90file.open(q90filename); m1file.open(m1filename); tfile.open(tfilename, std::ofstream::app); lfile.open(lfilename); ESSfile.open(ESSfilename, std::ofstream::app); thetafile.open(thetafilename, std::ofstream::app);

  tfile << stf_time << " ";

  for (int t = 0; t < sobj.ntimes; t++) {
    for (int k = 0; k < sobj.nvar; k++) {
      mfile << state_mean[t][k] << " ";
      m1file << f1_mean[t][k] << " ";
      q10file << q10[t][k] << " ";
      q50file << q50[t][k] << " ";
      q90file << q90[t][k] << " ";
    }
    mfile << endl; m1file << endl; q10file << endl; q50file << endl; q90file << endl;
    
    for (int s = 0; s < N; s++) {
      for (int r = 0; r < R; r++) 
	ESSfile << ESS[r][t][s] << " ";
      ESSfile << endl;
    }
  }

  for (int m=0; m<M; m++)
    for (int t=0; t<sobj.ntimes; t++) {
      lfile << lestimate[m][t] << endl;    
      if (run_IF2) {
	for (std::size_t theta_dim=0; theta_dim<sobj.theta.size(); theta_dim++)
	  thetafile << thetaTrajectory[m][t][theta_dim] << " ";
	thetafile << endl;
      }
    }
  
  mfile.close(); q10file.close(); q50file.close(); q90file.close(); m1file.close(); tfile.close(); lfile.close(); ESSfile.close(); thetafile.close();
  
  return 0;
}
  
