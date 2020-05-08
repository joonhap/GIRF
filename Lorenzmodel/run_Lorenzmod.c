#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cmath>
#include "../../stf_parallel.h"
#include "../Lorenzmod.h"

using std::vector;
using std::cout; 
using std::endl;

double square(double x) {
  return x*x;
}

double obs_int; // observation time interval

// // // MAIN FUNCTION // // //
int main(int argc, char** args)
{

  // create sobj
  /////////////////////////////////
  int dim = std::stoi(args[1]);
  int ndatapoint = 200;
  stf<double, double> sobj(dim, dim, ndatapoint); // state_dim, measurement_dim, ntimes
  obs_int = 0.1; // make sure the following line is changed accordingly
  std::string obs_int_string = "0.1";

  int theta_size = 3;
  sobj.theta.resize(theta_size);

  const double F = std::stof(args[2]); // external forcing of Lorenz 96 model
  sobj.theta[0] = F;
  const double F_datagen = 8.0; // the value of F used to generate the data
  const double sig_p = 1.0; // state process noise (standard deviation of Brownian motion at time 1)
  sobj.theta[1] = log(sig_p);
  const double sig_m = 1.0; // measurement process noise
  sobj.theta[2] = log(sig_m);

  sobj.pos_ivp.resize(0); // there are no initial value parameters. (in accordance of stf_parallel.tpp, ivp's should come after the non-ivp's in principle.)

  //////////////////////////////////////

  // set initial state
  vector<double> state0(dim, 0.0);
  state0[dim-1] = 0.01; // initial state is origin, except the last coordinate which is at 0.01.

  char write_or_read = 'w'; // data file: read('r'), write('w')

  if (write_or_read == 'w') {
    int simulation_seed = 9152 + dim*7245;
    sobj.simulate(state0, simulation_seed);
    
    // write obs into a file
    std::ofstream ofile;
    std::string ofilename = "../OBSINT_" + obs_int_string + "/obs_dim" + std::to_string(dim) + "F" + std::to_string(F_datagen) + "sig_p" + std::to_string(sig_p) + "sig_m" + std::to_string(sig_m) + ".txt";
    ofile.open(ofilename);
    for (int t = 0; t < sobj.ntimes; t++) {
      for (int k = 0; k < sobj.nvar_obs; k++)
	ofile << sobj.obs[t][k] << " ";
      ofile << endl;
    }
    ofile.close();
    
    // write states into a file
    std::ofstream sfile;
    std::string sfilename = "../OBSINT_" + obs_int_string + "/state_dim" + std::to_string(dim) + "F" + std::to_string(F_datagen) + "sig_p" + std::to_string(sig_p) + "sig_m" + std::to_string(sig_m) + ".txt";
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
    std::string ifilename = "../OBSINT_" + obs_int_string + "/obs_dim" + std::to_string(dim) + "F" + std::to_string(F_datagen) + "sig_p" + std::to_string(sig_p) + "sig_m" + std::to_string(sig_m) + ".txt";
    std::ifstream ifile(ifilename);
    if (!ifile.is_open()) {
      cout << "Data file does not exist." << endl;
      exit(1);
    }

    for (int t = 0; t < sobj.ntimes; t++)
      for (int k = 0; k < sobj.nvar_obs; k++) {
	ifile >> read_double;
	sobj.obs[t][k] = read_double;
      }
  }

  // // // // // // // //
  // space-time filter //
  // // // // // // // //

  int R = std::stoi(args[6]);
  int J = std::stoi(args[7]);
  
  int S = std::stoi(args[9]);

  std::string save_directory = "";

  omp_set_dynamic(0); // Explicitly turn off dynamic threads to make sure global threadprivate variables work properly
    
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
  
  // filenames
  ////////////////////////////////////////
  bool profile = std::stoi(args[4]); // slice or profile estimates? (0 for slice, 1 for profile)
  bool run_IF2 = std::stoi(args[5]); // run IF2? (1 for TRUE, 0 for FALSE)
  bool weighted_interaction = true; // should weighted island interaction be used?
  int max_lookahead = std::stoi(args[8]); // how many observations ahead in time should be considered in the intermediate measuremenet model?
  //  std::function<bool (int, int)> island_interact_condition = [S] (int t, int s) -> bool { return (t%2==1 && s==(S-1)); };
  std::function<bool (int, int)> island_interact_condition = [S] (int t, int s) -> bool { return (s==(S-1)); };

  if (argc==1) {
    cout << "no run id provided." << endl;
    exit(1);
  }
  int exp_no = std::stoi(args[3]);
  std::string data_spec = "_F" + std::to_string(F) + (profile? "_profile_" : "_slice_");
  std::string common_string = "_K" + std::to_string(sobj.nvar_obs) + data_spec + "R" + std::to_string(R) + "J" + std::to_string(J) + "S" + std::to_string(S) + "T" + std::to_string(sobj.ntimes) + "if_" + (run_IF2 ? "T" : "F") + "_wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + "_dtObs" + obs_int_string + ".txt";

  //////////////////////////////////////////
  std::clock_t stf_start = clock();

  // IF2
  if (run_IF2) {
    int start_m = std::stoi(args[10]); // starting IF2 iteration number (if a continued IF2, m is non-zero)
    int M = std::stoi(args[11]); // number of IF2 iterations
    int seed_addition = std::stoi(args[12]); // the number to add to the random seed to start with a new seed in case the IF2 search fails

    vector<vector<vector<double> > > initThetaSwarm(R, vector<vector<double> >(J, sobj.theta));
    vector<vector<vector<double> > > finalThetaSwarm;

    vector<double> perturb_size(theta_size, 0.0);
    perturb_size[0] = 0.00;
    perturb_size[1] = 0.01;
    perturb_size[2] = 0.01;

    for (int m = start_m; m < start_m+M; m+=M) {

      bool continued = (m != 0); // indicator of whether IF2 is continued from a previous one

      // Set up initial parameter swarm
      if (!continued) {
	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	std::mt19937 rng(2382432 + 75231*exp_no + 875335*int(F) + 61539*seed_addition);
	for (int i=0; i<R; i++)
	  for(int j=0; j<J; j++) {
	    double unif = uniform(rng);
	    //initThetaSwarm[i][j][0] = 0.0 + unif*20.0; // F
	    //unif = uniform(rng);
	    initThetaSwarm[i][j][1] = log(.3 + unif*4.7); // sig_p
	    unif = uniform(rng);
	    initThetaSwarm[i][j][2] = log(.3 + unif*4.7); // sig_m
	  }
      }
      else if (continued) {
	std::string read_string;
	std::string iswarmfilename = save_directory + "thetaswarm_m" + std::to_string(m) + common_string;
	std::ifstream iswarmfile(iswarmfilename);
	if (!iswarmfile.is_open()) {
	  cout << "The theta swarm file does not exist." << endl;
	  exit(1);
	}

	for (int i=0; i<R; i++)
	  for (int j=0; j<J; j++)
	    for (std::size_t theta_dim=0; theta_dim < sobj.theta.size(); theta_dim++) {
	      iswarmfile >> read_string;
	      initThetaSwarm[i][j][theta_dim] = std::stof(read_string);
	    }
      }

      std::function<double (int)> coolschedule = [m] (int iter) -> double { return pow(0.92, m + iter); };

      state_mean.clear();
      q10.clear();
      q50.clear();
      q90.clear();
      f1_mean.clear();
      lestimate.clear();
      ESS.clear();
      thetaTrajectory.clear();
      finalThetaSwarm.clear();

      int initial_seed = 590118 + 89267*exp_no + SEED_INC_PER_IF_ITER*m + 352359*int(F) + 87345*seed_addition;

      sobj.imfilter(initstates, S, R, J, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, initThetaSwarm, perturb_size, M, coolschedule, thetaTrajectory, finalThetaSwarm, weighted_interaction, max_lookahead, island_interact_condition, initial_seed);

      // Write parameter estimation results
      bool save_theta_swarm = true;
      if (save_theta_swarm) {
	std::ofstream oswarmfile;
	std::string oswarmfilename = save_directory + "thetaswarm_m" + std::to_string(m+M) + common_string; 
	oswarmfile.open(oswarmfilename);
	for (int i=0; i<R; i++)
	  for (int j=0; j<J; j++) {
	    for (std::size_t theta_dim=0; theta_dim<sobj.theta.size(); theta_dim++)
	      oswarmfile << finalThetaSwarm[i][j][theta_dim] << " ";
	    oswarmfile << endl;
	  }
	oswarmfile.close();
      }

      // Write run results
      std::ofstream lfile;
      std::string lfilename = save_directory + "im_est_likelihood" + common_string;
      std::ofstream thetafile;
      std::string thetafilename = save_directory + "im_theta" + common_string;
      lfile.open(lfilename, std::ofstream::app);
      thetafile.open(thetafilename, std::ofstream::app);

      for (size_t mm=0; mm<lestimate.size(); mm++)
	for (size_t tt=0; tt<lestimate[mm].size(); tt++) 
	  lfile << lestimate[mm][tt] << endl;
      
      for (size_t mm=0; mm<thetaTrajectory.size(); mm++)
	for (int t=0; t<sobj.ntimes; t++) {
	  for (std::size_t theta_dim=0; theta_dim<sobj.theta.size(); theta_dim++)
	    thetafile << thetaTrajectory[mm][t][theta_dim] << " ";
	  thetafile << endl;
	}
      lfile.close(); thetafile.close();
    }
  }
  // non-IF2
  if (!run_IF2) {
    vector<double> theta_point = sobj.theta; // the point of theta at which the likelihood will be evaluated

    if (profile) {
      std::string input_theta_filename = save_directory + "im_theta_K" + std::to_string(sobj.nvar_obs) + "_F" + std::to_string(F) + "_profile_" + "R" + "5" + "J" + "2000" + "S" + std::to_string(S) + "T" + std::to_string(sobj.ntimes) + "if_" + "T" + "_wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + ".txt";
      std::ifstream itf;
      itf.open(input_theta_filename);

      if(itf.is_open()) {
	itf.seekg(-1,std::ios_base::end);                // go to one spot before the EOF
	int newline_char = 0;                     // count the number of newline characters from the end
	while(newline_char < 2) {
	  char ch;
	  itf.get(ch);                            // Get current byte's data
	  if((int)itf.tellg() <= 1) {             // If the data was at or before the 0th byte
	    itf.seekg(0);                       // The first line is the last line
	    break;                // So stop there
	  }
	  else if(ch == '\n') {                   // If the data was a newline
	    itf.seekg(-2,std::ios_base::cur);    // Move to the front of that data, then to the front of the data before it
	    newline_char++; 
	  }
	  else {                                  // If the data was neither a newline nor at the 0 byte
	    itf.seekg(-2,std::ios_base::cur);     // Move to the front of that data, then to the front of the data before it
	  }
	}
	itf.seekg(2,std::ios_base::cur);
	std::string lastLine;            
	getline(itf,lastLine);                      // Read the current line
	std::istringstream iss(lastLine);
	for (size_t theta_dim = 0; theta_dim < sobj.theta.size(); theta_dim++) { // fill in theta vector
	  std::string theta_read_string;
	  iss >> theta_read_string;
	  theta_point[theta_dim] = std::stof(theta_read_string);
	}
	itf.close();
      }
      else {
	cout << "Theta file not found." << endl;
	exit(1);
      }
    }

    int initial_seed = 6548 + 487434*exp_no + 198131*int(F);

    sobj.imfilter(initstates, S, R, J, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, theta_point, weighted_interaction, max_lookahead, island_interact_condition, initial_seed);

    // Write run results
    std::ofstream lfile;
    std::string lfilename = save_directory + "im_est_likelihood" + common_string;
    lfile.open(lfilename, std::ofstream::app);
    for (size_t mm=0; mm<lestimate.size(); mm++)
      for (size_t tt=0; tt<lestimate[mm].size(); tt++) 
	lfile << lestimate[mm][tt] << endl;
    lfile.close(); 
  }

  std::clock_t stf_end = clock();

  double stf_time = (double) (stf_end - stf_start) / CLOCKS_PER_SEC;
  cout << "Total time elapsed (s): " << stf_time << endl;

  std::ofstream mfile;
  std::string mfilename = save_directory + "im_est_state_mean" + common_string;
    
  std::ofstream tfile;
  std::string tfilename = save_directory + "im_time" + common_string;

  std::ofstream ESSfile;
  std::string ESSfilename = save_directory + "im_ESS" + common_string;
    
  mfile.open(mfilename); tfile.open(tfilename, std::ofstream::app); ESSfile.open(ESSfilename); 

  tfile << stf_time << " ";
  for (int t = 0; t < sobj.ntimes; t++) {
    for (int k = 0; k < sobj.nvar; k++) {
      mfile << state_mean[t][k] << " ";
    }
    mfile << endl; 
    for (int s = 0; s < S; s++) {
      for (int r = 0; r < R; r++) 
	ESSfile << ESS[r][t][s] << " ";
      ESSfile << endl;
    }
  }

  mfile.close(); tfile.close(); ESSfile.close(); 

  return 0;
}

