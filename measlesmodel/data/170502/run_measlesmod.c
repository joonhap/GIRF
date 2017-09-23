#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cmath>
#include "../../../stf_parallel.h"
#include "../../measlesmod.h"

using std::vector;
using std::cout; 
using std::endl;

double square(double x) {
  return x*x;
}

double dist[20][20];
double dist_mean; // average of all pairwise distances. Used for standardizing distances.
double pop_mean; // average of all populations. Used for standardizing populations.
int birth_data[21][20];
int birth_data_first_year;
int birth_data_last_year;
int pop_data[21][20];
int pop_data_first_year;
int pop_data_last_year;
int cityindex[20]; // index for which cities will be analyzed (index starting from 0; London=19)
double school_prop;

#pragma omp threadprivate(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)


// read the matrix of distances between cities
void read_distance_matrix() {  
  std::ifstream infile("../../DistMatrix.txt");
  if (!infile.is_open()) {
    cout << "Distance data file does not exist." << endl;
    exit(1);
  }
  for(int k=0; k<20; k++) {
    for(int l=0; l<20; l++)
      infile >> dist[k][l];
  }
}

// read the matrix of birth data
void read_birth_data(void) {
  std::ifstream infile("../../20births.csv");
  if (!infile.is_open()) {
    cout << "Birth data file does not exist." << endl;
    exit(1);
  }
  std::string line;
  std::getline(infile, line); // the first row contains city names; discard it.
  bool is_first_year = true;
  int nyear = 0;
  while(std::getline(infile, line)) {
    std::istringstream ss(line);
    std::string token;
    std::getline(ss, token, ','); // first entry is year
    if (is_first_year) {
      birth_data_first_year = std::stoi(token);
      is_first_year = false;
    }
    birth_data_last_year = std::stoi(token);
    nyear++;
    for (int k=0; k<20; k++) {
      std::getline(ss,token,',');
      birth_data[nyear-1][k] = std::stoi(token); // record birth
    }
  }
}

// read population size data
void read_pop_data(void) {
  std::ifstream infile("../../20pop.csv");
  if (!infile.is_open()) {
    cout << "Population data file does not exist." << endl;
    exit(1);
  }
  std::string line;
  std::getline(infile, line); // the first row contains city names; discard it.
  bool is_first_year = true;
  int nyear = 0;
  while(std::getline(infile, line)) {
    std::istringstream ss(line);
    std::string token;
    std::getline(ss, token, ','); // first entry is year
    if (is_first_year) {
      pop_data_first_year = std::stoi(token);
      is_first_year = false;
    }
    pop_data_last_year = std::stoi(token);
    nyear++;
    for (int k=0; k<20; k++) {
      std::getline(ss,token,',');
      pop_data[nyear-1][k] = std::stoi(token); // record population
    }
  }
}



// // // MAIN FUNCTION // // //
int main(int argc, char** args)
{

  read_distance_matrix();  // read the matrix of distances between cities
  read_birth_data();  // read the matrix of birth data
  read_pop_data(); // read population data matrix

  // compute average distance and population size
  double dist_sum = 0.0; 
  for (int k=0; k<20; k++)
    for (int l=0; l<k; l++)
      dist_sum += dist[k][l];
  dist_mean = dist_sum / (20*19/2 + 0.0);
  double pop_sum = 0.0;
  for (int k=0; k<20; k++)
    pop_sum += pop_data[0][k];
  pop_mean = pop_sum / 20.0;

  // create sobj
  /////////////////////////////////
  int ncity = 20;
  stf<int, int> sobj(5*ncity, ncity, 832); // state_dim, measurement_dim, ntimes

  // define which cities will be analyzed
  for (int k=0; k<ncity; k++)
    cityindex[k] = k;
  //ncity=2
  //  cityindex[0] = 8; // Hastings
  //  cityindex[1] = 19; // London
  //ncity=5
  //cityindex[0] = 0; // Halesworth
  //cityindex[1] = 4; // Oswestry
  //cityindex[2] = 8; // Hastings
  //cityindex[3] = 14; // Sheffield
  //cityindex[4] = 19; // London

  int theta_size = 18;
  sobj.theta.resize(theta_size);

  school_prop = .739; // proportion of the year taken up by school term (270/365.25)

  const double R0 = 25.0; // basic reproduction number (matches beta_s = 33.0)
  sobj.theta[0] = log(R0);
  const double amp = .163; // seasonal amplitude (matches holiday effect=.7)  (for formula, see He et al (2009))  0 < amp < .68
  sobj.theta[1] = log(amp/(.68-amp));
  const double alpha = 0.97; // mixing exponent
  sobj.theta[2] = alpha;
  const double mu = 0.0166/52; // mortality rate
  sobj.theta[3] = mu;
  const double gen = 2.0; // length of one generation cycle of measles transmission (defined as 1/nu_EI + 1/nu_IR)
  sobj.theta[4] = log(gen);
  const double infec_period_prop = .5; // infectious period proportion compared to the generation cycle (=1/nu_EI + 1/nu_IR) (logit-transformed)
  sobj.theta[5] = log(infec_period_prop/(1.0-infec_period_prop));
  const double sigma2 = 0.08; // variability constant (for gamma process, sigma^2. see Box 1 in He et al 2009) (log-transformed)
  sobj.theta[6] = log(sigma2);
  const double rep = std::stof(args[2]); // reporting probability (logit-transformed)
  sobj.theta[7] = log(rep/(1.0-rep));
  const double repOD = 0.25; // reporting overdispersion (log-transformed)
  sobj.theta[8] = log(repOD);
  const double G = std::stof(args[1]); // gravity model coefficient (log-transformed)
  sobj.theta[9] = log(G);
  const double source_pow = 1.0; // exponent to source population size in gravity model
  sobj.theta[10] = source_pow;
  const double dest_pow = 1.0; // exponent to destination population size in gravity model
  sobj.theta[11] = dest_pow;
  const double dist_pow = 1.0; // exponent to distance in gravity model
  sobj.theta[12] = dist_pow;
  const double cohort_entry = 0.4; // cohort entry fraction
  sobj.theta[13] = log(cohort_entry/(1.0-cohort_entry));
  const int schoolstart = 251; // school start day 
  sobj.theta[14] = double(schoolstart);
  const int entryage = 4; // school entry age 
  sobj.theta[15] = double(entryage);
  const int startyear = 1949; // start year (the year at which t=0)
  sobj.theta[16] = double(startyear);
  double Sprop = 0.04; // initial proportion of S population
  sobj.theta[17] = log(Sprop/(1.0-Sprop));

  sobj.pos_ivp.resize(1); // there is one initial value parameter.
  sobj.pos_ivp[0] = 17; 

  //////////////////////////////////////

  // set initial state
  vector<int> state0(5*ncity, 0);

  char write_or_read = 'r'; // data file: read('r'), write('w'), or read real case data(20measles.csv) ('c')

  for (int k = 0; k < ncity; k++) {
    state0[5*k+4] = pop_data[startyear-pop_data_first_year][cityindex[k]]; // total population
    state0[5*k] = int(Sprop * state0[5*k+4]); // initial susceptible population
  }
  if (write_or_read == 'w' || write_or_read == 'r') {
    for (int k=0; k<ncity; k++) {
      state0[5*k+1] = int(0.00027*state0[5*k+4]);
      state0[5*k+2] = int(0.00032*state0[5*k+4]);
    }
  }

  // read data file
  if (write_or_read == 'w') {
    int simulation_seed = 9152 + (20-ncity)*7245;
    sobj.simulate(state0, simulation_seed);
    
    // write obs into a file
    
    std::ofstream ofile;
    std::string ofilename = "../obs_ncity" + std::to_string(ncity) + "_ps9.txt";
    ofile.open(ofilename);
    for (int t = 0; t < sobj.ntimes; t++) {
      for (int k = 0; k < sobj.nvar_obs; k++)
	ofile << sobj.obs[t][k] << " ";
      ofile << endl;
    }
    ofile.close();
    
    // write states into a file
    std::ofstream sfile;
    std::string sfilename = "../state_ncity" + std::to_string(ncity) + "_ps9.txt";
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
    int read_int = 0;
    std::string ifilename = "../obs_ncity" + std::to_string(ncity) + "_ps9.txt";
    std::ifstream ifile(ifilename);
    if (!ifile.is_open()) {
      cout << "Data file does not exist." << endl;
      exit(1);
    }

    for (int t = 0; t < sobj.ntimes; t++)
      for (int k = 0; k < sobj.nvar_obs; k++) {
	ifile >> read_int;
	sobj.obs[t][k] = read_int;
      }
  }

  // // // // // // // //
  // space-time filter //
  // // // // // // // //

  int R = 1;
  int J = 4000;

  int N = ncity;

  std::string save_directory = "";

  omp_set_dynamic(0); // Explicitly turn off dynamic threads to make sure global threadprivate variables work properly

#pragma omp parallel num_threads(R) copyin(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)
  { } // empty block to create R copies of threadprivate variables.
  
  // storage matrix
  vector<vector<double> > state_mean;
  vector<vector<int> > q10;
  vector<vector<int> > q50;
  vector<vector<int> > q90;
  vector<vector<double> > f1_mean;
  vector<vector<double> > lestimate;
  vector<vector<vector<double> > > ESS;
  vector<vector<vector<double> > > thetaTrajectory;
  
  vector<vector<vector<int> > > initstates(R, vector<vector<int> >(J, state0));
  
  // filenames
  ////////////////////////////////////////
  bool run_IF2 = true;
  bool weighted_interaction = true; // should weighted island interaction be used?
  int max_lookahead = 2; // how many observations ahead in time should be considered in the intermediate measuremenet model?
  //  std::function<bool (int, int)> island_interact_condition = [N] (int t, int s) -> bool { return (t%2==1 && s==(N-1)); };
  std::function<bool (int, int)> island_interact_condition = [N] (int t, int s) -> bool { return (s==(N-1)); };

  if (argc==1) {
    cout << "no run id provided." << endl;
    exit(1);
  }
  int exp_no = std::stoi(args[3]);
  int seed_addition = std::stoi(args[4]); // the number to add to the seed (in case the run without addition fails.)
  std::string data_spec = "G" + std::to_string(G) + "rep" + std::to_string(rep);
  std::string common_string = "_K" + std::to_string(sobj.nvar_obs) + data_spec + "R" + std::to_string(R) + "J" + std::to_string(J) + "N" + std::to_string(N) + "T" + std::to_string(sobj.ntimes) + "if_" + (run_IF2 ? "T" : "F") + "_wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + ".txt";

  //////////////////////////////////////////
  int M = std::stoi(args[6]);
    
  std::clock_t stf_start = clock();

  if (run_IF2) {
    vector<vector<vector<double> > > initThetaSwarm(R, vector<vector<double> >(J, sobj.theta));
    vector<vector<vector<double> > > finalThetaSwarm;

    int start_m = std::stoi(args[5]); // starting IF2 iteration number (if a continued IF2, m is non-zero)
    
    vector<double> perturb_size(theta_size, 0.0);
    perturb_size[0] = 0.01;
    perturb_size[1] = 0.01;
    perturb_size[2] = 0.001;
    perturb_size[4] = 0.01;
    perturb_size[5] = 0.01;
    perturb_size[6] = 0.01;
    //perturb_size[7] = 0.01;
    perturb_size[8] = 0.01;
    perturb_size[13] = 0.01;
    
    for (int m = start_m; m < start_m+M; m+=M) {

      bool continued = (m != 0); // indicator of whether IF2 is continued from a previous one

      if (!continued) {
	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	std::mt19937 rng(2382432 + 75231*exp_no + 875335*int(G) + 8157 * seed_addition); 
	for (int i=0; i<R; i++)
	  for(int j=0; j<J; j++) {
	    initThetaSwarm[i][j][0] = log(5.0) + (log(50.0)-log(5.0))*uniform(rng); // R0
	    initThetaSwarm[i][j][1] = log(0.05/(.68-0.05)) + (log(.5/(.68-.5))-log(.05/(.68-.05)))*uniform(rng); // amp
	    initThetaSwarm[i][j][2] = 0.92 + (1.05-0.92)*uniform(rng); // alpha
	    initThetaSwarm[i][j][4] = log(0.7) + (log(3.5)-log(0.7))*uniform(rng); // gen
	    initThetaSwarm[i][j][5] = log(.2/.8) + (log(.8/.2)-log(.2/.8))*uniform(rng); // infec_period_prop
	    initThetaSwarm[i][j][6] = log(0.02) + (log(0.20)-log(0.02))*uniform(rng); // sigma2
	    //initThetaSwarm[i][j][7] = log(.2/.8) + (log(.8/.2)-log(.2/.8))*uniform(rng); // rep
	    initThetaSwarm[i][j][8] = log(.1) + (log(.5)-log(.1))*uniform(rng); // repOD
	    initThetaSwarm[i][j][13] = log(.1/.9) + (log(.8/.2)-log(.1/.9))*uniform(rng); // cohort_entry
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

      int initial_seed = 590118 + 89267*exp_no + SEED_INC_PER_IF_ITER*m + 352359*int(G) + 7127*seed_addition;
      sobj.imfilter(initstates, N, R, J, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, initThetaSwarm, perturb_size, M, coolschedule, thetaTrajectory, finalThetaSwarm, weighted_interaction, max_lookahead, island_interact_condition, initial_seed); // run IF2

      bool save_theta_swarm = true;
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

  std::clock_t stf_end = clock();

  double stf_time = (double) (stf_end - stf_start) / CLOCKS_PER_SEC;
  cout << "Total time elapsed (s): " << stf_time << endl;

  std::ofstream mfile;
  std::string mfilename = save_directory + "im_est_state_mean" + common_string;
    
  std::ofstream q10file;
  std::string q10filename = save_directory + "im_est_q10" + common_string;
    
  std::ofstream q50file;
  std::string q50filename = save_directory + "im_est_q50" + common_string;

  std::ofstream q90file;
  std::string q90filename = save_directory + "im_est_q90" + common_string;

  std::ofstream m1file;
  std::string m1filename = save_directory + "im_est_f1_mean" + common_string;
    
  std::ofstream tfile;
  std::string tfilename = save_directory + "im_time" + common_string;
    
  std::ofstream lfile;
  std::string lfilename = save_directory + "im_est_likelihood" + common_string;
    
  std::ofstream ESSfile;
  std::string ESSfilename = save_directory + "im_ESS" + common_string;
    
  std::ofstream thetafile;
  std::string thetafilename = save_directory + "im_theta" + common_string;
    
  mfile.open(mfilename); q10file.open(q10filename); q50file.open(q50filename); q90file.open(q90filename); m1file.open(m1filename); tfile.open(tfilename, std::ofstream::app); lfile.open(lfilename, std::ofstream::app); ESSfile.open(ESSfilename); thetafile.open(thetafilename, std::ofstream::app);

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

  mfile.close(); q10file.close(); q50file.close(); q90file.close(); m1file.close(); tfile.close(); lfile.close(); ESSfile.close(); thetafile.close();

  return 0;
}

