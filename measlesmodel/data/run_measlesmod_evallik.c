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

double dist[NCITY][NCITY];
double dist_mean; // average of all pairwise distances. Used for standardizing distances.
double pop_mean; // average of all populations. Used for standardizing populations.
int birth_data[NYEAR][NCITY];
int birth_data_first_year;
int birth_data_last_year;
int pop_data[NYEAR][NCITY];
int pop_data_first_year;
int pop_data_last_year;
int cityindex[NCITY]; // index for which cities will be analyzed (index ranging from 0 to NCITY-1) (in alphabetical order for this dataset)
double school_prop;

//#pragma omp threadprivate(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)


// read the matrix of distances between cities
void read_distance_matrix() {  
  std::ifstream infile("../../DistMatrix.txt");
  if (!infile.is_open()) {
    cout << "Distance data file does not exist." << endl;
    exit(1);
  }
  for(int k=0; k<NCITY; k++) {
    for(int l=0; l<NCITY; l++)
      infile >> dist[k][l];
  }
}

// read the matrix of birth data
void read_birth_data(void) {
  std::ifstream infile("../../UKbirths.csv");
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
    for (int k=0; k<NCITY; k++) {
      std::getline(ss,token,',');
      birth_data[nyear-1][k] = std::stoi(token); // record birth
    }
  }
}

// read population size data
void read_pop_data(void) {
  std::ifstream infile("../../UKpop.csv");
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
    for (int k=0; k<NCITY; k++) {
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
  for (int k=0; k<NCITY; k++)
    for (int l=0; l<k; l++)
      dist_sum += dist[k][l];
  dist_mean = dist_sum / (NCITY*(NCITY-1)/2 + 0.0);
  double pop_sum = 0.0;
  for (int k=0; k<NCITY; k++)
    pop_sum += pop_data[0][k];
  pop_mean = pop_sum / double(NCITY);

  // create sobj
  /////////////////////////////////
  int ncity = 40;
  int ndatapoint = 415;
  stf<int, int> sobj(5*ncity, ncity, ndatapoint); // state_dim, measurement_dim, ntimes

  // define which cities will be analyzed
  for (int k=0; k<ncity; k++)
    cityindex[k] = k;

  int non_ivp_size = 17;
  int theta_size = non_ivp_size+3*ncity;
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
  const double rep = 0.5; // reporting probability (logit-transformed)
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

  sobj.pos_ivp.resize(3*ncity); // there are 3*ncity initial value parameters. ivp's should come after the non-ivp's.
  for (int k=0; k<3*ncity; k++)
    sobj.pos_ivp[k] = non_ivp_size+k; 

  //////////////////////////////////////

  // set initial state
  vector<int> state0(5*ncity, 0);

  for (int k = 0; k < ncity; k++) {
    state0[5*k+4] = pop_data[startyear-pop_data_first_year][cityindex[k]]; // total population
  }
  // read the matrix of measles case data
  vector<vector<int> > case_data(sobj.ntimes, vector<int>(NCITY));

  std::ifstream infile("../../UKmeasles.csv");
  if (!infile.is_open()) {
    cout << "Measles case data file does not exist." << endl;
    exit(1);
  }
  std::string line;
  std::getline(infile, line); // the first row is column names; discard it.

  int biweek_no = 0;
  double tol = 1e-6;
  while(std::getline(infile, line)) {
    std::istringstream ss(line);
    std::string token;
    std::getline(ss, token, ','); // first entry is time
    if (std::stof(token) < startyear+tol) continue; // if the time of the row is not greater than startyear, move to the next line
    
    biweek_no++; // increase week number
    if (biweek_no > sobj.ntimes) break; // if the sobj.obs matrix is full, exit the loop
    for (int k=0; k<NCITY; k++) {
      std::getline(ss,token,',');
      case_data[biweek_no-1][k] = std::stoi(token); // record case data
    }
  }
  
  for (int t = 0; t < sobj.ntimes; t++) 
    for (int k = 0; k < sobj.nvar_obs; k++) 
      sobj.obs[t][k] = case_data[t][cityindex[k]];

  for (int k = 0; k < ncity; k++) {
    state0[5*k+4] = pop_data[startyear-pop_data_first_year][cityindex[k]]; // total population
  }

  // create another sobj with much shorter ntimes, for ivp estimation
  int ntimes_ivp = 3; // the number of observations to be included
  stf<int, int> sobj_ivp(5*ncity, ncity, ntimes_ivp); // state_dim, measurement_dim, ntimes
  sobj_ivp.theta = sobj.theta; // copy theta from sobj
  sobj_ivp.pos_ivp = sobj.pos_ivp; // copy pos_ivp from sobj
  for (int t=0; t<ntimes_ivp; t++)
    sobj_ivp.obs[t] = sobj.obs[t]; // copy obs from sobj


  // // // // // // // //
  // space-time filter //
  // // // // // // // //

  int R = 10;
  int J = 1000;
  
  int S = ncity;

  std::string save_directory = "";

  omp_set_dynamic(0); // Explicitly turn off dynamic threads to make sure global threadprivate variables work properly
    
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
  bool run_IF2 = false;
  bool weighted_interaction = true; // should weighted island interaction be used?
  int max_lookahead = 3; // how many observations ahead in time should be considered in the intermediate measuremenet model?
  //  std::function<bool (int, int)> island_interact_condition = [S] (int t, int s) -> bool { return (t%2==1 && s==(S-1)); };
  std::function<bool (int, int)> island_interact_condition = [S] (int t, int s) -> bool { return (s==(S-1)); };

  if (argc==1) {
    cout << "no run id provided." << endl;
    exit(1);
  }
  int exp_no = std::stoi(args[2]);
  int seed_addition = std::stoi(args[3]); // the number to add to the seed (in case the run without addition fails.)
  std::string data_spec = "_case_G" + std::to_string(G);
  std::string common_string = "_K" + std::to_string(sobj.nvar_obs) + data_spec + "R" + std::to_string(R) + "J" + std::to_string(J) + "S" + std::to_string(S) + "T" + std::to_string(sobj.ntimes) + "if_" + (run_IF2 ? "T" : "F") + "_wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + ".txt";

  //////////////////////////////////////////

  // read the theta value from the last IF2 iteration (non-IVP)
  vector<double> theta_point;
  theta_point.resize(sobj.theta.size());

  std::string input_theta_filename = save_directory + "im_theta_K" + std::to_string(sobj.nvar_obs) + data_spec + "R" + "5" + "J" + "600" + "S" + std::to_string(S) + "T" + std::to_string(sobj.ntimes) + "if_" + "T" + "_wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + ".txt";
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

  bool ivp_exist = std::stoi(args[4]); // ivp_exist=0: do only ivp estimation, ivp_exist=1: do only likelihood evaluation

  if (!ivp_exist) {

    int initial_seed = 8183571 + 766319*exp_no + 98134*int(G) + seed_addition;

    // run IF2 for ivp estimation
    std::function<double (int)> coolschedule_ivp = [] (int iter) -> double { return pow(0.92, 10); };
    vector<double> perturb_size_ivp(theta_size, 0.0); // perturb only the ivp's.
    for (int k=0; k<ncity; k++) {
      perturb_size_ivp[non_ivp_size+3*k] = 0.05; // initial S proportion
      perturb_size_ivp[non_ivp_size+3*k+1] = 0.1; // initial E proportion
      perturb_size_ivp[non_ivp_size+3*k+2] = 0.1; // initial I proportion
    }
    int M_ivp = 80; // mif repetitions for ivp estimation
    int R_ivp = R*10; // for estimation of ivp, run parallel particle filters to keep some amount of variability. Otherwise, all particles will collapse at one point.
    int J_ivp = J/10; 

    // construct new thetaSwarm array with different island blocks
    vector<vector<vector<double> > > initThetaSwarm_ivp(R_ivp, vector<vector<double> >(J_ivp, theta_point));
    vector<vector<vector<int> > > initstates_ivp(R_ivp, vector<vector<int> >(J_ivp, state0));
    vector<vector<vector<double> > > finalThetaSwarm_ivp;

#pragma omp parallel num_threads(R_ivp) copyin(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)
    { } // empty block to create R copies of threadprivate variables.

    sobj_ivp.imfilter(initstates_ivp, S, R_ivp, J_ivp, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, initThetaSwarm_ivp, perturb_size_ivp, M_ivp, coolschedule_ivp, thetaTrajectory, finalThetaSwarm_ivp, weighted_interaction, max_lookahead, island_interact_condition, initial_seed+34681);

    // write ivp estimation results
    std::ofstream ivpfile;
    std::string ivpfilename = save_directory + "im_ivp" + common_string;

    ivpfile.open(ivpfilename, std::ofstream::app);
    for (size_t m_ivp=0; m_ivp < thetaTrajectory.size(); m_ivp++)
      for (int t=0; t<ntimes_ivp; t++) {
	for (size_t ivp_dim = 0; ivp_dim < sobj_ivp.pos_ivp.size(); ivp_dim++) 
	  ivpfile << 1.0/(1.0+exp(-thetaTrajectory[m_ivp][t][sobj_ivp.pos_ivp[ivp_dim]])) << " "; // print ivp estimates
	ivpfile << endl;
      }
    ivpfile.close();

    // write the likelihood estimates for the first few observations used for the ivp estimation
    std::ofstream livpfile;
    std::string livpfilename = save_directory + "im_ivp_est_likelihood" + common_string;

    livpfile.open(livpfilename, std::ofstream::app);
    for (size_t m_ivp=0; m_ivp < lestimate.size(); m_ivp++)
      for (int t=0; t<ntimes_ivp; t++)
	livpfile << lestimate[m_ivp][t] << endl; // print likelihood estimate for debugging
    livpfile.close();
  }
  
  if (ivp_exist) { // likelihood evaluation
#pragma omp parallel num_threads(R) copyin(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)
    { } // empty block to create R copies of threadprivate variables.

    std::string input_ivp_filename = save_directory + "im_ivp" + common_string;
    std::ifstream itf_ivp;
    itf_ivp.open(input_ivp_filename);

    if(itf_ivp.is_open()) {
      itf_ivp.seekg(-1,std::ios_base::end);                // go to one spot before the EOF
      int newline_char = 0;                     // count the number of newline characters from the end
      while(newline_char < 2) {
	char ch;
	itf_ivp.get(ch);                            // Get current byte's data
	if((int)itf_ivp.tellg() <= 1) {             // If the data was at or before the 0th byte
	  itf_ivp.seekg(0);                       // The first line is the last line
	  break;                // So stop there
	}
	else if(ch == '\n') {                   // If the data was a newline
	  itf_ivp.seekg(-2,std::ios_base::cur);    // Move to the front of that data, then to the front of the data before it
	  newline_char++; 
	}
	else {                                  // If the data was neither a newline nor at the 0 byte
	  itf_ivp.seekg(-2,std::ios_base::cur);     // Move to the front of that data, then to the front of the data before it
	}
      }
      itf_ivp.seekg(2,std::ios_base::cur);
      std::string lastLine;            
      getline(itf_ivp,lastLine);                      // Read the current line
      std::istringstream iss(lastLine);
      for (size_t ivp_dim = 0; ivp_dim < sobj.pos_ivp.size(); ivp_dim++) { // fill in the ivp components
	std::string ivp_read_string;
	iss >> ivp_read_string;
	theta_point[non_ivp_size + ivp_dim] = -log(1/std::stof(ivp_read_string) -1);
      }
      itf_ivp.close();
    }
    
    int initial_seed = 8183571 + 766319*exp_no + 98134*int(G) + seed_addition;

    vector<vector<vector<double> > > initThetaSwarm(R, vector<vector<double> >(J, theta_point));

    vector<double> perturb_size(theta_size, 0.0); // no perturbation
    std::function<double (int)> coolschedule = [] (int iter) -> double { return 1.0; };
    vector<vector<vector<double> > > finalThetaSwarm;
    
    sobj.imfilter(initstates, S, R, J, lestimate, state_mean, q10, q50, q90, square, f1_mean, ESS, initThetaSwarm, perturb_size, 1, coolschedule, thetaTrajectory, finalThetaSwarm, weighted_interaction, max_lookahead, island_interact_condition, initial_seed); // run filter

    // write results
    std::ofstream lfile;
    std::string lfilename = save_directory + "im_est_likelihood" + common_string;
    
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
    
    std::ofstream ESSfile;
    std::string ESSfilename = save_directory + "im_ESS" + common_string;
    
    lfile.open(lfilename); mfile.open(mfilename); q10file.open(q10filename); q50file.open(q50filename); q90file.open(q90filename); m1file.open(m1filename); ESSfile.open(ESSfilename); 

    for (int m=0; m<1; m++)
      for (size_t t=0; t<lestimate[m].size(); t++) 
	lfile << lestimate[m][t] << endl;

    for (int t = 0; t < sobj.ntimes; t++) {
      for (int k = 0; k < sobj.nvar; k++) {
	mfile << state_mean[t][k] << " ";
	m1file << f1_mean[t][k] << " ";
	q10file << q10[t][k] << " ";
	q50file << q50[t][k] << " ";
	q90file << q90[t][k] << " ";
      }
      mfile << endl; m1file << endl; q10file << endl; q50file << endl; q90file << endl;
      for (int s = 0; s < S; s++) {
	for (int r = 0; r < R; r++) 
	  ESSfile << ESS[r][t][s] << " ";
	ESSfile << endl;
      }
    }

    lfile.close(); mfile.close(); q10file.close(); q50file.close(); q90file.close(); m1file.close(); ESSfile.close(); 
  }
  
  return 0;
}
