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

//#pragma omp threadprivate(dist, birth_data, birth_data_first_year, birth_data_last_year, pop_data, pop_data_first_year, pop_data_last_year, dist_mean, pop_mean, cityindex)


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

  int theta_size = 17+3*ncity;
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
  const double repOD = 0.25; //reporting overdispersion (log-transformed)
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

  sobj.pos_ivp.resize(3*ncity); // there are 3*ncity initial value parameters.
  for (int k=0; k<3*ncity; k++)
    sobj.pos_ivp[k] = 17+k; 

  //////////////////////////////////////

  // read theta from a estimated theta file
  vector<double> theta_point;
  theta_point.resize(sobj.theta.size());

  std::string data_spec = "_case_G" + std::to_string(300.0);
  bool weighted_interaction = true;
  int max_lookahead = 2;
  int exp_no = 4;
  int N = ncity;
  std::string save_directory = "";

  std::string input_theta_filename = save_directory + "im_theta_K" + std::to_string(sobj.nvar_obs) + data_spec + "R" + "1" + "J" + "4000" + "N" + std::to_string(N) + "T" + std::to_string(sobj.ntimes) + "if_" + "T" + "_wi_" + (weighted_interaction ? "T" : "F") + "_la" + std::to_string(max_lookahead) + "_" + std::to_string(exp_no) + ".txt";
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
    for (size_t theta_dim = 0; theta_dim < sobj.theta.size(); theta_dim++) {
      std::string theta_read_string;
      iss >> theta_read_string;
      theta_point[theta_dim] = std::stof(theta_read_string);
    }
      
    itf.close();
  }
  sobj.theta = theta_point;
  sobj.theta[9] = log(G); // simulate at the given G
  
  // set initial state
  vector<int> state0(5*ncity, 0);

  for (int k=0; k<ncity; k++) {
    state0[5*k+4] = pop_data[startyear-pop_data_first_year][cityindex[k]]; // total population
    state0[5*k+0] = state0[5*k+4] / (1.0+exp(-theta_point[sobj.pos_ivp[3*k+0]]));
    state0[5*k+1] = state0[5*k+4] / (1.0+exp(-theta_point[sobj.pos_ivp[3*k+1]]));
    state0[5*k+2] = state0[5*k+4] / (1.0+exp(-theta_point[sobj.pos_ivp[3*k+2]]));
  }

  // simulate POMP
  int simulation_seed = 9152 + (20-ncity)*7245 + 18719*int(G);
  sobj.simulate(state0, simulation_seed);
    
  // write obs into a file
    
  std::ofstream ofile;
  std::string ofilename = "../obs_ncity" + std::to_string(ncity) + "_G" + std::to_string(int(G)) + "ps(G300_n4).txt";
  ofile.open(ofilename);
  for (int t = 0; t < sobj.ntimes; t++) {
    for (int k = 0; k < sobj.nvar_obs; k++)
      ofile << sobj.obs[t][k] << " ";
    ofile << endl;
  }
  ofile.close();
    
  // write states into a file
  std::ofstream sfile;
  std::string sfilename = "../state_ncity" + std::to_string(ncity) + "_G" + std::to_string(int(G)) + "ps(G300_n4).txt";
  sfile.open(sfilename);
  for (int t = 0; t < sobj.ntimes+1 ; t++) {
    for (int k = 0; k < sobj.nvar; k++)
      sfile << sobj.hidden_state[t][k] << " ";
    sfile << endl;
  }
  sfile.close();
    
  cout << "simulation finished." << endl;


  return 0;
}

