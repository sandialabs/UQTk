#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

typedef std::vector<double> RealVector;

/*************************************************
Define Likelihood function
*************************************************/
class Likelihood{
public:
  Likelihood(){};
  ~Likelihood(){}; 
  double eval(RealVector&);
};

// 3d bimodal likelihood
double Likelihood::eval(RealVector& x){
  double loglik = log(exp(-.5*((x[0]-0.5)*(x[0]-0.5)/.01 + 
    (x[1]-0.5)*(x[1]-0.5)/.01 +
    (x[2]-0.5)*(x[2]-0.5)/.01)) +
    exp(-.5*((x[0]+0.5)*(x[0]+0.5)/.01 + 
      (x[1]+0.5)*(x[1]+0.5)/.01 +
      (x[2]+0.5)*(x[2]+0.5)/.01)));
  return loglik; 
}

/*************************************************
Define Prior function
*************************************************/
class Prior{
public:
  Prior(){};
  ~Prior(){}; 
  double eval(RealVector&);
};


// Gaussian prior
double Prior::eval(RealVector& x){
  double logp = -.5*(x[0]*x[0]/1.0 + x[1]*x[1]/1.0 + x[2]*x[2]/1.0);
  return logp;
}

int main (int argc, const char *argv[]) {

  std::ifstream input_file("mcmcstates_local.dat");
  std::ofstream output_file;
  std::ofstream flag_file;
  std::string line, token;
  std::stringstream iss;
  Likelihood L; 
  Prior P; 

  if (argc > 1) {
    output_file.open("tmcmc_lp.dat");
  } else {
    output_file.open("tmcmc_ll.dat");
  }

  remove( "done.txt" );

  while (std::getline(input_file, line)) {
    RealVector x;
    iss << line;
    double d;
    while (iss >> d) {
      x.push_back(d);
    }

    if (argc > 1) {
      output_file << std::setprecision(18) << P.eval(x) << "\n";
    } else {
      output_file << std::setprecision(18) << L.eval(x) << "\n";
    }
    
    x.clear();
    iss.clear();
  }
  output_file.close();

  flag_file.open("done.txt");
  flag_file << "Done" << std::endl;
  flag_file.close();
  return 0;
}