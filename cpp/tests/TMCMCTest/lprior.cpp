#include "lprior.hpp"

double normPDF(double x, double mean, double var) {
  double pval = 1.0 / sqrt(2 * 3.1415926535897 * var *var) *
  exp(- 0.5 * (x - mean) * (x - mean) * (1.0 / (var * var)));

  return pval;
}

double unifPDF(double x, double a, double b) {
  if (a <= x && x <= b) {
    return (1.0 / (b - a));
  } else {
    return pow(10, -300);
  }
}

double priorLL(RealVector &x, CharVector &distr,
  RealVector &paramA, RealVector &paramB) {
  double lp = 0;

  // assert(x.size() == distr.size());

  for (size_t i = 0; i < x.size(); ++i) {
    switch (distr[i]) {
      case 'n':
        lp += log(normPDF(x[i], paramA[i], paramB[i]));
        break;

      case 'u':
        lp += log(unifPDF(x[i], paramA[i], paramB[i]));
        break;
    }
  }

  return lp;
}

void parseLP(std::string setup_fname, std::string sample_fname) {
  std::stringstream iss;
  std::string line;
  std::string token;
  std::ofstream DAT;
  DAT.open("tmcmc_lp.dat");

  CharVector distr;
  RealVector paramA;
  RealVector paramB;

  std::ifstream setup_file(setup_fname);
  std::ifstream sample_file(sample_fname);
  // Check file is open
  if (!setup_file.is_open()) {
    std::cout << "Couldn't open setup file\n";
    throw;
  } else if (!sample_file.is_open()) {
    std::cout << "Couldn't open sample file\n";
    throw;
  }
  // Parse Setup_file
  while (std::getline(setup_file, line)) {
    int step = 0;
    iss << line;
    while (std::getline(iss, token, ' ')) {
      switch (step) {
        case 0:
          if (token[0] != 'n' && token[0] != 'u') {
            std::cout << "Too many parameters for distribution" << std::endl;
            throw("Improperly formatted setup file");
          } else {
            distr.push_back(token[0]);
            step++;
          }
          break;

        case 1:
          paramA.push_back(std::stod(token));
          step++;
          break;

        case 2:
          paramB.push_back(std::stod(token));
          step++;
          if (distr.back() == 'n' || distr.back() == 'u') {
            step = 0;
          }
          break;

        default:
          break;
      }
    }
    line.clear();
    iss.clear();
  }

  // Parse Sample File
  while (std::getline(sample_file, line)) {
    RealVector x(0);
    iss << line;
    while (std::getline(iss, token, ' ')) {
      x.push_back(std::stod(token));
    }
    // assert(x.size() == distr.size());

    DAT << std::setprecision(20) << priorLL(x, distr, paramA, paramB) << '\n';

    x.clear();
    line.clear();
    iss.clear();
  }
}

int main(int argc, const char* argv[]) {
  if (argc != 5) {
    std::cout << "Usage: " << argv[0] << " -f <setup> -s <sample>" <<std::endl;
  } else {
    parseLP(argv[2], argv[4]);
  }
  return 0;
}
