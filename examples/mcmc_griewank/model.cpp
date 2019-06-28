#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
typedef std::vector<double> RealVector;

double rosenbrockFunc(RealVector &x);
double himmelblauFunc(RealVector &x);
double schwefelFunc(RealVector &x);
double griewankFunc(RealVector &x);
double trimodal(RealVector &x);
double multimodalBivariateNorm(RealVector &x);
double HighDMVN(RealVector &x);
double wingWeight(RealVector &x);

int main () {
  std::ifstream input_file("mcmcstates_1.dat");
  std::ofstream output_file("tmcmc_ll.dat");
  std::string line, token;
  std::stringstream iss;

  while (std::getline(input_file, line)) {
    RealVector x;
    iss << line;
    double d;
    while (iss >> d) {
      x.push_back(d);
    }

    /* Change Function HERE */
    output_file << std::setprecision(18) << himmelblauFunc(x) << "\n";
    x.clear();
    iss.clear();
  }
  output_file.close();
  return 0;
}

/* MODEL CODE */

double rosenbrockFunc(RealVector &x) {
  double sum = pow((1 - x[0]),2) + 100 * pow((x[1] - pow(x[0], 2)), 2);

  return -sum;
}
double himmelblauFunc(RealVector &x) {
  double s1 = pow((pow(x[0], 2) + x[1] - 11), 2);
  double s2 = pow((x[0] + pow(x[1], 2) - 7), 2);

  return -(s1 + s2);
}

double schwefelFunc(RealVector &x) {
  double sum = 0.0;

  for (size_t i = 0; i < x.size(); ++i) {
    sum += x[i] * sin(sqrt(abs(x[i])));
  }

  return -(418.9829*x.size() - sum);
}

double griewankFunc(RealVector &x) {
  double prod = 1.0;
  double sum = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    prod *= cos(x[i] / sqrt(i + 1));
  }

  for (size_t i = 0; i < x.size(); ++i) {
    sum += pow(x[i], 2) / 4000;
  }
  return -(sum - prod + 1);
}

double trimodal(RealVector &x) {
  double muA = 0;
  double sigmaA = 1;
  double muB = -50;
  double sigmaB = 5;
  double muC = 50;
  double sigmaC = 3;

  double p1 = 1/(sqrt(2 * 3.1415926535897 * pow(sigmaA,2))) *
                exp(-0.5 * pow((x[0] - muA), 2) / (2 * pow(sigmaA, 2)));

  double p2 = 1/(sqrt(2 * 3.1415926535897 * pow(sigmaB, 2))) *
                exp(-0.5 * pow((x[0] - muB), 2) / (2 * pow(sigmaB, 2)));

  double p3 =  1/(sqrt(2 * 3.1415926535897 * pow(sigmaC, 2))) *
                exp(-0.5 * pow((x[0] - muC), 2) / (2 * pow(sigmaC, 2)));
  return log(0.3 * p1 + 0.5 * p2 + 0.2 * p3);
}

double bivariatePDF(RealVector &mu, RealVector &sigma, RealVector &x,
                      size_t dim1 = 0, size_t dim2 = 1) {
  double p1 = 1/(sqrt(2 * 3.1415926535897 * pow(sigma[dim1], 2))) *
                exp(-0.5 * pow((x[0] - mu[dim1]), 2) / (2 * pow(sigma[dim1], 2)));

  double p2 = 1/(sqrt(2 * 3.1415926535897 * pow(sigma[dim2], 2))) *
                exp(-0.5 * pow((x[0] - mu[dim2]), 2) / (2 * pow(sigma[dim2], 2)));

  return p1 * p2;
}

double multimodalBivariateNorm(RealVector &x) {
  RealVector sigma = {0.5, 0.5};
  RealVector mu1 = {0.5, 0.5};
  RealVector mu2 = {9.5, 0.5};
  RealVector mu3 = {9.5, 9.5};
  RealVector mu4 = {0.5, 9.5};
  RealVector mu5 = {5.0, 5.0};
  RealVector mu6 = {5.0, 9.5};
  RealVector mu7 = {9.5, 5.0};
  RealVector mu8 = {5.0, 0.5};
  RealVector mu9 = {0.5, 5.0};
  RealVector mu10 = {4.5, 4.5};
  RealVector mu11 = {4.5, 5.5};
  RealVector mu12 = {5.5, 4.5};
  RealVector mu13 = {5.5, 5.5};

  RealVector PDFs;
  PDFs.push_back(bivariatePDF(mu1, sigma, x));// 5 5
  PDFs.push_back(bivariatePDF(mu2, sigma, x));// 95 5
  PDFs.push_back(bivariatePDF(mu3, sigma, x));// 95 95
  PDFs.push_back(bivariatePDF(mu4, sigma, x));// 5 95
  PDFs.push_back(bivariatePDF(mu5, sigma, x));// 50 50
  PDFs.push_back(bivariatePDF(mu6, sigma, x));// 50 95
  PDFs.push_back(bivariatePDF(mu7, sigma, x));// 95 50
  PDFs.push_back(bivariatePDF(mu8, sigma, x));// 50 5
  PDFs.push_back(bivariatePDF(mu9, sigma, x));// 5 50
  PDFs.push_back(bivariatePDF(mu10, sigma, x)); // 45 45
  PDFs.push_back(bivariatePDF(mu11, sigma, x)); // 45 55
  PDFs.push_back(bivariatePDF(mu12, sigma, x)); // 55 45
  PDFs.push_back(bivariatePDF(mu13, sigma, x)); // 55 55


  double res = *std::max_element(PDFs.begin(), PDFs.end());
  return log(res);
}

double HighDMVN(RealVector &x) {
  double res = 0;
  double mu = 20;
  double sigma = 4;
  for (size_t i = 0; i < x.size(); ++i) {
    res += log(1.0/sigma) - (pow((x[i] - mu), 2) /(2 * pow(sigma, 2))) - (0.5 * log(2 * 3.1415926535897));
  }

  return res;
}

double wingWeight(RealVector &x) {
  double Sw = x[0]; // 152
  double Wfw = x[1]; // 280
  double A = x[2]; // 6.0
  double LamCaps = x[3]; // -10
  double q = x[4]; // 42
  double lam = x[5]; // 0.5
  double tc = x[6]; // 0.08
  double Nz = x[7]; // 2.6
  double Wdg = x[8]; // 1700
  double Wp = x[9]; // 0.045

  double fact1 = 0.036 * pow(Sw, 0.758) * pow(Wfw, 0.0035);
  double fact2 = pow((A / pow((cos(LamCaps)), 2)), 0.6);
  double fact3 = pow(q, 0.006) * pow(lam, 0.04);
  double fact4 = pow(((100 * tc) / cos(LamCaps)), 2);
  double fact5 = pow((Nz * Wdg), 0.49);

  double term1 = Sw * Wp;

  double res = fact1 * fact2 * fact3 * fact4 * fact5 + term1;

  return -log(res);
}
