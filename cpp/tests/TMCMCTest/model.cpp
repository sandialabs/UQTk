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

double normPDF(double x, double mean, double delta);
double normLogPDF(double x, double mean, double delta);
double CauchyDistr(double x, double x0, double gamma);
double bivariatePDF(RealVector &x, double m1, double s1, double m2, double s2);

double lineNormal(RealVector &x, double delta);
double lineUniform(RealVector &x, double delta);
double hoopUniform(RealVector &x, double delta);
double hoopNormal(RealVector &x, double delta);
double sphereUniform(RealVector &x, double delta);
double sphereNormal(RealVector &x, double delta);
double nonlinearUniform(RealVector &x, double delta);
double nonlinearNormal(RealVector &x, double delta);
double sineUniform(RealVector &x, double delta);
double discontNorm(RealVector &x, double delta);
double fourSquareUnif(RealVector &x, double delta);
double nSphere(RealVector &x, double delta);
double nDataset(RealVector &x, double delta);
double paraboloidConstraints(RealVector &x, double delta);
double torus(RealVector &x, double delta);
double sumConstraints(RealVector &x, double delta);
double nonLinearConstraints(RealVector &x, double delta);
double UQTk(RealVector &x, double delta);
double OneDTrimodalNorm(RealVector &x, double delta);
double MixNormUnif(RealVector &x);
double TwoDMultimodal(RealVector &x);
double schwefelFunc(RealVector &x);
double MultiDMultimodal(RealVector &x);
double MultiDUnimodal(RealVector &x);
double MultiDUnimodalGradLog(double x);
double MultiDBimodal(RealVector &x);
double MultiDBimodalGradLog(double x, double ndim);

int main (int argc, const char *argv[]) {
  std::ifstream input_file("mcmcstates_local.dat");
  std::ofstream output_file;
  std::ofstream gradLog_file;
  std::ofstream flag_file;
  std::string line, token;
  std::stringstream iss;
  std::ifstream delta_file("delta.dat");
  double delta;

  delta_file >> delta;

  if (argc > 2) {
    std::string deltaStr(argv[1]);
    delta = std::stod(deltaStr);
    output_file.open("tmcmc_lp.dat");
  } else {
    output_file.open("tmcmc_ll.dat");
  }

  remove( "done.txt" );

  gradLog_file.open("gradlog.dat");
  while (std::getline(input_file, line)) {
    RealVector x;
    iss << line;
    double d;
    while (iss >> d) {
      x.push_back(d);
    }

    double ndim = x.size();
    for (size_t i = 0; i < ndim; ++i) {
      gradLog_file << std::setprecision(18) << MultiDBimodalGradLog(x[i],ndim) << " ";
    }

    gradLog_file << "\n";
    output_file << std::setprecision(18) << MultiDBimodal(x) << "\n";
    x.clear();
    iss.clear();
  }
  output_file.close();

  flag_file.open("done.txt");
  flag_file << "Done" << std::endl;
  flag_file.close();
  return 0;
}

/* MODEL CODE */

double normPDF(double x, double mean = 0, double delta = 1) {
  double p1 = 1.0 / (sqrt(2 * 3.14159265358979323 * pow(delta, 2))) *
              exp(-0.5 * pow((x - mean), 2) / pow(delta, 2));
  return p1;
}

double normLogPDF(double x, double mean = 0, double delta = 1) {
  return log(1.0/delta) - (pow((x - mean), 2) / (2 * pow(delta, 2)))
        - (0.5 * log(2 * 3.14159265358979323));
}

double CauchyDistr(double x, double x0, double gamma) {
  double p1 = (1.0 / (3.14159265358979323 * gamma)) *
              (pow(gamma, 2) / (pow((x - x0), 2) + pow(gamma, 2)));

  return p1;
}

double lineNormal(RealVector &x, double delta = 1) {
  /* Manifold is a straight y = 2 line from x = 0 to x = 10 */
  double mean = 2;

  double pval;
  if (x[0] < 0 || x[0] > 10) {
    pval = pow(10, -50);
  } else {
    pval = normLogPDF(x[1], mean, delta);
  }

  if (x[0] > 10) {
    // std::cout << "PVAL: " << pval << std::endl;
  }
  return pval;
}

double lineUniform(RealVector &x, double delta = 1) {
  /* Same manifold as above, squeezed by uniform */
  double mean = 5;

  double pval;
  if (x[0] < 0 || x[0] > 10 || x[1] < (mean - delta) || x[1] > mean + delta) {
    pval = pow(10, -50);
  } else {
    pval = log(0.1 / (2*delta));
  }
  return pval;
}

double hoopUniform(RealVector &x, double delta = 1) {
  /* A hoop with thickness 2 delta, radius of 10, center at (5, 10) */
  double radius = 10;
  double centerx = 5;
  double centery = 10;

  double r = sqrt(pow((x[0] - centerx), 2) + pow((x[1] - centery), 2));

  double pval;
  if (r > radius + delta || r < radius - delta) {
    pval = 0.00000000000001;
  } else {
    pval = 1 / (4* radius * delta * 3.1415926535897);
  }

  return log(pval);
}

double hoopNormal(RealVector &x, double delta = 1) {
  /* A hoop with thickness 2*delta, radius of 10, center at (5, 10) */
  double radius = 10;
  double centerx = 5;
  double centery = 10;

  double r = sqrt(pow((x[0] - centerx), 2) + pow((x[1] - centery), 2));

  double pval;

  pval = normLogPDF(r, radius, delta);

  return pval;
}

double sphereUniform(RealVector &x, double delta = 1) {
  /* A sphere with thickness 2*delta, radius of 10, center at (5, 5, 5) */
  double radius = 10;
  double centerx = 5;
  double centery = 5;
  double centerz = 5;

  double r = sqrt(pow((x[0] - centerx), 2) + pow((x[1] - centery), 2) + pow((x[2] - centerz), 2));

  double pval;

  if (r > radius + delta || r < radius + delta) {
    pval = 0.00000000000001;
  } else {
    pval = 1.0 /(6 * pow(r, 2) * delta + 2 * pow(delta, 2));
  }

  return log(pval);
}

double sphereNormal(RealVector &x, double delta = 1) {
  /* Sphere with thickness 2*delta, radius of 10, center at (5, -5, 5) */
  double radius = 10;
  double centerx = 5;
  double centery = -5;
  double centerz = 5;

  double r = sqrt(pow((x[0] - centerx), 2) + pow((x[1] - centery), 2) + pow((x[2] - centerz), 2));

  double pval;
  pval = normLogPDF(r, radius, delta);

  return pval;
}

double nonlinearUniform(RealVector &x, double delta = 1) {
  /* Nonlinear function following xy = 4 */
  double constant = 4;

  double pval;
  if (x[0]*x[1] > constant + delta || x[0]*x[1] < constant - delta) {
    pval = -100;
  } else {
    pval = log(0.25);
  }

  return pval;
}

double nonlinearNormal(RealVector &x, double delta = 1) {
  /* Nonlinear function following xy = 4 */
  double constant = 4;

  double deviate = constant - x[0]*x[1];
  double pval;

  pval = log(1.0/delta) - (pow((deviate), 2) /(2 * pow(delta, 2))) - (0.5 * log(2 * 3.1415926535897));

  return pval;
}


double sineUniform(RealVector &x, double delta = 1) {
  /* Sine wave */

  double deviate = x[1] - sin(x[0]);

  double pval;
  if (deviate < -delta || deviate > delta) {
    pval = 0.00000000000001;
  } else {
    pval = 0.25;
  }

  return log(pval);
}

double discontNorm(RealVector &x, double delta = 1) {
  /* Discontinuous linear function */

  double A1 = -10;
  double A2 = -5;
  double Amean = 0;
  double B1 = 0;
  double B2 = 5;
  double Bmean = 5;

  double pval;
  if (x[0] < A1 || (x[0] > A2 && x[0] < B1) || x[0] > B2) {
    pval = -20;
  } else {
    if (x[0] >= A1 && x[0] <= A2) {
      pval = normLogPDF(x[1], Amean, delta);
    } else {
      pval = normLogPDF(x[1], Bmean, delta);
    }
  }

  return pval;
}

double fourSquareUnif(RealVector &x, double delta = 1) {
  /* Four squares */
  double X1 = 0;
  double X2 = 5;
  double Y1 = 0;
  double Y2 = 5;

  delta = 0.5;
  double pval;
  if (x[0] < X1 - delta || (x[0] > X1 + delta && x[0] < X2 - delta) || x[0] > X2 + delta) {
    pval = -20;
  } else if (x[1] < Y1 - delta || (x[1] > Y1 + delta && x[1] < Y2 - delta) || x[1] > Y2 + delta) {
    pval = -20;
  } else {
    pval = -0.01;
  }

  return pval;
}

double nSphere(RealVector &x, double delta = 1) {
  /* N-dim sphere */
  double radius = 5;

  // L2 Norm
  double norm = sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));
  double deviate = norm - radius;

  double pval = normLogPDF(deviate, 0, delta);

  return pval;
}

double nDataset(RealVector &x, double delta = 1) {
  /* Sample from "datasets" with specific 1st, 2nd moments */
  double mean = 10;
  double variance = 10;
  double max = 50;
  double min = -50;

  double sampleMean = std::accumulate(x.begin(), x.end(), 0.0) / (double) x.size();

  std::transform(x.begin(), x.end(), x.begin(),
              [sampleMean](double d) {return pow((d - sampleMean), 2);});

  double sampleVar = std::accumulate(x.begin(), x.end(), 0.0) / (double) (x.size() - 1);

  double p1 = normLogPDF(sampleMean, mean, delta);

  double p2 = normLogPDF(sampleVar, variance, delta);

  for (size_t i = 0; i < x.size(); ++i) {
    if (x[i] > max || x[i] < min) {
      return pow(10, -20);
    } else {
      return p1 + p2;
    }
  }
}

double paraboloidConstraints(RealVector &x, double delta = 1) {
  /* 3 or 6 dim paraboloid constraint */
  double a1 = 1.0;
  double b1 = 3.0;
  double a2 = 2.0;
  double b2 = 2.0;

  double deviate1 = pow(x[0],2) / pow(a1, 2) + pow(x[1], 2) / pow(b1, 2) - x[2];

  double p1 = normLogPDF(deviate1, 0, delta);

  if (x.size() == 6) {
    double deviate2 = pow(x[3], 2) / pow(a2, 2) + pow(x[4], 2) / pow(b2, 2) - x[5];

    double p2 = normLogPDF(deviate2, 0, delta);

    return p1 + p2;
  }

  return p1;
}

double torus(RealVector &x, double delta = 1) {
  /* Torus, a is inner circle radius, c is torus tube radius */
  double a = 2;
  double c = 5;

  double deviate = pow((c - sqrt(pow(x[0], 2) + pow(x[1], 2))), 2)
                   + pow(x[2], 2) - pow(a, 2);

  double p1 = normLogPDF(deviate, 0, delta);

  return p1;
}

double sumConstraints(RealVector &x, double delta = 1) {
  /* Sum of pairs of dimensions have to add up */
  RealVector goal = {10.0, -15.0, 5.0, -8.0, 15.0};
  double res = 0;

  size_t ndims = x.size();
  size_t dimCut = ndims/5;

  for (size_t i = 0; i < 5; ++i) {
    double sum = 0;
    for (size_t j = dimCut * i; j < dimCut * (i + 1); ++j) {
      sum = sum + x[j];
    }
    double deviate = sum - goal[i];
    res = res + normLogPDF(deviate, 0, delta);
  }

  return res;
}

double nonLinearConstraints(RealVector &x, double delta = 1) {
  /* Pairs of dimensions in ellipses */
  RealVector centersX = {0.0, 5.0, -5.0, 2.5, -2.5};
  RealVector centersY = {0.0, 0.0, 0.0, -2.5, -2.5};
  RealVector semiMajor = {2.0, 2.0, 2.0, 2.0, 2.0};
  RealVector semiMinor = {2.0, 2.0, 2.0, 2.0, 2.0};
  double res = 0;

  for (size_t i = 0; i < centersX.size(); ++i) {
    double deviate = pow((x[2*i] - centersX[i]), 2) / pow(semiMajor[i], 2)
                     + pow((x[2*i + 1] - centersY[i]), 2) / pow(semiMinor[i], 2) - 1;
    res += normLogPDF(deviate, 0, delta);

  }
  return res;
}


double UQTk(RealVector &x, double delta = 1) {
  /* Draw UQTk */
  double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;

  // U
  if (x[0] < -7 || x[0] > -4) {
    p1 = pow(10, -50);
  } else {
    p1 = normLogPDF(x[1], 0, delta);
  }
  if (x[3] > 5 || x[3] < 0) {
    p2 = pow(10, -50);
  } else {
    p2 = normLogPDF(x[2], -7, delta);
  }
  if (x[5] > 5 || x[5] < 0) {
    p3 = pow(10, -50);
  } else {
    p3 = normLogPDF(x[4], -4, delta);
  }

  // Q
  double p4dev = pow((x[6] + 1.5), 2) / pow(1.5, 2) + pow((x[7] - 2.5), 2) / pow(2.5, 2) - 1;
  p4 = normLogPDF(p4dev, 0, delta);

  if (x[8] < -2 || x[8] > 0) {
    p5 = pow(10, -50);
  } else {
    double p5dev = x[9] + (4.0/3.0) * x[8];
    p5 = normLogPDF(p5dev, 0, delta);
  }

  // T
  if (x[10] < 1 || x[10] > 4) {
    p6 = pow(10, -50);
  } else {
    p6 = normLogPDF(x[11], 5, delta);
  }
  if (x[13] > 5 || x[13] < 0) {
    p7 = pow(10, -50);
  } else {
    p7 = normLogPDF(x[12], 2.5, delta);
  }

  // k
  if (x[15] > 5 || x[15] < 0) {
    p8 = pow(10, -50);
  } else {
    p8 = normLogPDF(x[14], 5, delta);
  }
  if (x[16] < 5 || x[16] > 7.5) {
    p9 = pow(10, -50);
  } else {
    double p9dev = x[17] + (2.0 / 2.5) * x[16] - 6;
    p9 = normLogPDF(p9dev, 0, delta);
  }

  if (x[18] < 5 || x[18] > 7.5) {
    p10 = pow(10, -50);
  } else {
    double p10dev = x[19] - 0.4 * x[18];
    p10 = normLogPDF(p10dev, 0, delta);
  }

  return p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10;
}


double OneDTrimodalNorm(RealVector &x, double delta = 1) {
  RealVector means = {-20.0, 0.0, 15.0};
  delta = 1;
  double p1 = 1.0 / sqrt(2 * 3.14159265358979323 * pow(2, 2)) *
              exp(-0.5 * (1.0 / pow(2, 2)) * pow((x[0] - means[0]), 2));
  double p2 = 1.0 / sqrt(2 * 3.14159265358979323 * pow(1, 2)) *
              exp(-0.5 * (1.0 / pow(1, 2)) * pow((x[0] - means[1]), 2));
  double p3 = 1.0 / sqrt(2 * 3.14159265358979323 * pow(4, 2)) *
              exp(-0.5 * (1.0 / pow(4, 2)) * pow((x[0] - means[2]), 2));

  return log(0.5 * p1 + 0.25 * p2 + 0.25 * p3);
}


double MixNormUnif(RealVector &x) {
  double a = -5;
  double b = 0;
  double mean = 20;
  double sigma = 1;
  double p1 = 0.0;

  if (a <= x[0] && x[0] <= b) {
    p1 = 1.0 / (b - a);
  }
  double p2 = 1.0 / sqrt(2 * 3.14159265358979323 * pow(sigma, 2)) *
              exp(-0.5 * (1.0 / pow(sigma, 2)) * pow((x[0] - mean), 2));

  return log(p1 + p2);
}

double bivariatePDF(RealVector &x, double m1, double s1, double m2, double s2) {
  double p1 = 1.0 /(sqrt(2 * 3.1415926535897 * pow(s1, 2))) *
                exp(-0.5 * pow((x[0] - m1), 2) / (pow(s1, 2)));

  double p2 = 1.0 /(sqrt(2 * 3.1415926535897 * pow(s2, 2))) *
                exp(-0.5 * pow((x[1] - m2), 2) / (pow(s2, 2)));

  return p1 * p2;
}

double TwoDMultimodal(RealVector &x) {
  RealVector meanX = {0.2, 0.2, 10.0, 10.0, 0.0};
  RealVector meanY = {0.8, 0.8, -10.0, 10.0, 0.0};
  RealVector sigmaX = {0.05, 0.05, 0.25, 0.25, 0.25};
  RealVector sigmaY = {0.05, 0.05, 1.0, 1.0, 1.0};

  double p1 = bivariatePDF(x, meanX[0], sigmaX[0], meanY[0], sigmaY[0]);
  double p2 = bivariatePDF(x, meanX[1], sigmaX[1], meanY[1], sigmaY[1]);
  /*double p3 = bivariatePDF(x, meanX[2], sigmaX[2], meanY[2], sigmaY[2]);
  double p4 = bivariatePDF(x, meanX[3], sigmaX[3], meanY[3], sigmaY[3]);
  double p5 = bivariatePDF(x, meanX[4], sigmaX[4], meanY[4], sigmaY[4]);
  */
  return log(0.7 * p1 + 0.3 * p2);
}

double schwefelFunc(RealVector &x) {
  double sum = 0.0;

  for (size_t i = 0; i < x.size(); ++i) {
    sum += x[i] * sin(sqrt(abs(x[i])));
  }

  return -(418.9829*x.size() - sum);
}


double MultiDMultimodal(RealVector &x) {
  double p1 = 1.0;
  double p2 = 1.0;
  double p3 = 1.0;

  for (size_t i = 0; i < x.size(); ++i) {
    p1 *= normPDF(x[i], 0.8, 0.05);
    p2 *= normPDF(x[i], 0.2, 0.05);
    if (i % 2 == 0) {
      p3 *= normPDF(x[i], 0.2, 0.05);
    } else {
      p3 *= normPDF(x[i], 0.8, 0.05);
    }
  }
  return log(0.25 * p1 + 0.5 * p2 + 0.25 * p3);
}

double MultiDUnimodal(RealVector &x) {
  double p1 = 0.0;

  for (size_t i = 0; i < x.size(); ++i) {
    p1 += normLogPDF(x[i], 0, 1);
  }
  return p1;
}

double MultiDUnimodalGradLog(double x) {
  return -((x - 0.8) / 0.05);
}

double MultiDBimodal(RealVector &x) {
  double p1 = 1.0;
  double p2 = 1.0;
  for (size_t i = 0; i < x.size(); ++i) {
    p1 *= normPDF(x[i], 0.25, 0.05);
    p2 *= normPDF(x[i], 0.75, 0.05);
  }

  return log(p1 + p2);
}

double MultiDBimodalGradLog(double x, double ndim) {
  double mean1 = 0.25;
  double mean2 = 0.75;
  double sig1 = 0.05;
  double sig2 = 0.05;

  double detSig1 = sqrt(pow(sig1, ndim));
  double detSig2 = sqrt(pow(sig2, ndim));
  double gradA = -(1.0 / sig1) * (x - mean1);
  double gradB = -(1.0 / sig2) * (x - mean2);
  double expA = exp(-0.5 * (x - mean1) * (1.0 / sig1) * (x - mean1));
  double expB = exp(-0.5 * (x - mean2) * (1.0 / sig2) * (x - mean2));

  double num = detSig1 * gradA * expA + detSig2 * gradB * expB;
  double denom = detSig1 * expA + detSig2 * expB;

  return (num / denom);
}
