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

double normPDF(double x, double mean = 0, double delta = 1) {
  return log(1.0/delta) - (pow((x - mean), 2) / (2 * pow(delta, 2)))
        - (0.5 * log(2 * 3.1415926535897));
}

double lineNormal(RealVector &x, double delta = 1) {
  /* Manifold is a straight y = 2 line from x = 0 to x = 10 */
  double mean = 2;

  double pval;
  if (x[0] < 0 || x[0] > 10) {
    pval = pow(10, -50);
  } else {
    pval = normPDF(x[1], mean, delta);
  }

  if (x[0] > 10) {
    std::cout << "PVAL: " << pval << std::endl;
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

  pval = normPDF(r, radius, delta);

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
  pval = normPDF(r, radius, delta);

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
  /* Nonlinear function following xy = 10 */
  double constant = 10;

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
      pval = normPDF(x[1], Amean, delta);
    } else {
      pval = normPDF(x[1], Bmean, delta);
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

  double pval = normPDF(deviate, 0, delta);

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

  double p1 = normPDF(sampleMean, mean, delta);

  double p2 = normPDF(sampleVar, variance, delta);

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

  double p1 = normPDF(deviate1, 0, delta);

  if (x.size() == 6) {
    double deviate2 = pow(x[3], 2) / pow(a2, 2) + pow(x[4], 2) / pow(b2, 2) - x[5];

    double p2 = normPDF(deviate2, 0, delta);

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

  double p1 = normPDF(deviate, 0, delta);

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
    res = res + normPDF(deviate, 0, delta);
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
    res += normPDF(deviate, 0, delta);

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
    p1 = normPDF(x[1], 0, delta);
  }
  if (x[3] > 5 || x[3] < 0) {
    p2 = pow(10, -50);
  } else {
    p2 = normPDF(x[2], -7, delta);
  }
  if (x[5] > 5 || x[5] < 0) {
    p3 = pow(10, -50);
  } else {
    p3 = normPDF(x[4], -4, delta);
  }

  // Q
  double p4dev = pow((x[6] + 1.5), 2) / pow(1.5, 2) + pow((x[7] - 2.5), 2) / pow(2.5, 2) - 1;
  p4 = normPDF(p4dev, 0, delta);

  if (x[8] < -2 || x[8] > 0) {
    p5 = pow(10, -50);
  } else {
    double p5dev = x[9] + (4.0/3.0) * x[8];
    p5 = normPDF(p5dev, 0, delta);
  }

  // T
  if (x[10] < 1 || x[10] > 4) {
    p6 = pow(10, -50);
  } else {
    p6 = normPDF(x[11], 5, delta);
  }
  if (x[13] > 5 || x[13] < 0) {
    p7 = pow(10, -50);
  } else {
    p7 = normPDF(x[12], 2.5, delta);
  }

  // k
  if (x[15] > 5 || x[15] < 0) {
    p8 = pow(10, -50);
  } else {
    p8 = normPDF(x[14], 5, delta);
  }
  if (x[16] < 5 || x[16] > 7.5) {
    p9 = pow(10, -50);
  } else {
    double p9dev = x[17] + (2.0 / 2.5) * x[16] - 6;
    p9 = normPDF(p9dev, 0, delta);
  }

  if (x[18] < 5 || x[18] > 7.5) {
    p10 = pow(10, -50);
  } else {
    double p10dev = x[19] - 0.4 * x[18];
    p10 = normPDF(p10dev, 0, delta);
  }

  return p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10;
}




int main (int argc, const char *argv[]) {
  std::ifstream input_file("mcmcstates_1.dat");
  std::ofstream output_file;
  std::string line, token;
  std::stringstream iss;
  std::ifstream delta_file("delta.dat");
  double delta;

  delta_file >> delta;

  if (argc > 1) {
    std::string deltaStr(argv[1]);
    delta = std::stod(deltaStr);
    output_file.open("tmcmc_lp.dat");
  } else {
    output_file.open("tmcmc_ll.dat");
  }

  while (std::getline(input_file, line)) {
    RealVector x;
    iss << line;
    double d;
    while (iss >> d) {
      x.push_back(d);
    }

    output_file << std::setprecision(18) << nonlinearNormal(x, delta) << "\n";
    x.clear();
    iss.clear();
  }
  //std::cout << std::setprecision(18) << "DELTA: " << delta << std::endl;
  output_file.close();
  return 0;
}
