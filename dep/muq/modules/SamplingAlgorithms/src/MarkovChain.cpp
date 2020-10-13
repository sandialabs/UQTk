#include "MUQ/SamplingAlgorithms/MarkovChain.h"

#include <unsupported/Eigen/FFT>

using namespace muq::SamplingAlgorithms;

Eigen::VectorXd MarkovChain::ESS(int blockDim) const
{
  Eigen::MatrixXd sampMat = AsMatrix(blockDim);

  Eigen::VectorXd ess(sampMat.rows());
  for(int row=0; row<sampMat.rows(); ++row)
    ess(row) = SingleComponentESS(sampMat.row(row));

  return ess;
}

double MarkovChain::SingleComponentESS(Eigen::Ref<const Eigen::VectorXd> const& trace)
{
  int size = trace.size();
  // must have a positive number of samples
  assert( size>=0 );

  // ESS needs at least 2 samples; just return ESS=0.0 if there aren't enough
  if( size<2 ) {
    return 0.0;
  }

  double traceMean = trace.mean();

  Eigen::FFT<double> fft;

  int tmax    = floor(0.5*size);
  double Stau = 1.5;

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> freqVec;
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                1> timeVec = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>::Zero(size + tmax);
  for (int i = 0; i < size; ++i) {
    timeVec(i) = std::complex<double>(trace(i) - traceMean, 0.0);
  }


  fft.fwd(freqVec, timeVec);


  // compute out1*conj(out2) and store in out1 (x+yi)*(x-yi)
  double real;

  for (int i = 0; i < size + tmax; ++i) {
    real       = freqVec(i).real() * freqVec(i).real() + freqVec(i).imag() * freqVec(i).imag();
    freqVec(i) = std::complex<double>(real, 0.0);
  }

  // now compute the inverse fft to get the autocorrelation (stored in timeVec)
  fft.inv(timeVec, freqVec);

  for (int i = 0; i < tmax + 1; ++i) {
    timeVec(i) = std::complex<double>(timeVec(i).real() / double(size - i), 0.0);
  }

  // the following loop uses ideas from "Monte Carlo errors with less errors." by Ulli Wolff to figure out how far we
  // need to integrate
  //int MaxLag = 0;
  double Gint = 0;
  int    Wopt = 0;
  for (int i = 1; i < tmax + 1; ++i) {
    Gint += timeVec(i).real() / timeVec(0).real(); // in1[i][0] /= scale;

    double tauW;
    if (Gint <= 0) {
      tauW = 1.0e-15;
    } else {
      tauW = Stau / log((Gint + 1) / Gint);
    }
    double gW = exp(-double(i) / tauW) - tauW / sqrt(double(i) * size);

    if (gW < 0) {
      Wopt = i;
      tmax = std::min(tmax, 2 * i);
      break;
    }
  }


  // correct for bias
  double CFbbopt = timeVec(0).real();
  for (int i = 0; i < Wopt + 1; ++i) {
    CFbbopt += 2 * timeVec(i + 1).real();
  }

  CFbbopt = CFbbopt / size;

  for (int i = 0; i < tmax + 1; ++i) {
    timeVec(i) += std::complex<double>(CFbbopt, 0.0);
  }

  // compute the normalized autocorrelation
  double scale = timeVec(0).real();
  for (int i = 0; i < Wopt; ++i) {
    timeVec(i) = std::complex<double>(timeVec(i).real() / scale, 0.0);
  }

  double tauint = 0;
  for (int i = 0; i < Wopt; ++i) {
    tauint += timeVec(i).real();
  }

  tauint -= 0.5;

  // return the effective sample size
  double frac = 1.0 / (2.0 * tauint);
  frac = fmin(1.0, frac);
  return size * frac;
}
