
#include "MUQ/Utilities/RandomGenerator.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <array>


using namespace muq::Utilities;
using namespace std;

array<seed_seq::result_type, RandomGenerator::GeneratorType::state_size> UrandomRead()
{
  // use a random device instead of /dev/urandom (may result in same behavior if not on windows)
  random_device rd;

  //how much data the generator needs
  array<seed_seq::result_type, RandomGenerator::GeneratorType::state_size> seed_data;

  //generate all the randomness
  for (int i = 0; i < seed_data.size(); ++i) {
    seed_data.at(i) = rd();
  }

  return seed_data;
}

SeedGenerator::SeedGenerator(const array<seed_seq::result_type,
                             RandomGenerator::GeneratorType::state_size>& seed_data) : seed_seq(begin(seed_data),
                                                                                                end(seed_data))
{}

SeedGenerator::SeedGenerator() : SeedGenerator(UrandomRead())
{}

double RandomGenerator::GetNormal()
{
  static std::normal_distribution<> Gauss_RNG(0, 1);

  Gauss_RNG.reset();
  return Gauss_RNG(GetGenerator());
}

double RandomGenerator::GetUniform()
{
  static  std::uniform_real_distribution<> Uniform_RNG(0, 1);

  return Uniform_RNG(GetGenerator());
}

double RandomGenerator::GetGamma(double const alpha, double const beta)
{
  std::gamma_distribution<> Gamma_RNG(alpha, beta);
  return Gamma_RNG(GetGenerator());
}

Eigen::MatrixXd RandomGenerator::GetGamma(double const alpha, double const beta, int rows, int cols)
{

    Eigen::MatrixXd output(rows,cols);

    for(int j=0; j<cols; ++j)
    {
        for(int i=0; i<rows; ++i)
        {
            output(i,j) = GetGamma(alpha,beta);
        }
    }

    return output;
}



int RandomGenerator::GetUniformInt(int lb, int ub)
{
  assert(ub >= lb);

  static  std::uniform_real_distribution<> Uniform_RNG(0, 1);
  return round(Uniform_RNG(GetGenerator()) * (ub - lb) + lb);
}

/** Generate N unique random integers from the range lb, ub */
Eigen::MatrixXi RandomGenerator::GetUniformInt(int lb, int ub, int rows, int cols, bool unique)
{
    const int N = rows*cols;

    assert(ub >= lb);

    const int maxN = ub - lb + 1;
    if(unique){
        assert(N <= maxN);
    }

    if(unique)
    {
        // fill a set with all numbers between lb and ub
        vector<int> allInts;
        allInts.reserve(ub - lb + 1);
        for (int i = lb; i <= ub; ++i) {
            allInts.push_back(i);
        }

        // create the output vector and fill it
        Eigen::VectorXi output(N);
        for (int i = 0; i < N; ++i) {
            // generate a random index into allInts
            int ind = GetUniformInt(i, maxN - 1);
            std::swap(allInts[i], allInts[ind]);
        }

        return Eigen::Map<Eigen::MatrixXi>(&allInts[0], rows, cols);
    }
    else
    {
        Eigen::MatrixXi output(rows,cols);
        for(int j=0; j<cols; ++j){
            for(int i=0; i<rows; ++i)
                output(i,j) = GetUniformInt(lb,ub);
        }

        return output;
    }
}

Eigen::MatrixXd RandomGenerator::GetUniform(int rows, int cols)
{
  Eigen::MatrixXd result(rows,cols);

  for(int j=0; j<cols; ++j){
      for (int i = 0; i < rows; ++i) {
          result(i,j) = GetUniform();
      }
  }

  return result;
}


Eigen::MatrixXd RandomGenerator::GetNormal(int const m, int const n)
{
  Eigen::MatrixXd result(m, n);

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      result(i, j) = GetNormal();
    }
  }
  return result;
}

void RandomGenerator::SetSeed(int seedval)
{
  GetGenerator().seed(seedval);
}

RandomGenerator::GeneratorType RandomGenerator::CopyGenerator()
{
  return GetGenerator();
}

void RandomGenerator::SetGenerator(GeneratorType state)
{
  GetGenerator() = state;
}

int RandomGenerator::GetDiscrete(Eigen::VectorXd const& discProbs)
{
  double u = RandomGenerator::GetUniform();
  double cumSum = 0.0;

  for(int i=0; i<discProbs.size(); ++i){
    cumSum += discProbs(i);
    if(u<cumSum)
      return  i;
  }

  return 0;
}

Eigen::MatrixXi RandomGenerator::GetDiscrete(Eigen::VectorXd const& discProbs, int rows, int cols)
{
  std::vector<int> indices(discProbs.size());
  for(int i=0; i<discProbs.size(); ++i)
    indices[i] = i;

  // Sort in descending order so that we're more likely to terminate early in the cumulative sum loop below
  std::sort(indices.begin(), indices.end(), [discProbs](int const& a, int const& b) -> bool{ return discProbs(a) > discProbs(b); });

  Eigen::MatrixXi output(rows,cols);
  for(int j=0; j<cols; ++j){
    for(int i=0; i<rows; ++i){

      double cumSum = 0;
      double u = RandomGenerator::GetUniform();

      for(int k=0; k<discProbs.size(); ++k){
        cumSum += discProbs(indices.at(k));
        if(u<cumSum){
          output(i,j) = indices.at(k);
          break;
        }
      }
    }
  }

  return output;
}


RandomGenerator::GeneratorType& RandomGenerator::GetGenerator()
{

#if defined(__has_feature) && __has_feature(cxx_thread_local)
#  define MUQ_NATIVE_TLS thread_local
#endif

#if defined(MUQ_NATIVE_TLS)
  static MUQ_NATIVE_TLS SeedGenerator seedGen;
  static MUQ_NATIVE_TLS RandomGenerator::GeneratorType BaseGenerator(seedGen.seed_seq);
#else
  static SeedGenerator seedGen;
  static RandomGenerator::GeneratorType BaseGenerator(seedGen.seed_seq);
#endif

  return BaseGenerator;
}


RandomGeneratorTemporarySetSeed::RandomGeneratorTemporarySetSeed(int seed)
{
  oldGenerator = RandomGenerator::GetGenerator();
  RandomGenerator::SetSeed(seed);
}

RandomGeneratorTemporarySetSeed::~RandomGeneratorTemporarySetSeed()
{
  RandomGenerator::GetGenerator() = oldGenerator;
}
