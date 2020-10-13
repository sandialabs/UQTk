
#ifndef _RandomGenerator_h
#define _RandomGenerator_h

#include <random>

#include <Eigen/Core>
#include <array>
#include <mutex>


namespace muq {
namespace Utilities {

/**
 * Call this class to temporarily set the seed. The previous state is restored when the object
 * is destructed. This is a fairly low quality way to do this, but is fine for testing.
 */
class RandomGeneratorTemporarySetSeed
{
public:
	RandomGeneratorTemporarySetSeed(int seed);

	virtual ~RandomGeneratorTemporarySetSeed();
private:
	std::mt19937 oldGenerator;
};



/** The RandGen generator class is a wrapper around the std::random number generation library.  The point of this class
 *  is to provide a super easy to use interface for generating Random numbers of canonical scalar distributions (e.g., Uniform, Gaussian, Gamma, Uniform Integer, etc...
 */
class RandomGenerator {
public:

  typedef std::mt19937 GeneratorType;

  /** Get a standard Gaussian distribution sample. */
  static double          GetNormal();
  static Eigen::MatrixXd GetNormal(int rows){return GetNormal(rows,1);};
  static Eigen::MatrixXd GetNormal(int rows, int cols);

  /** Get a uniformly distributed real sample in [0,1). */
  static double          GetUniform();
  static Eigen::MatrixXd GetUniform(int rows){return GetUniform(rows,1);};
  static Eigen::MatrixXd GetUniform(int rows, int cols);

  static double          GetGamma(double const alpha, double const beta);
  static Eigen::MatrixXd GetGamma(double const alpha, double const beta, int rows){return GetGamma(alpha, beta, rows, 1);};
  static Eigen::MatrixXd GetGamma(double const alpha, double const beta, int rows, int cols);

  /** Get a random integer, distributed uniformly, between the bounds. */
  static int             GetUniformInt(int lb, int ub);
  static Eigen::MatrixXi GetUniformInt(int lb, int ub, int rows){return GetUniformInt(lb, ub, rows, 1, true);};
  static Eigen::MatrixXi GetUniformInt(int lb, int ub, int rows, bool unique){return GetUniformInt(lb, ub, rows, 1, unique);};
  static Eigen::MatrixXi GetUniformInt(int lb, int ub, int rows, int cols){return GetUniformInt(lb, ub, rows, cols, true);};
  static Eigen::MatrixXi GetUniformInt(int lb, int ub, int rows, int cols, bool unique);


  /** Sample a discrete distribution with replacement. */
  static int GetDiscrete(Eigen::VectorXd const& discProbs);
  static Eigen::MatrixXi GetDiscrete(Eigen::VectorXd const& discProbs, int rows, int cols=1);

  /// Set the seed for the random number generator. This is a fairly low quality way to do this, but is fine for testing.
  static void SetSeed(int seedval);

  ///Store a copy of the generator, designed for use with SetGenerator. Be very careful in using this, or you could wreck the library's access to random numbers.
  static GeneratorType CopyGenerator();

  ///Set the state of the generator, designed for use with CopyGenerator. Be very careful in using this, or you could wreck the library's access to random numbers.
  static void SetGenerator(GeneratorType);

private:

  friend RandomGeneratorTemporarySetSeed;
  static GeneratorType & GetGenerator();
};


class SeedGenerator {
public:

  SeedGenerator();
  SeedGenerator(const std::array<std::seed_seq::result_type, RandomGenerator::GeneratorType::state_size> &seed_data);
  //draw a seed for the entire state, as per: http://stackoverflow.com/questions/15509270/does-stdmt19937-require-warmup
  std::seed_seq seed_seq;
  //virtual ~SeedGenerator() = default;
};



} // namespace GeneralUtilties
} //namespace muq

#endif // ifndef _RandomGenerator_h
