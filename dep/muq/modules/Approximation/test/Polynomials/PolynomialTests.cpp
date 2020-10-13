#include <gtest/gtest.h>

#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"
#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/Laguerre.h"
#include "MUQ/Approximation/Polynomials/Jacobi.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Approximation;

TEST(Polynomial, Monomial) {
  // create a monomial object
  auto mono = std::make_shared<Monomial>();
  EXPECT_EQ(mono->numInputs, 2);
  EXPECT_EQ(mono->numOutputs, 1);

  // a point to evaluate the monomial
  const double x = 2.0;

  // the order of the polynomial
  const unsigned int max_order = 200;

  Eigen::VectorXd allEvals = mono->EvaluateAllTerms(max_order, x);

  for( unsigned int p=0; p<=max_order; ++p ) {

    EXPECT_DOUBLE_EQ(allEvals(p), std::pow(x, p));

    // evaluate the monomial
    const std::vector<boost::any>& result = mono->Evaluate(p, x);
    EXPECT_EQ(result.size(), 1);
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(result[0]), std::pow(x, p));

    // Get the first derivative
    double deriv = mono->DerivativeEvaluate(p,1,x);
    EXPECT_DOUBLE_EQ(double(p)*std::pow(x,p-1.0), deriv);

    // Get the second derivative
    deriv = mono->DerivativeEvaluate(p,2,x);
    EXPECT_DOUBLE_EQ(double(p)*double(p-1)*std::pow(x,p-2.0), deriv);

    // Third derivative
    deriv = mono->DerivativeEvaluate(p,3,x);
    EXPECT_DOUBLE_EQ(double(p)*double(p-1)*double(p-2)*std::pow(x,p-3.0), deriv);

  }

}

TEST(Polynomial, PhysicistHermite) {
  // create a Hermite object
  auto hermite = std::make_shared<PhysicistHermite>();


  // Evaluations
  EXPECT_DOUBLE_EQ(1.0, boost::any_cast<double const>(hermite->Evaluate((unsigned int)0, 0.4) [0]));
  EXPECT_DOUBLE_EQ(0.6, boost::any_cast<double const>(hermite->Evaluate((unsigned int)1, 0.3) [0]));
  EXPECT_DOUBLE_EQ(4.0*std::pow(0.6, 2.0)-2.0, boost::any_cast<double const>(hermite->Evaluate((unsigned int)2, 0.6) [0]));
  EXPECT_DOUBLE_EQ(33.6235290625, boost::any_cast<double const>(hermite->Evaluate((unsigned int)5, 0.325) [0]));
  EXPECT_NEAR(6219.5581337600015, boost::any_cast<double const>(hermite->Evaluate((unsigned int)8, 1.6) [0]), 2.0e-12);
  EXPECT_NEAR(6.075804453410837e11, boost::any_cast<double const>(hermite->Evaluate((unsigned int)20, -0.845) [0]), 3.0e-4);

  // First derivatives
  const double x = 0.23;
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(2.0, hermite->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(8.0*x, hermite->DerivativeEvaluate(2,1,x));
  EXPECT_DOUBLE_EQ(24.0*x*x - 12.0, hermite->DerivativeEvaluate(3,1,x));
  EXPECT_DOUBLE_EQ(64.0*std::pow(x,3.0) -96.0*x, hermite->DerivativeEvaluate(4,1,x));

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(8.0, hermite->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(48.0*x, hermite->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(192.0*std::pow(x,2.0) - 96.0, hermite->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(48.0, hermite->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(384.0*x, hermite->DerivativeEvaluate(4,3,x));

  // Check the normalization constant
  EXPECT_DOUBLE_EQ(sqrt(M_PI), hermite->Normalization(0));
  EXPECT_DOUBLE_EQ(sqrt(M_PI)*2.0, hermite->Normalization(1));
  EXPECT_DOUBLE_EQ(sqrt(M_PI)*8.0, hermite->Normalization(2));
  EXPECT_DOUBLE_EQ(sqrt(M_PI)*48.0, hermite->Normalization(3));

}

TEST(Polynomial, ProbabilistHermite) {
  // create a Hermite object
  auto hermite = std::make_shared<ProbabilistHermite>();

  const double x = 1.32;
  EXPECT_DOUBLE_EQ(1.0, hermite->BasisEvaluate(0, x));
  EXPECT_DOUBLE_EQ(x, hermite->BasisEvaluate(1, x));
  EXPECT_DOUBLE_EQ(x*x-1.0, hermite->BasisEvaluate(2, x));
  EXPECT_DOUBLE_EQ(x*x*x - 3.0*x, hermite->BasisEvaluate(3, x));
  EXPECT_DOUBLE_EQ(x*x*x*x - 6.0*x*x + 3.0, hermite->BasisEvaluate(4, x));
  EXPECT_NEAR(std::pow(x,5) - 10*std::pow(x,3) + 15.0*x, hermite->BasisEvaluate(5, x), 1e-10);

  // First derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(1.0, hermite->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(2.0*x, hermite->DerivativeEvaluate(2,1,x));
  EXPECT_DOUBLE_EQ(3.0*x*x - 3.0, hermite->DerivativeEvaluate(3,1,x));
  EXPECT_DOUBLE_EQ(4.0*std::pow(x,3.0) - 12.0*x, hermite->DerivativeEvaluate(4,1,x));

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(2.0, hermite->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(6.0*x, hermite->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(12.0*std::pow(x,2.0) - 12.0, hermite->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, hermite->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(6.0, hermite->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(24.0*x, hermite->DerivativeEvaluate(4,3,x));

  // Check the normalization constant
  EXPECT_DOUBLE_EQ(sqrt(2.0*M_PI), hermite->Normalization(0));
  EXPECT_DOUBLE_EQ(sqrt(2.0*M_PI), hermite->Normalization(1));
  EXPECT_DOUBLE_EQ(sqrt(2.0*M_PI)*2.0, hermite->Normalization(2));
  EXPECT_DOUBLE_EQ(sqrt(2.0*M_PI)*6.0, hermite->Normalization(3));

}

TEST(Polynomial, Legendre) {

  // create a Legendre object
  auto legendre = std::make_shared<Legendre>();

  // A sample of points tested against mathematica
  EXPECT_DOUBLE_EQ(1.0, boost::any_cast<double const>(legendre->Evaluate((unsigned int)0, 0.3) [0]));
  EXPECT_DOUBLE_EQ(0.3, boost::any_cast<double const>(legendre->Evaluate((unsigned int)1, 0.3) [0]));
  EXPECT_DOUBLE_EQ(0.5*(3.0*std::pow(0.3, 2.0)-1.0), boost::any_cast<double const>(legendre->Evaluate((unsigned int)2, 0.3) [0]));
  EXPECT_DOUBLE_EQ(0.3375579333496094, boost::any_cast<double const>(legendre->Evaluate((unsigned int)5, 0.325) [0]));
  EXPECT_NEAR(-0.05346106275520913, boost::any_cast<double const>(legendre->Evaluate((unsigned int)20, -0.845) [0]), 1.0e-14);
  EXPECT_NEAR(-0.1119514835092105, boost::any_cast<double const>(legendre->Evaluate((unsigned int)50, 0.1264) [0]), 1e-14);
  EXPECT_NEAR(-0.001892916076323403, boost::any_cast<double const>(legendre->Evaluate((unsigned int)200, -0.3598) [0]), 1e-14);
  EXPECT_NEAR(0.01954143166718206, boost::any_cast<double const>(legendre->Evaluate((unsigned int)1000, 0.4587) [0]), 1e-14);

  // First derivatives
  const double x = 0.23;
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(1.0, legendre->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(3.0*x, legendre->DerivativeEvaluate(2,1,x));
  EXPECT_DOUBLE_EQ(7.5*x*x - 1.5, legendre->DerivativeEvaluate(3,1,x));
  EXPECT_DOUBLE_EQ(17.5*std::pow(x,3.0) - 7.5*x, legendre->DerivativeEvaluate(4,1,x));

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(3.0, legendre->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(15.0*x, legendre->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(52.5*std::pow(x,2.0) - 7.5, legendre->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, legendre->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(15.0, legendre->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(105.0*x, legendre->DerivativeEvaluate(4,3,x));

  // Check the normalization constant
  EXPECT_DOUBLE_EQ(2.0, legendre->Normalization(0));
  EXPECT_DOUBLE_EQ(2.0/3.0, legendre->Normalization(1));
  EXPECT_DOUBLE_EQ(2.0/5.0, legendre->Normalization(2));
  EXPECT_DOUBLE_EQ(2.0/7.0, legendre->Normalization(3));
}



TEST(Polynomial, Laguerre) {

  // create a Legendre object
  auto poly = std::make_shared<Laguerre>(0.0);

  const double x = 1.32;
  EXPECT_DOUBLE_EQ(1.0, poly->BasisEvaluate(0, x));
  EXPECT_DOUBLE_EQ(-x+1.0, poly->BasisEvaluate(1, x));
  EXPECT_DOUBLE_EQ(0.5*(x*x-4.0*x+2.0), poly->BasisEvaluate(2, x));
  EXPECT_DOUBLE_EQ((1.0/6.0)*(-x*x*x + 9.0*x*x - 18.0*x + 6), poly->BasisEvaluate(3, x));
  EXPECT_DOUBLE_EQ((1.0/24.0)*(std::pow(x,4) - 16.0*std::pow(x,3.0) + 72.0*std::pow(x,2.0) - 96.0*x + 24.0), poly->BasisEvaluate(4, x));
  EXPECT_NEAR((1.0/120)*(-1.0*std::pow(x,5) +25.0*std::pow(x,4) - 200.0*std::pow(x,3.0) + 600.0*std::pow(x,2.0) - 600.0*x + 120.0), poly->BasisEvaluate(5, x), 1e-10);

  // First derivatives
  EXPECT_DOUBLE_EQ(0.0, poly->DerivativeEvaluate(0,1,x));
  EXPECT_DOUBLE_EQ(-1.0, poly->DerivativeEvaluate(1,1,x));
  EXPECT_DOUBLE_EQ(x - 2.0, poly->DerivativeEvaluate(2,1,x));
  EXPECT_NEAR(-0.5*x*x +3.0*x -3.0, poly->DerivativeEvaluate(3,1,x), 1e-10);
  EXPECT_NEAR((1.0/6.0)*std::pow(x,3.0) - 2.0*x*x + 6.0*x - 4.0, poly->DerivativeEvaluate(4,1,x), 1e-10);

  // Second derivatives
  EXPECT_DOUBLE_EQ(0.0, poly->DerivativeEvaluate(0,2,x));
  EXPECT_DOUBLE_EQ(0.0, poly->DerivativeEvaluate(1,2,x));
  EXPECT_DOUBLE_EQ(1.0, poly->DerivativeEvaluate(2,2,x));
  EXPECT_DOUBLE_EQ(-x + 3.0, poly->DerivativeEvaluate(3,2,x));
  EXPECT_DOUBLE_EQ(0.5*x*x -4.0*x + 6.0, poly->DerivativeEvaluate(4,2,x));

  // Third derivatives
  EXPECT_DOUBLE_EQ(0.0, poly->DerivativeEvaluate(0,3,x));
  EXPECT_DOUBLE_EQ(0.0, poly->DerivativeEvaluate(1,3,x));
  EXPECT_DOUBLE_EQ(0.0, poly->DerivativeEvaluate(2,3,x));
  EXPECT_DOUBLE_EQ(-1.0, poly->DerivativeEvaluate(3,3,x));
  EXPECT_DOUBLE_EQ(x-4.0, poly->DerivativeEvaluate(4,3,x));

  // Check the normalization constant
  EXPECT_DOUBLE_EQ(1.0, poly->Normalization(0));
  EXPECT_DOUBLE_EQ(1.0, poly->Normalization(1));
  EXPECT_DOUBLE_EQ(1.0, poly->Normalization(2));
  EXPECT_DOUBLE_EQ(1.0, poly->Normalization(3));
}


TEST(Polynomial, Jacobi) {

  // Compare the jacobi with a=b=0 to a Legendre (which should be the same thing)
  const double a = 0.0;
  const double b = 0.0;
  auto poly = std::make_shared<Jacobi>(a,b);
  auto legendre = std::make_shared<Legendre>();

  const double x = 0.32;
  for(unsigned i=0; i<5; ++i){
      EXPECT_DOUBLE_EQ(legendre->BasisEvaluate(i,x), poly->BasisEvaluate(i,x));
      EXPECT_DOUBLE_EQ(legendre->DerivativeEvaluate(i,1,x), poly->DerivativeEvaluate(i,1,x));
      EXPECT_DOUBLE_EQ(legendre->DerivativeEvaluate(i,2,x), poly->DerivativeEvaluate(i,2,x));
      EXPECT_DOUBLE_EQ(legendre->DerivativeEvaluate(i,3,x), poly->DerivativeEvaluate(i,3,x));
      EXPECT_DOUBLE_EQ(legendre->DerivativeEvaluate(i,4,x), poly->DerivativeEvaluate(i,4,x));
  }

}

TEST(Polynomial, Factory)
{

    std::shared_ptr<IndexedScalarBasis> monomial = IndexedScalarBasis::Construct("Monomial");
    EXPECT_TRUE(std::dynamic_pointer_cast<Monomial>(monomial));

    std::shared_ptr<IndexedScalarBasis> hermite1 = IndexedScalarBasis::Construct("PhysicistHermite");
    EXPECT_TRUE(std::dynamic_pointer_cast<PhysicistHermite>(hermite1));

    std::shared_ptr<IndexedScalarBasis> hermite2 = IndexedScalarBasis::Construct("ProbabilistHermite");
    EXPECT_TRUE(std::dynamic_pointer_cast<ProbabilistHermite>(hermite2));

    std::shared_ptr<IndexedScalarBasis> legendre = IndexedScalarBasis::Construct("Legendre");
    EXPECT_TRUE(std::dynamic_pointer_cast<Legendre>(legendre));

    std::shared_ptr<IndexedScalarBasis> laguerre = IndexedScalarBasis::Construct("Laguerre");
    EXPECT_TRUE(std::dynamic_pointer_cast<Laguerre>(laguerre));

    std::shared_ptr<IndexedScalarBasis> jacobi = IndexedScalarBasis::Construct("Jacobi");
    EXPECT_TRUE(std::dynamic_pointer_cast<Jacobi>(jacobi));

    EXPECT_THROW(IndexedScalarBasis::Construct("CowPoly"), muq::NotRegisteredError);

}
