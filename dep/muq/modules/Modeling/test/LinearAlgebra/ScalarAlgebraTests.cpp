#include <gtest/gtest.h>

#include <memory>

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Modeling;

TEST(ScalarAlgebraTests, Size) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  // should have size 1
  EXPECT_EQ(alg->Size(xd), 1);
  EXPECT_EQ(alg->Size(xf), 1);
  EXPECT_EQ(alg->Size(xi), 1);
  EXPECT_EQ(alg->Size(xui), 1);
}

TEST(ScalarAlgebraTests, SquareRoot) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  double xd=4.0; float xf=3.0; int xi=2; unsigned int xui=8;

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->SquareRoot(xd)), std::sqrt(xd));
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->SquareRoot(xf)), std::sqrt(xf));
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->SquareRoot(xi)), std::sqrt(xi));
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->SquareRoot(xui)), std::sqrt(xui));
}

TEST(ScalarAlgebraTests, LogDeterminate) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  double xd=4.0; float xf=3.0; int xi=2; unsigned int xui=8;

  // test scalar norm
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(xd), std::log(std::fabs(xd)));
  EXPECT_FLOAT_EQ(alg->LogDeterminate(xf), std::log(std::fabs(xf)));
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(xi), std::log(std::fabs(xi)));
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(xui), std::log(xui));
}

TEST(ScalarAlgebraTests, AccessElement) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  // test scalar
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(xd)), xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(xf)), xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(xi)), xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->AccessElement(xui)), xui);
}

TEST(ScalarAlgebraTests, Zero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar zero
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Zero(typeid(double))), 0.0);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Zero(typeid(float))), 0.0);
  EXPECT_EQ(boost::any_cast<int const>(alg->Zero(typeid(int))), 0);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Zero(typeid(unsigned int))), 0);
}

TEST(ScalarAlgebraTests, IsZero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar zero
  EXPECT_TRUE(alg->IsZero((double)0.0));
  EXPECT_FALSE(alg->IsZero((double)1.0));
  EXPECT_TRUE(alg->IsZero((float)0.0));
  EXPECT_FALSE(alg->IsZero((float)-2.3));
  EXPECT_TRUE(alg->IsZero((int)0));
  EXPECT_FALSE(alg->IsZero((int)-6));
  EXPECT_TRUE(alg->IsZero((unsigned int)0));
  EXPECT_FALSE(alg->IsZero((unsigned int)18));
}

TEST(ScalarAlgebraTests, Norm) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  // test scalar norm
  EXPECT_DOUBLE_EQ(alg->Norm(xd), std::fabs(xd));
  EXPECT_FLOAT_EQ(alg->Norm(xf), std::fabs(xf));
  EXPECT_DOUBLE_EQ(alg->Norm(xi), std::fabs(xi));
  EXPECT_DOUBLE_EQ(alg->Norm(xui), std::fabs(xui));
}

TEST(ScalarAlgebraTests, InnerProduct) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xd, xd), xd*xd);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xd, xf), xd*xf);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xd, xi), xd*xi);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xd, xui), xd*xui);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xf, xd), xf*xd);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xf, xf), xf*xf);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xf, xi), xf*xi);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xf, xui), xf*xui);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xi, xd), xi*xd);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xi, xf), xi*xf);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xi, xi), xi*xi);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xi, xui), xi*xui);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xui, xd), xui*xd);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xui, xf), xui*xf);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xui, xi), xui*xi);
  EXPECT_DOUBLE_EQ(alg->InnerProduct(xui, xui), xui*xui);
}

TEST(ScalarAlgebraTests, OuterProduct) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->OuterProduct(xd, xd)), xd*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<double const>(alg->OuterProduct(xd, xf)), xd*xf);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->OuterProduct(xd, xi)), xd*xi);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->OuterProduct(xd, xui)), xd*xui);
  EXPECT_FLOAT_EQ(boost::any_cast<double const>(alg->OuterProduct(xf, xd)), xf*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->OuterProduct(xf, xf)), xf*xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->OuterProduct(xf, xi)), xf*xi);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->OuterProduct(xf, xui)), xf*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->OuterProduct(xi, xd)), xi*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->OuterProduct(xi, xf)), xi*xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->OuterProduct(xi, xi)), xi*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->OuterProduct(xi, xui)), xi*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->OuterProduct(xui, xd)), xui*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->OuterProduct(xui, xf)), xui*xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->OuterProduct(xui, xi)), xui*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->OuterProduct(xui, xui)), xui*xui);
}

TEST(ScalarAlgebraTests, Identity) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // scalars
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Identity(typeid(double))), (double)1.0);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Identity(typeid(float))), (float)1.0);
  EXPECT_EQ(boost::any_cast<int const>(alg->Identity(typeid(int))), (int)1);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Identity(typeid(unsigned int))), (unsigned int)1);
}

TEST(ScalarAlgebraTests, Inverse) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Inverse(xd)), 1.0/xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Inverse(xf)), 1.0/xf);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Inverse(xi)), 1.0/xi);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Inverse(xui)), 1.0/xui);
}

TEST(ScalarAlgebraTests, Add) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xd, xd)), xd+xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xd, xf)), xd+xf);
  EXPECT_EQ(boost::any_cast<double const>(alg->Add(xd, xi)), xd+xi);
  EXPECT_EQ(boost::any_cast<double const>(alg->Add(xd, xui)), xd+xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xf, xd)), xf+xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Add(xf, xf)), xf+xf);
  EXPECT_EQ(boost::any_cast<float const>(alg->Add(xf, xi)), xf+xi);
  EXPECT_EQ(boost::any_cast<float const>(alg->Add(xf, xui)), xf+xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xi, xd)), xi+xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Add(xi, xf)), xi+xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->Add(xi, xi)), xi+xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Add(xi, xui)), xi+xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xui, xd)), xui+xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Add(xui, xf)), xui+xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Add(xui, xi)), xui+xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Add(xui, xui)), xui+xui);
}

TEST(ScalarAlgebraTests, Subtract) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xd, xd)), xd-xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xd, xf)), xd-xf);
  EXPECT_EQ(boost::any_cast<double const>(alg->Subtract(xd, xi)), xd-xi);
  EXPECT_EQ(boost::any_cast<double const>(alg->Subtract(xd, xui)), xd-xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xf, xd)), xf-xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Subtract(xf, xf)), xf-xf);
  EXPECT_EQ(boost::any_cast<float const>(alg->Subtract(xf, xi)), xf-xi);
  EXPECT_EQ(boost::any_cast<float const>(alg->Subtract(xf, xui)), xf-xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xi, xd)), xi-xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Subtract(xi, xf)), xi-xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->Subtract(xi, xi)), xi-xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Subtract(xi, xui)), xi-xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xui, xd)), xui-xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Subtract(xui, xf)), xui-xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Subtract(xui, xi)), xui-xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Subtract(xui, xui)), xui-xui);
}

TEST(ScalarAlgebraTests, Multiply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xd, xd)), xd*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xd, xf)), xd*xf);
  EXPECT_EQ(boost::any_cast<double const>(alg->Multiply(xd, xi)), xd*xi);
  EXPECT_EQ(boost::any_cast<double const>(alg->Multiply(xd, xui)), xd*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xf, xd)), xf*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Multiply(xf, xf)), xf*xf);
  EXPECT_EQ(boost::any_cast<float const>(alg->Multiply(xf, xi)), xf*xi);
  EXPECT_EQ(boost::any_cast<float const>(alg->Multiply(xf, xui)), xf*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xi, xd)), xi*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Multiply(xi, xf)), xi*xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->Multiply(xi, xi)), xi*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Multiply(xi, xui)), xi*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xui, xd)), xui*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Multiply(xui, xf)), xui*xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Multiply(xui, xi)), xui*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Multiply(xui, xui)), xui*xui);
}

TEST(ScalarAlgebraTests, Apply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xd)), xd*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xf)), xf*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xi)), xi*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xui)), xui*xd);

  EXPECT_FLOAT_EQ(boost::any_cast<double const>(alg->Apply(xf, xd)), xd*xf); // use FLOAT_EQ because the precision is lower
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xf, xf)), xf*xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xf, xi)), xi*xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xf, xui)), xui*xf);

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xi, xd)), xd*xi);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xi, xf)), xf*xi);
  EXPECT_EQ(boost::any_cast<int const>(alg->Apply(xi, xi)), xi*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Apply(xi, xui)), xui*xi);

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xui, xd)), xd*xui);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xui, xf)), xf*xui);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Apply(xui, xi)), xi*xui);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Apply(xui, xui)), xui*xui);
}

TEST(ScalarAlgebraTests, ApplyInverse) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xd)), 1.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xf)), xf/xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xi)), xi/xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xui)), xui/xd);

  EXPECT_FLOAT_EQ(boost::any_cast<double const>(alg->ApplyInverse(xf, xd)), xd/xf); // use FLOAT_EQ because the precision is lower
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->ApplyInverse(xf, xf)), 1.0);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->ApplyInverse(xf, xi)), xi/xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->ApplyInverse(xf, xui)), xui/xf);

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xd)), xd/xi);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xf)), xf/xi);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xi)), 1.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xui)), (double)xui/(double)xi);

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xd)), xd/xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xf)), xf/xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xi)), (double)xi/(double)xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xui)), 1);
}
