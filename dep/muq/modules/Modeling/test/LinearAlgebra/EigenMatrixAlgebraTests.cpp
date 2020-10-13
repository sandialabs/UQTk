#include <gtest/gtest.h>

#include <memory>

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Modeling;

TEST(EigenMatrixAlgebraTests, Size) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // check for Eigen::Matrix's
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Constant(2.0);
  EXPECT_EQ(alg->Size(mat2d), 4); // total number of entries
  EXPECT_EQ(alg->Size(mat2d, 0), 2); // number of rows
  EXPECT_EQ(alg->Size(mat2d, 1), 2); // number of cols
  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Constant(2.0);
  EXPECT_EQ(alg->Size(mat3d), 9); // total number of entries
  EXPECT_EQ(alg->Size(mat3d, 0), 3); // number of rows
  EXPECT_EQ(alg->Size(mat3d, 1), 3); // number of cols
  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Constant(2.0);
  EXPECT_EQ(alg->Size(mat4d), 16); // total number of entries
  EXPECT_EQ(alg->Size(mat4d, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(mat4d, 1), 4); // number of cols
  const Eigen::MatrixXd matXd = Eigen::MatrixXd::Constant(4, 5, 2.0);
  EXPECT_EQ(alg->Size(matXd), 20); // total number of entries
  EXPECT_EQ(alg->Size(matXd, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(matXd, 1), 5); // number of cols

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Constant(2.0);
  EXPECT_EQ(alg->Size(mat2f), 4); // total number of entries
  EXPECT_EQ(alg->Size(mat2f, 0), 2); // number of rows
  EXPECT_EQ(alg->Size(mat2f, 1), 2); // number of cols
  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Constant(2.0);
  EXPECT_EQ(alg->Size(mat3f), 9); // total number of entries
  EXPECT_EQ(alg->Size(mat3f, 0), 3); // number of rows
  EXPECT_EQ(alg->Size(mat3f, 1), 3); // number of cols
  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Constant(2.0);
  EXPECT_EQ(alg->Size(mat4f), 16); // total number of entries
  EXPECT_EQ(alg->Size(mat4f, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(mat4f, 1), 4); // number of cols
  const Eigen::MatrixXf matXf = Eigen::MatrixXf::Constant(4, 5, 2.0);
  EXPECT_EQ(alg->Size(matXf), 20); // total number of entries
  EXPECT_EQ(alg->Size(matXf, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(matXf, 1), 5); // number of cols

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Constant(2.0);
  EXPECT_EQ(alg->Size(mat2i), 4); // total number of entries
  EXPECT_EQ(alg->Size(mat2i, 0), 2); // number of rows
  EXPECT_EQ(alg->Size(mat2i, 1), 2); // number of cols
  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Constant(2.0);
  EXPECT_EQ(alg->Size(mat3i), 9); // total number of entries
  EXPECT_EQ(alg->Size(mat3i, 0), 3); // number of rows
  EXPECT_EQ(alg->Size(mat3i, 1), 3); // number of cols
  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Constant(2.0);
  EXPECT_EQ(alg->Size(mat4i), 16); // total number of entries
  EXPECT_EQ(alg->Size(mat4i, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(mat4i, 1), 4); // number of cols
  const Eigen::MatrixXi matXi = Eigen::MatrixXi::Constant(4, 5, 2.0);
  EXPECT_EQ(alg->Size(matXi), 20); // total number of entries
  EXPECT_EQ(alg->Size(matXi, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(matXi, 1), 5); // number of cols
}

TEST(EigenMatrixAlgebraTests, AccessElement) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test Eigen::Matrix's
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  for( unsigned int i=0; i<2; ++i ) {
    for( unsigned int j=0; j<2; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(mat2d, i, j)), mat2d(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(mat2f, i, j)), mat2f(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mat2i, i, j)), mat2i(i, j));
    }
  }

  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  for( unsigned int i=0; i<3; ++i ) {
    for( unsigned int j=0; j<3; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(mat3d, i, j)), mat3d(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(mat3f, i, j)), mat3f(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mat3i, i, j)), mat3i(i, j));
    }
  }

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  for( unsigned int i=0; i<4; ++i ) {
    for( unsigned int j=0; j<4; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(mat4d, i, j)), mat4d(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(mat4f, i, j)), mat4f(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mat4i, i, j)), mat4i(i, j));
    }
  }

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(5, 8);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(5, 8);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5, 8);
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<8; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(matd, i, j)), matd(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(matf, i, j)), matf(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mati, i, j)), mati(i, j));
    }
  }
}

TEST(EigenMatrixAlgebraTests, Zero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Matrix zero
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix2d const&>(alg->Zero(typeid(Eigen::Matrix2d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix2f const&>(alg->Zero(typeid(Eigen::Matrix2f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix2i const&>(alg->Zero(typeid(Eigen::Matrix2i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix3d const&>(alg->Zero(typeid(Eigen::Matrix3d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix3f const&>(alg->Zero(typeid(Eigen::Matrix3f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix3i const&>(alg->Zero(typeid(Eigen::Matrix3i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix4d const&>(alg->Zero(typeid(Eigen::Matrix4d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix4f const&>(alg->Zero(typeid(Eigen::Matrix4f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix4i const&>(alg->Zero(typeid(Eigen::Matrix4i))).norm(), 0.0);

  const Eigen::MatrixXd matd = boost::any_cast<Eigen::MatrixXd const&>(alg->Zero(typeid(Eigen::MatrixXd), 13, 5));
  const Eigen::MatrixXf matf = boost::any_cast<Eigen::MatrixXf const&>(alg->Zero(typeid(Eigen::MatrixXf), 13, 5));
  const Eigen::MatrixXi mati = boost::any_cast<Eigen::MatrixXi const&>(alg->Zero(typeid(Eigen::MatrixXi), 13, 5));
  EXPECT_DOUBLE_EQ(matd.norm(), 0.0);
  EXPECT_EQ(matd.rows(), 13);
  EXPECT_EQ(matd.cols(), 5);
  EXPECT_DOUBLE_EQ(matf.norm(), 0.0);
  EXPECT_EQ(matf.rows(), 13);
  EXPECT_EQ(matf.cols(), 5);
  EXPECT_DOUBLE_EQ(mati.norm(), 0.0);
  EXPECT_EQ(mati.rows(), 13);
  EXPECT_EQ(mati.cols(), 5);
}

TEST(EigenMatrixAlgebraTests, IsZero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  EXPECT_TRUE(alg->IsZero((Eigen::Matrix2d)Eigen::Matrix2d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix2d)Eigen::Matrix2d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix2f)Eigen::Matrix2f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix2f)Eigen::Matrix2f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix2i)Eigen::Matrix2i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix2i)Eigen::Matrix2i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Matrix3d)Eigen::Matrix3d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix3d)Eigen::Matrix3d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix3f)Eigen::Matrix3f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix3f)Eigen::Matrix3f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix3i)Eigen::Matrix3i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix3i)Eigen::Matrix3i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Matrix4d)Eigen::Matrix4d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix4d)Eigen::Matrix4d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix4f)Eigen::Matrix4f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix4f)Eigen::Matrix4f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix4i)Eigen::Matrix4i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix4i)Eigen::Matrix4i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::MatrixXd)Eigen::MatrixXd::Zero(21, 5)));
  EXPECT_FALSE(alg->IsZero((Eigen::MatrixXd)Eigen::MatrixXd::Random(38, 13)));
  EXPECT_TRUE(alg->IsZero((Eigen::MatrixXf)Eigen::MatrixXf::Zero(19, 58)));
  EXPECT_FALSE(alg->IsZero((Eigen::MatrixXf)Eigen::MatrixXf::Random(22, 87)));
  EXPECT_TRUE(alg->IsZero((Eigen::MatrixXi)Eigen::MatrixXi::Zero(36, 23)));
  EXPECT_FALSE(alg->IsZero((Eigen::MatrixXi)Eigen::MatrixXi::Random(5, 9)));
}

TEST(EigenMatrixAlgebraTests, Identity) {
    auto alg = std::shared_ptr<AnyAlgebra>();

    // Eigen::Matrix
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d>(alg->Identity(typeid(Eigen::Matrix2d))).array()==Eigen::Matrix2d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f>(alg->Identity(typeid(Eigen::Matrix2f))).array()==Eigen::Matrix2f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i>(alg->Identity(typeid(Eigen::Matrix2i))).array()==Eigen::Matrix2i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d>(alg->Identity(typeid(Eigen::Matrix3d))).array()==Eigen::Matrix3d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f>(alg->Identity(typeid(Eigen::Matrix3f))).array()==Eigen::Matrix3f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i>(alg->Identity(typeid(Eigen::Matrix3i))).array()==Eigen::Matrix3i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d>(alg->Identity(typeid(Eigen::Matrix4d))).array()==Eigen::Matrix4d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f>(alg->Identity(typeid(Eigen::Matrix4f))).array()==Eigen::Matrix4f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i>(alg->Identity(typeid(Eigen::Matrix4i))).array()==Eigen::Matrix4i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd>(alg->Identity(typeid(Eigen::MatrixXd), 4, 9)).array()==Eigen::MatrixXd::Identity(4, 9).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf>(alg->Identity(typeid(Eigen::MatrixXf), 43, 21)).array()==Eigen::MatrixXf::Identity(43, 21).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi>(alg->Identity(typeid(Eigen::MatrixXi), 23, 23)).array()==Eigen::MatrixXi::Identity(23, 23).array()).all());
}

TEST(EigenMatrixAlgebraTests, Norm) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test Eigen::Matrix norm
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd matXd = Eigen::MatrixXd::Random(8, 4);
  EXPECT_DOUBLE_EQ(alg->Norm(mat2d), mat2d.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(mat3d), mat3d.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(mat4d), mat4d.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(matXd), matXd.norm());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf matXf = Eigen::MatrixXf::Random(8, 4);
  EXPECT_FLOAT_EQ(alg->Norm(mat2f), mat2f.norm());
  EXPECT_FLOAT_EQ(alg->Norm(mat3f), mat3f.norm());
  EXPECT_FLOAT_EQ(alg->Norm(mat4f), mat4f.norm());
  EXPECT_FLOAT_EQ(alg->Norm(matXf), matXf.norm());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi matXi = Eigen::MatrixXi::Random(8, 4);
  EXPECT_EQ(alg->Norm(mat2i), mat2i.norm());
  EXPECT_EQ(alg->Norm(mat3i), mat3i.norm());
  EXPECT_EQ(alg->Norm(mat4i), mat4i.norm());
  EXPECT_EQ(alg->Norm(matXi), matXi.norm());
}

TEST(EigenMatrixAlgebraTests, Add) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Matrices
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(mat2d, mat2d)).array()==(mat2d+mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(mat2d, mat2Xd)).array()==(mat2d+mat2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(mat2Xd, mat2d)).array()==(mat2d+mat2Xd).array()).all());
  const Eigen::Vector2d vec2d(1.0, 1.0);
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  Eigen::Matrix2d result2d = mat2d;
  result2d.diagonal() += vec2d;
  Eigen::MatrixXd result2Xd = mat2d;
  result2Xd.diagonal() += vec2Xd;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(mat2d, vec2d)).array()==result2d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(mat2d, vec2Xd)).array()==result2Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(vec2d, mat2d)).array()==result2d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(vec2Xd, mat2d)).array()==result2Xd.array()).all());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(mat2f, mat2f)).array()==(mat2f+mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(mat2f, mat2Xf)).array()==(mat2f+mat2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(mat2Xf, mat2f)).array()==(mat2f+mat2Xf).array()).all());
  const Eigen::Vector2f vec2f(1.0, 1.0);
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  Eigen::Matrix2f result2f = mat2f;
  result2f.diagonal() += vec2f;
  Eigen::MatrixXf result2Xf = mat2f;
  result2Xf.diagonal() += vec2Xf;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(mat2f, vec2f)).array()==result2f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(mat2f, vec2Xf)).array()==result2Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(vec2f, mat2f)).array()==result2f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(vec2Xf, mat2f)).array()==result2Xf.array()).all());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(mat2i, mat2i)).array()==(mat2i+mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(mat2i, mat2Xi)).array()==(mat2i+mat2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mat2Xi, mat2i)).array()==(mat2i+mat2Xi).array()).all());
  const Eigen::Vector2i vec2i(1.0, 1.0);
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  Eigen::Matrix2i result2i = mat2i;
  result2i.diagonal() += vec2i;
  Eigen::MatrixXi result2Xi = mat2i;
  result2Xi.diagonal() += vec2Xi;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(mat2i, vec2i)).array()==result2i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(mat2i, vec2Xi)).array()==result2Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(vec2i, mat2i)).array()==result2i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(vec2Xi, mat2i)).array()==result2Xi.array()).all());

  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(mat3d, mat3d)).array()==(mat3d+mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(mat3d, mat3Xd)).array()==(mat3d+mat3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(mat3Xd, mat3d)).array()==(mat3d+mat3Xd).array()).all());
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  Eigen::Matrix3d result3d = mat3d;
  result3d.diagonal() += vec3d;
  Eigen::MatrixXd result3Xd = mat3d;
  result3Xd.diagonal() += vec3Xd;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(mat3d, vec3d)).array()==result3d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(mat3d, vec3Xd)).array()==result3Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(vec3d, mat3d)).array()==result3d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(vec3Xd, mat3d)).array()==result3Xd.array()).all());

  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(mat3f, mat3f)).array()==(mat3f+mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(mat3f, mat3Xf)).array()==(mat3f+mat3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(mat3Xf, mat3f)).array()==(mat3f+mat3Xf).array()).all());
  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  Eigen::Matrix3f result3f = mat3f;
  result3f.diagonal() += vec3f;
  Eigen::MatrixXf result3Xf = mat3f;
  result3Xf.diagonal() += vec3Xf;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(mat3f, vec3f)).array()==result3f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(mat3f, vec3Xf)).array()==result3Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(vec3f, mat3f)).array()==result3f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(vec3Xf, mat3f)).array()==result3Xf.array()).all());

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(mat3i, mat3i)).array()==(mat3i+mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(mat3i, mat3Xi)).array()==(mat3i+mat3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mat3Xi, mat3i)).array()==(mat3i+mat3Xi).array()).all());
  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  Eigen::Matrix3i result3i = mat3i;
  result3i.diagonal() += vec3i;
  Eigen::MatrixXi result3Xi = mat3i;
  result3Xi.diagonal() += vec3Xi;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(mat3i, vec3i)).array()==result3i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(mat3i, vec3Xi)).array()==result3Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(vec3i, mat3i)).array()==result3i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(vec3Xi, mat3i)).array()==result3Xi.array()).all());

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(mat4d, mat4d)).array()==(mat4d+mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(mat4d, mat4Xd)).array()==(mat4d+mat4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(mat4Xd, mat4d)).array()==(mat4d+mat4Xd).array()).all());
  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  Eigen::Matrix4d result4d = mat4d;
  result4d.diagonal() += vec4d;
  Eigen::MatrixXd result4Xd = mat4d;
  result4Xd.diagonal() += vec4Xd;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(mat4d, vec4d)).array()==result4d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(mat4d, vec4Xd)).array()==result4Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(vec4d, mat4d)).array()==result4d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(vec4Xd, mat4d)).array()==result4Xd.array()).all());

  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(mat4f, mat4f)).array()==(mat4f+mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(mat4f, mat4Xf)).array()==(mat4f+mat4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(mat4Xf, mat4f)).array()==(mat4f+mat4Xf).array()).all());
  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  Eigen::Matrix4f result4f = mat4f;
  result4f.diagonal() += vec4f;
  Eigen::MatrixXf result4Xf = mat4f;
  result4Xf.diagonal() += vec4Xf;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(mat4f, vec4f)).array()==result4f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(mat4f, vec4Xf)).array()==result4Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(vec4f, mat4f)).array()==result4f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(vec4Xf, mat4f)).array()==result4Xf.array()).all());

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(mat4i, mat4i)).array()==(mat4i+mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(mat4i, mat4Xi)).array()==(mat4i+mat4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mat4Xi, mat4i)).array()==(mat4i+mat4Xi).array()).all());
  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  Eigen::Matrix4i result4i = mat4i;
  result4i.diagonal() += vec4i;
  Eigen::MatrixXi result4Xi = mat4i;
  result4Xi.diagonal() += vec4Xi;
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(mat4i, vec4i)).array()==result4i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(mat4i, vec4Xi)).array()==result4Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(vec4i, mat4i)).array()==result4i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(vec4Xi, mat4i)).array()==result4Xi.array()).all());

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(8,8);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(9,9);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5,5);
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(matd, matd)).array()==(matd+matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(matf, matf)).array()==(matf+matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mati, mati)).array()==(mati+mati).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(8);
  Eigen::MatrixXd resultd = matd;
  resultd.diagonal() += vecd;
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(matd, vecd)).array()==resultd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(vecd, matd)).array()==resultd.array()).all());

  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  Eigen::MatrixXf resultf = matf;
  resultf.diagonal() += vecf;
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(matf, vecf)).array()==resultf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(vecf, matf)).array()==resultf.array()).all());

  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  Eigen::MatrixXi resulti = mati;
  resulti.diagonal() += veci;
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mati, veci)).array()==resulti.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(veci, mati)).array()==resulti.array()).all());
}

TEST(EigenMatrixAlgebraTests, Subtract) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Matrices
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Subtract(mat2d, mat2d)).array()==(mat2d-mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Subtract(mat2d, mat2Xd)).array()==(mat2d-mat2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(mat2Xd, mat2d)).array()==(mat2Xd-mat2d).array()).all());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Subtract(mat2f, mat2f)).array()==(mat2f-mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Subtract(mat2f, mat2Xf)).array()==(mat2f-mat2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(mat2Xf, mat2f)).array()==(mat2Xf-mat2f).array()).all());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Subtract(mat2i, mat2i)).array()==(mat2i-mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Subtract(mat2i, mat2Xi)).array()==(mat2i-mat2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mat2Xi, mat2i)).array()==(mat2Xi-mat2i).array()).all());

  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Subtract(mat3d, mat3d)).array()==(mat3d-mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Subtract(mat3d, mat3Xd)).array()==(mat3d-mat3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(mat3Xd, mat3d)).array()==(mat3Xd-mat3d).array()).all());

  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Subtract(mat3f, mat3f)).array()==(mat3f-mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Subtract(mat3f, mat3Xf)).array()==(mat3f-mat3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(mat3Xf, mat3f)).array()==(mat3Xf-mat3f).array()).all());

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Subtract(mat3i, mat3i)).array()==(mat3i-mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Subtract(mat3i, mat3Xi)).array()==(mat3i-mat3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mat3Xi, mat3i)).array()==(mat3Xi-mat3i).array()).all());

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Subtract(mat4d, mat4d)).array()==(mat4d-mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Subtract(mat4d, mat4Xd)).array()==(mat4d-mat4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(mat4Xd, mat4d)).array()==(mat4Xd-mat4d).array()).all());

  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Subtract(mat4f, mat4f)).array()==(mat4f-mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Subtract(mat4f, mat4Xf)).array()==(mat4f-mat4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(mat4Xf, mat4f)).array()==(mat4Xf-mat4f).array()).all());

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Subtract(mat4i, mat4i)).array()==(mat4i-mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Subtract(mat4i, mat4Xi)).array()==(mat4i-mat4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mat4Xi, mat4i)).array()==(mat4Xi-mat4i).array()).all());

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(8,62);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(9,2);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5,5);
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(matd, matd)).array()==(matd-matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(matf, matf)).array()==(matf-matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mati, mati)).array()==(mati-mati).array()).all());
}

TEST(EigenMatrixAlgebraTests, Multiply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  // test the Eigen::Matrices
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(xd, mat2d)).array()==(xd*mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(mat2d, xd)).array()==(xd*mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(mat2d, mat2d)).array()==(mat2d*mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(mat2d, mat2Xd)).array()==(mat2d*mat2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(mat2Xd, mat2d)).array()==(mat2Xd*mat2d).array()).all());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(xf, mat2f)).array()==(xf*mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(mat2f, xf)).array()==(xf*mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(mat2f, mat2f)).array()==(mat2f*mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(mat2f, mat2Xf)).array()==(mat2f*mat2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(mat2Xf, mat2f)).array()==(mat2Xf*mat2f).array()).all());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(xi, mat2i)).array()==(xi*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, xi)).array()==(xi*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(xui, mat2i)).array()==(xui*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, xui)).array()==(xui*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, mat2i)).array()==(mat2i*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, mat2Xi)).array()==(mat2i*mat2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mat2Xi, mat2i)).array()==(mat2Xi*mat2i).array()).all());

  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(xd, mat3d)).array()==(xd*mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(mat3d, xd)).array()==(xd*mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(mat3d, mat3d)).array()==(mat3d*mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(mat3d, mat3Xd)).array()==(mat3d*mat3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(mat3Xd, mat3d)).array()==(mat3Xd*mat3d).array()).all());

  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(xf, mat3f)).array()==(xf*mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(mat3f, xf)).array()==(xf*mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(mat3f, mat3f)).array()==(mat3f*mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(mat3f, mat3Xf)).array()==(mat3f*mat3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(mat3Xf, mat3f)).array()==(mat3Xf*mat3f).array()).all());

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(xi, mat3i)).array()==(xi*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, xi)).array()==(xi*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(xui, mat3i)).array()==(xui*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, xui)).array()==(xui*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, mat3i)).array()==(mat3i*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, mat3Xi)).array()==(mat3i*mat3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mat3Xi, mat3i)).array()==(mat3Xi*mat3i).array()).all());

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(xd, mat4d)).array()==(xd*mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(mat4d, xd)).array()==(xd*mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(mat4d, mat4d)).array()==(mat4d*mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(mat4d, mat4Xd)).array()==(mat4d*mat4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(mat4Xd, mat4d)).array()==(mat4Xd*mat4d).array()).all());

  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(xf, mat4f)).array()==(xf*mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(mat4f, xf)).array()==(xf*mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(mat4f, mat4f)).array()==(mat4f*mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(mat4f, mat4Xf)).array()==(mat4f*mat4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(mat4Xf, mat4f)).array()==(mat4Xf*mat4f).array()).all());

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(xi, mat4i)).array()==(xi*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, xi)).array()==(xi*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(xui, mat4i)).array()==(xui*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, xui)).array()==(xui*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, mat4i)).array()==(mat4i*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, mat4Xi)).array()==(mat4i*mat4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mat4Xi, mat4i)).array()==(mat4Xi*mat4i).array()).all());

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(8,8);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(9,9);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5,5);
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(xd, matd)).array()==(xd*matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(matd, xd)).array()==(xd*matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(matd, matd)).array()==(matd*matd).array()).all());

  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(xf, matf)).array()==(xf*matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(matf, xf)).array()==(xf*matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(matf, matf)).array()==(matf*matf).array()).all());

  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(xi, mati)).array()==(xi*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mati, xi)).array()==(xi*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(xui, mati)).array()==(xui*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mati, xui)).array()==(xui*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mati, mati)).array()==(mati*mati).array()).all());
}

TEST(EigenMatrixAlgebraTests, Apply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->Apply(mat2d, vec2d)).array()==(mat2d*vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(mat2d, vec2Xd)).array()==(mat2d*vec2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->Apply(mat2Xd, vec2d)).array()==(mat2Xd*vec2d).array()).all());

  mat2d = Eigen::Matrix2d::Identity() + 1.0e-5*mat2d*mat2d.transpose();
  Eigen::LLT<Eigen::Matrix2d> chol2d;
  chol2d.compute(mat2d);
  mat2Xd = Eigen::MatrixXd::Identity(2,2) + 1.0e-5*mat2Xd*mat2Xd.transpose();
  Eigen::LLT<Eigen::Matrix2d> chol2Xd;
  chol2Xd.compute(mat2Xd);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector2d>(alg->Apply(chol2d, vec2d))-(mat2d*vec2d)).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXd>(alg->Apply(chol2d, vec2Xd))-(mat2d*vec2Xd)).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector2d>(alg->Apply(chol2Xd, vec2d))-(mat2Xd*vec2d)).norm(), 0.0, 1.0e-10);

  Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->Apply(mat2f, vec2f)).array()==(mat2f*vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(mat2f, vec2Xf)).array()==(mat2f*vec2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->Apply(mat2Xf, vec2f)).array()==(mat2Xf*vec2f).array()).all());

  mat2f = Eigen::Matrix2f::Identity() + 1.0e-5*mat2f*mat2f.transpose();
  Eigen::LLT<Eigen::Matrix2f> chol2f;
  chol2f.compute(mat2f);
  mat2Xf = Eigen::MatrixXf::Identity(2,2) + 1.0e-5*mat2Xf*mat2Xf.transpose();
  Eigen::LLT<Eigen::Matrix2f> chol2Xf;
  chol2Xf.compute(mat2Xf);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector2f>(alg->Apply(chol2f, vec2f))-(mat2f*vec2f)).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXf>(alg->Apply(chol2f, vec2Xf))-(mat2f*vec2Xf)).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector2f>(alg->Apply(chol2Xf, vec2f))-(mat2Xf*vec2f)).norm(), 0.0, 1.0e-6);

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->Apply(mat2i, vec2i)).array()==(mat2i*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(mat2i, vec2Xi)).array()==(mat2i*vec2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->Apply(mat2Xi, vec2i)).array()==(mat2Xi*vec2i).array()).all());

  Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->Apply(mat3d, vec3d)).array()==(mat3d*vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(mat3d, vec3Xd)).array()==(mat3d*vec3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->Apply(mat3Xd, vec3d)).array()==(mat3Xd*vec3d).array()).all());

  mat3d = Eigen::Matrix3d::Identity() + 1.0e-5*mat3d*mat3d.transpose();
  Eigen::LLT<Eigen::Matrix3d> chol3d;
  chol3d.compute(mat3d);
  mat3Xd = Eigen::MatrixXd::Identity(3,3) + 1.0e-5*mat3Xd*mat3Xd.transpose();
  Eigen::LLT<Eigen::Matrix3d> chol3Xd;
  chol3Xd.compute(mat3Xd);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector3d>(alg->Apply(chol3d, vec3d))-(mat3d*vec3d)).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXd>(alg->Apply(chol3d, vec3Xd))-(mat3d*vec3Xd)).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector3d>(alg->Apply(chol3Xd, vec3d))-(mat3Xd*vec3d)).norm(), 0.0, 1.0e-10);

  Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->Apply(mat3f, vec3f)).array()==(mat3f*vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(mat3f, vec3Xf)).array()==(mat3f*vec3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->Apply(mat3Xf, vec3f)).array()==(mat3Xf*vec3f).array()).all());

  mat3f = Eigen::Matrix3f::Identity() + 1.0e-5*mat3f*mat3f.transpose();
  Eigen::LLT<Eigen::Matrix3f> chol3f;
  chol3f.compute(mat3f);
  mat3Xf = Eigen::MatrixXf::Identity(3,3) + 1.0e-5*mat3Xf*mat3Xf.transpose();
  Eigen::LLT<Eigen::Matrix3f> chol3Xf;
  chol3Xf.compute(mat3Xf);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector3f>(alg->Apply(chol3f, vec3f))-(mat3f*vec3f)).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXf>(alg->Apply(chol3f, vec3Xf))-(mat3f*vec3Xf)).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector3f>(alg->Apply(chol3Xf, vec3f))-(mat3Xf*vec3f)).norm(), 0.0, 1.0e-6);

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->Apply(mat3i, vec3i)).array()==(mat3i*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(mat3i, vec3Xi)).array()==(mat3i*vec3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->Apply(mat3Xi, vec3i)).array()==(mat3Xi*vec3i).array()).all());

  Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->Apply(mat4d, vec4d)).array()==(mat4d*vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(mat4d, vec4Xd)).array()==(mat4d*vec4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->Apply(mat4Xd, vec4d)).array()==(mat4Xd*vec4d).array()).all());

  mat4d = Eigen::Matrix4d::Identity() + 1.0e-5*mat4d*mat4d.transpose();
  Eigen::LLT<Eigen::Matrix4d> chol4d;
  chol4d.compute(mat4d);
  mat4Xd = Eigen::MatrixXd::Identity(4,4) + 1.0e-5*mat4Xd*mat4Xd.transpose();
  Eigen::LLT<Eigen::Matrix4d> chol4Xd;
  chol4Xd.compute(mat4Xd);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector4d>(alg->Apply(chol4d, vec4d))-(mat4d*vec4d)).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXd>(alg->Apply(chol4d, vec4Xd))-(mat4d*vec4Xd)).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector4d>(alg->Apply(chol4Xd, vec4d))-(mat4Xd*vec4d)).norm(), 0.0, 1.0e-10);

  Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->Apply(mat4f, vec4f)).array()==(mat4f*vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(mat4f, vec4Xf)).array()==(mat4f*vec4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->Apply(mat4Xf, vec4f)).array()==(mat4Xf*vec4f).array()).all());

  mat4f = Eigen::Matrix4f::Identity() + 1.0e-5*mat4f*mat4f.transpose();
  Eigen::LLT<Eigen::Matrix4f> chol4f;
  chol4f.compute(mat4f);
  mat4Xf = Eigen::MatrixXf::Identity(4,4) + 1.0e-5*mat4Xf*mat4Xf.transpose();
  Eigen::LLT<Eigen::Matrix4f> chol4Xf;
  chol4Xf.compute(mat4Xf);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector4f>(alg->Apply(chol4f, vec4f))-(mat4f*vec4f)).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXf>(alg->Apply(chol4f, vec4Xf))-(mat4f*vec4Xf)).norm(), 0.0, 1.0e-6);
  EXPECT_NEAR((boost::any_cast<Eigen::Vector4f>(alg->Apply(chol4Xf, vec4f))-(mat4Xf*vec4f)).norm(), 0.0, 1.0e-6);

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->Apply(mat4i, vec4i)).array()==(mat4i*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(mat4i, vec4Xi)).array()==(mat4i*vec4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->Apply(mat4Xi, vec4i)).array()==(mat4Xi*vec4i).array()).all());

  Eigen::MatrixXd matXd = Eigen::MatrixXd::Random(5,8);
  const Eigen::VectorXd vecXd = Eigen::VectorXd::Random(8);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(matXd, vecXd)).array()==(matXd*vecXd).array()).all());
  matXd = Eigen::MatrixXd::Random(8,8);
  matXd = Eigen::MatrixXd::Identity(8,8) + 1.0e-5*matXd*matXd.transpose();
  Eigen::LLT<Eigen::MatrixXd> cholXd;
  cholXd.compute(matXd);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXd>(alg->Apply(cholXd, vecXd))-(matXd*vecXd)).norm(), 0.0, 1.0e-10);

  Eigen::MatrixXf matXf = Eigen::MatrixXf::Random(13,13);
  const Eigen::VectorXf vecXf = Eigen::VectorXf::Random(13);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(matXf, vecXf)).array()==(matXf*vecXf).array()).all());
  matXf = Eigen::MatrixXf::Identity(13,13) + 1.0e-5*matXf*matXf.transpose();
  Eigen::LLT<Eigen::MatrixXf> cholXf;
  cholXf.compute(matXf);
  EXPECT_NEAR((boost::any_cast<Eigen::VectorXf>(alg->Apply(cholXf, vecXf))-(matXf*vecXf)).norm(), 0.0, 1.0e-6);

  const Eigen::MatrixXi matXi = Eigen::MatrixXi::Random(4,2);
  const Eigen::VectorXi vecXi = Eigen::VectorXi::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(matXi, vecXi)).array()==(matXi*vecXi).array()).all());
}

TEST(EigenMatrixAlgebraTests, ApplyInverse) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  Eigen::VectorXd soln2d = boost::any_cast<Eigen::Vector2d const>(alg->ApplyInverse(mat2d, vec2d));
  Eigen::VectorXd soln2Xd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(mat2d, vec2Xd));
  Eigen::Vector2d soln2d_mat = boost::any_cast<Eigen::Vector2d const>(alg->ApplyInverse(mat2Xd, vec2d));
  EXPECT_NEAR((mat2d*soln2d-vec2d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat2Xd*soln2d_mat-vec2d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat2d*soln2Xd-vec2Xd).norm(), 0.0, 1.0e-10);

  mat2d = Eigen::Matrix2d::Identity() + 1.0e-2*mat2d*mat2d.transpose();
  Eigen::LLT<Eigen::Matrix2d> chol2d;
  chol2d.compute(mat2d);
  soln2d = boost::any_cast<Eigen::Vector2d const>(alg->ApplyInverse(chol2d, vec2d));
  soln2Xd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(chol2d, vec2Xd));
  EXPECT_NEAR((mat2d*soln2d-vec2d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat2d*soln2Xd-vec2Xd).norm(), 0.0, 1.0e-10);

  Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  Eigen::Vector2f soln2f = boost::any_cast<Eigen::Vector2f const>(alg->ApplyInverse(mat2f, vec2f));
  Eigen::VectorXf soln2Xf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(mat2f, vec2Xf));
  Eigen::Vector2f soln2f_mat = boost::any_cast<Eigen::Vector2f const>(alg->ApplyInverse(mat2Xf, vec2f));
  EXPECT_NEAR((mat2f*soln2f-vec2f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat2Xf*soln2f_mat-vec2f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat2f*soln2Xf-vec2Xf).norm(), 0.0, 1.0e-5);

  mat2f = Eigen::Matrix2f::Identity() + 1.0e-2*mat2f*mat2f.transpose();
  Eigen::LLT<Eigen::Matrix2f> chol2f;
  chol2f.compute(mat2f);
  soln2f = boost::any_cast<Eigen::Vector2f const>(alg->ApplyInverse(chol2f, vec2f));
  soln2Xf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(chol2f, vec2Xf));
  EXPECT_NEAR((mat2f*soln2f-vec2f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat2f*soln2Xf-vec2Xf).norm(), 0.0, 1.0e-5);

  Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  Eigen::VectorXd soln3d = boost::any_cast<Eigen::Vector3d const>(alg->ApplyInverse(mat3d, vec3d));
  Eigen::VectorXd soln3Xd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(mat3d, vec3Xd));
  Eigen::Vector3d soln3d_mat = boost::any_cast<Eigen::Vector3d const>(alg->ApplyInverse(mat3Xd, vec3d));
  EXPECT_NEAR((mat3d*soln3d-vec3d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat3Xd*soln3d_mat-vec3d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat3d*soln3Xd-vec3Xd).norm(), 0.0, 1.0e-10);

  mat3d = Eigen::Matrix3d::Identity() + 1.0e-3*mat3d*mat3d.transpose();
  Eigen::LLT<Eigen::Matrix3d> chol3d;
  chol3d.compute(mat3d);
  soln3d = boost::any_cast<Eigen::Vector3d const>(alg->ApplyInverse(chol3d, vec3d));
  soln3Xd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(chol3d, vec3Xd));
  EXPECT_NEAR((mat3d*soln3d-vec3d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat3d*soln3Xd-vec3Xd).norm(), 0.0, 1.0e-10);

  Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  Eigen::Vector3f soln3f = boost::any_cast<Eigen::Vector3f const>(alg->ApplyInverse(mat3f, vec3f));
  Eigen::VectorXf soln3Xf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(mat3f, vec3Xf));
  Eigen::Vector3f soln3f_mat = boost::any_cast<Eigen::Vector3f const>(alg->ApplyInverse(mat3Xf, vec3f));
  EXPECT_NEAR((mat3f*soln3f-vec3f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat3Xf*soln3f_mat-vec3f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat3f*soln3Xf-vec3Xf).norm(), 0.0, 1.0e-5);

  mat3f = Eigen::Matrix3f::Identity() + 1.0e-3*mat3f*mat3f.transpose();
  Eigen::LLT<Eigen::Matrix3f> chol3f;
  chol3f.compute(mat3f);
  soln3f = boost::any_cast<Eigen::Vector3f const>(alg->ApplyInverse(chol3f, vec3f));
  soln3Xf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(chol3f, vec3Xf));
  EXPECT_NEAR((mat3f*soln3f-vec3f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat3f*soln3Xf-vec3Xf).norm(), 0.0, 1.0e-5);

  Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  Eigen::VectorXd soln4d = boost::any_cast<Eigen::Vector4d const>(alg->ApplyInverse(mat4d, vec4d));
  Eigen::VectorXd soln4Xd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(mat4d, vec4Xd));
  Eigen::Vector4d soln4d_mat = boost::any_cast<Eigen::Vector4d const>(alg->ApplyInverse(mat4Xd, vec4d));
  EXPECT_NEAR((mat4d*soln4d-vec4d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat4Xd*soln4d_mat-vec4d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat4d*soln4Xd-vec4Xd).norm(), 0.0, 1.0e-10);

  mat4d = Eigen::Matrix4d::Identity() + 1.0e-4*mat4d*mat4d.transpose();
  Eigen::LLT<Eigen::Matrix4d> chol4d;
  chol4d.compute(mat4d);
  soln4d = boost::any_cast<Eigen::Vector4d const>(alg->ApplyInverse(chol4d, vec4d));
  soln4Xd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(chol4d, vec4Xd));
  EXPECT_NEAR((mat4d*soln4d-vec4d).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((mat4d*soln4Xd-vec4Xd).norm(), 0.0, 1.0e-10);

  Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  Eigen::Vector4f soln4f = boost::any_cast<Eigen::Vector4f const>(alg->ApplyInverse(mat4f, vec4f));
  Eigen::VectorXf soln4Xf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(mat4f, vec4Xf));
  Eigen::Vector4f soln4f_mat = boost::any_cast<Eigen::Vector4f const>(alg->ApplyInverse(mat4Xf, vec4f));
  EXPECT_NEAR((mat4f*soln4f-vec4f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat4Xf*soln4f_mat-vec4f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat4f*soln4Xf-vec4Xf).norm(), 0.0, 1.0e-5);

  mat4f = Eigen::Matrix4f::Identity() + 1.0e-4*mat4f*mat4f.transpose();
  Eigen::LLT<Eigen::Matrix4f> chol4f;
  chol4f.compute(mat4f);
  soln4f = boost::any_cast<Eigen::Vector4f const>(alg->ApplyInverse(chol4f, vec4f));
  soln4Xf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(chol4f, vec4Xf));
  EXPECT_NEAR((mat4f*soln4f-vec4f).norm(), 0.0, 1.0e-5);
  EXPECT_NEAR((mat4f*soln4Xf-vec4Xf).norm(), 0.0, 1.0e-5);

  Eigen::MatrixXd matXd = Eigen::MatrixXd::Random(5, 5);;
  const Eigen::VectorXd vecXd = Eigen::VectorXd::Random(5);
  Eigen::VectorXd solnXd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(matXd, vecXd));
  EXPECT_NEAR((matXd*solnXd-vecXd).norm(), 0.0, 1.0e-10);

  matXd = Eigen::MatrixXd::Identity(5,5) + 1.0e-4*matXd*matXd.transpose();
  Eigen::LLT<Eigen::MatrixXd> cholXd;
  cholXd.compute(matXd);
  solnXd = boost::any_cast<Eigen::VectorXd const>(alg->ApplyInverse(cholXd, vecXd));
  EXPECT_NEAR((matXd*solnXd-vecXd).norm(), 0.0, 1.0e-10);

  Eigen::MatrixXf matXf = Eigen::MatrixXf::Random(5, 5);;
  const Eigen::VectorXf vecXf = Eigen::VectorXf::Random(5);
  Eigen::VectorXf solnXf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(matXf, vecXf));
  EXPECT_NEAR((matXf*solnXf-vecXf).norm(), 0.0, 1.0e-5);

  matXf = Eigen::MatrixXf::Identity(5,5) + 1.0e-4*matXf*matXf.transpose();
  Eigen::LLT<Eigen::MatrixXf> cholXf;
  cholXf.compute(matXf);
  solnXf = boost::any_cast<Eigen::VectorXf const>(alg->ApplyInverse(cholXf, vecXf));
  EXPECT_NEAR((matXf*solnXf-vecXf).norm(), 0.0, 1.0e-5);
}

TEST(EigenMatrixAlgebraTests, SquareRoot) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  mat2d = Eigen::Matrix2d::Identity() + 1.0e-4*mat2d*mat2d.transpose();
  Eigen::LLT<Eigen::Matrix2d> chol2d;
  chol2d.compute(mat2d);
  const Eigen::Matrix2d L2d = boost::any_cast<Eigen::Matrix2d>(alg->SquareRoot(chol2d));
  EXPECT_NEAR((L2d*L2d.transpose()-mat2d).norm(), 0.0, 1.0e-10);

  Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  mat2f = Eigen::Matrix2f::Identity() + 1.0e-4*mat2f*mat2f.transpose();
  Eigen::LLT<Eigen::Matrix2f> chol2f;
  chol2f.compute(mat2f);
  const Eigen::Matrix2f L2f = boost::any_cast<Eigen::Matrix2f>(alg->SquareRoot(chol2f));
  EXPECT_NEAR((L2f*L2f.transpose()-mat2f).norm(), 0.0, 1.0e-6);

  Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  mat3d = Eigen::Matrix3d::Identity() + 1.0e-4*mat3d*mat3d.transpose();
  Eigen::LLT<Eigen::Matrix3d> chol3d;
  chol3d.compute(mat3d);
  const Eigen::Matrix3d L3d = boost::any_cast<Eigen::Matrix3d>(alg->SquareRoot(chol3d));
  EXPECT_NEAR((L3d*L3d.transpose()-mat3d).norm(), 0.0, 1.0e-10);

  Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  mat3f = Eigen::Matrix3f::Identity() + 1.0e-4*mat3f*mat3f.transpose();
  Eigen::LLT<Eigen::Matrix3f> chol3f;
  chol3f.compute(mat3f);
  const Eigen::Matrix3f L3f = boost::any_cast<Eigen::Matrix3f>(alg->SquareRoot(chol3f));
  EXPECT_NEAR((L3f*L3f.transpose()-mat3f).norm(), 0.0, 1.0e-6);

  Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  mat4d = Eigen::Matrix4d::Identity() + 1.0e-4*mat4d*mat4d.transpose();
  Eigen::LLT<Eigen::Matrix4d> chol4d;
  chol4d.compute(mat4d);
  const Eigen::Matrix4d L4d = boost::any_cast<Eigen::Matrix4d>(alg->SquareRoot(chol4d));
  EXPECT_NEAR((L4d*L4d.transpose()-mat4d).norm(), 0.0, 1.0e-10);

  Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  mat4f = Eigen::Matrix4f::Identity() + 1.0e-4*mat4f*mat4f.transpose();
  Eigen::LLT<Eigen::Matrix4f> chol4f;
  chol4f.compute(mat4f);
  const Eigen::Matrix4f L4f = boost::any_cast<Eigen::Matrix4f>(alg->SquareRoot(chol4f));
  EXPECT_NEAR((L4f*L4f.transpose()-mat4f).norm(), 0.0, 1.0e-6);

  Eigen::MatrixXd matXd = Eigen::MatrixXd::Random(8,8);
  matXd = Eigen::MatrixXd::Identity(8,8) + 1.0e-4*matXd*matXd.transpose();
  Eigen::LLT<Eigen::MatrixXd> cholXd;
  cholXd.compute(matXd);
  const Eigen::MatrixXd LXd = boost::any_cast<Eigen::MatrixXd>(alg->SquareRoot(cholXd));
  EXPECT_NEAR((LXd*LXd.transpose()-matXd).norm(), 0.0, 1.0e-10);

  Eigen::MatrixXf matXf = Eigen::MatrixXf::Random(13,13);
  matXf = Eigen::MatrixXf::Identity(13,13) + 1.0e-4*matXf*matXf.transpose();
  Eigen::LLT<Eigen::MatrixXf> cholXf;
  cholXf.compute(matXf);
  const Eigen::MatrixXf LXf = boost::any_cast<Eigen::MatrixXf>(alg->SquareRoot(cholXf));
  EXPECT_NEAR((LXf*LXf.transpose()-matXf).norm(), 0.0, 1.0e-6);

}

TEST(EigenMatrixAlgebraTests, Determinate) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  mat2d = Eigen::Matrix2d::Identity() + 1.0e-4*mat2d*mat2d.transpose();
  Eigen::LLT<Eigen::Matrix2d> chol2d;
  chol2d.compute(mat2d);
  Eigen::Matrix2d L2d = chol2d.matrixL();
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(chol2d), 2.0*L2d.diagonal().array().log().sum());

  Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  mat2f = Eigen::Matrix2f::Identity() + 1.0e-4*mat2f*mat2f.transpose();
  Eigen::LLT<Eigen::Matrix2f> chol2f;
  chol2f.compute(mat2f);
  Eigen::Matrix2f L2f = chol2f.matrixL();
  EXPECT_FLOAT_EQ(alg->LogDeterminate(chol2f), 2.0*L2f.diagonal().array().log().sum());

  Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  mat3d = Eigen::Matrix3d::Identity() + 1.0e-4*mat3d*mat3d.transpose();
  Eigen::LLT<Eigen::Matrix3d> chol3d;
  chol3d.compute(mat3d);
  Eigen::Matrix3d L3d = chol3d.matrixL();
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(chol3d), 2.0*L3d.diagonal().array().log().sum());

  Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  mat3f = Eigen::Matrix3f::Identity() + 1.0e-4*mat3f*mat3f.transpose();
  Eigen::LLT<Eigen::Matrix3f> chol3f;
  chol3f.compute(mat3f);
  Eigen::Matrix3f L3f = chol3f.matrixL();
  EXPECT_FLOAT_EQ(alg->LogDeterminate(chol3f), 2.0*L3f.diagonal().array().log().sum());

  Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  mat4d = Eigen::Matrix4d::Identity() + 1.0e-4*mat4d*mat4d.transpose();
  Eigen::LLT<Eigen::Matrix4d> chol4d;
  chol4d.compute(mat4d);
  Eigen::Matrix4d L4d = chol4d.matrixL();
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(chol4d), 2.0*L4d.diagonal().array().log().sum());

  Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  mat4f = Eigen::Matrix4f::Identity() + 1.0e-4*mat4f*mat4f.transpose();
  Eigen::LLT<Eigen::Matrix4f> chol4f;
  chol4f.compute(mat4f);
  Eigen::Matrix4f L4f = chol4f.matrixL();
  EXPECT_FLOAT_EQ(alg->LogDeterminate(chol4f), 2.0*L4f.diagonal().array().log().sum());

  Eigen::MatrixXd matXd = Eigen::MatrixXd::Random(5,5);
  matXd = Eigen::MatrixXd::Identity(5,5) + 1.0e-4*matXd*matXd.transpose();
  Eigen::LLT<Eigen::MatrixXd> cholXd;
  cholXd.compute(matXd);
  Eigen::MatrixXd LXd = cholXd.matrixL();
  EXPECT_DOUBLE_EQ(alg->LogDeterminate(cholXd), 2.0*LXd.diagonal().array().log().sum());

  Eigen::MatrixXf matXf = Eigen::MatrixXf::Random(6,6);
  matXf = Eigen::MatrixXf::Identity(6,6) + 1.0e-4*matXf*matXf.transpose();
  Eigen::LLT<Eigen::MatrixXf> cholXf;
  cholXf.compute(matXf);
  Eigen::MatrixXf LXf = cholXf.matrixL();
  EXPECT_FLOAT_EQ(alg->LogDeterminate(cholXf), 2.0*LXf.diagonal().array().log().sum());
}
