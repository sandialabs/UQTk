#include <gtest/gtest.h>


#include <memory>

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Modeling;

TEST(EigenVectorAlgebraTests, Size) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // check for the Eigen::Vector's
  const Eigen::Vector2d test2d = Eigen::Vector2d::Random();
  const Eigen::Vector3d test3d = Eigen::Vector3d::Random();
  const Eigen::Vector4d test4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd testXd = Eigen::VectorXd::Random(6);
  EXPECT_EQ(alg->Size(test2d), 2);
  EXPECT_EQ(alg->Size(test3d), 3);
  EXPECT_EQ(alg->Size(test4d), 4);
  EXPECT_EQ(alg->Size(testXd), 6);

  const Eigen::Vector2f test2f = Eigen::Vector2f::Random();
  const Eigen::Vector3f test3f = Eigen::Vector3f::Random();
  const Eigen::Vector4f test4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf testXf = Eigen::VectorXf::Random(6);
  EXPECT_EQ(alg->Size(test2f), 2);
  EXPECT_EQ(alg->Size(test3f), 3);
  EXPECT_EQ(alg->Size(test4f), 4);
  EXPECT_EQ(alg->Size(testXf), 6);

  const Eigen::Vector2i test2i = Eigen::Vector2i::Random();
  const Eigen::Vector3i test3i = Eigen::Vector3i::Random();
  const Eigen::Vector4i test4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi testXi = Eigen::VectorXi::Random(6);
  EXPECT_EQ(alg->Size(test2i), 2);
  EXPECT_EQ(alg->Size(test3i), 3);
  EXPECT_EQ(alg->Size(test4i), 4);
  EXPECT_EQ(alg->Size(testXi), 6);
}

TEST(EigenVectorAlgebraTests, AccessElement) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test Eigen::Vector's
  const Eigen::Vector2d test2d(2.0, 4.0);
  const Eigen::Vector2f test2f(-32.0, 12.5);
  const Eigen::Vector2i test2i(-5, 4);
  for( unsigned int i=0; i<2; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(test2d, i)), test2d(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(test2f, i)), test2f(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(test2i, i)), test2i(i));
  }

  const Eigen::Vector3d test3d(3.0, 4.0, -1.3);
  const Eigen::Vector3f test3f(-33.0, 13.5, 0.1);
  const Eigen::Vector3i test3i(-5, 4, 8);
  for( unsigned int i=0; i<3; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(test3d, i)), test3d(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(test3f, i)), test3f(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(test3i, i)), test3i(i));
  }

  const Eigen::Vector4d test4d(4.0, 4.0, -1.4, 1.3);
  const Eigen::Vector4f test4f(-44.0, 14.5, 0.1, 3.6);
  const Eigen::Vector4i test4i(-5, 4, 8, 3);
  for( unsigned int i=0; i<4; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(test4d, i)), test4d(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(test4f, i)), test4f(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(test4i, i)), test4i(i));
  }

  const Eigen::VectorXd testXd = Eigen::VectorXd::Random(13);
  const Eigen::VectorXf testXf = Eigen::VectorXf::Random(13);
  const Eigen::VectorXi testXi = Eigen::VectorXi::Random(13);
  for( unsigned int i=0; i<13; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(testXd, i)), testXd(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(testXf, i)), testXf(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(testXi, i)), testXi(i));
  }
}

TEST(EigenVectorAlgebraTests, Zero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Vector zero
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector2d const&>(alg->Zero(typeid(Eigen::Vector2d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector2f const&>(alg->Zero(typeid(Eigen::Vector2f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector2i const&>(alg->Zero(typeid(Eigen::Vector2i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector3d const&>(alg->Zero(typeid(Eigen::Vector3d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector3f const&>(alg->Zero(typeid(Eigen::Vector3f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector3i const&>(alg->Zero(typeid(Eigen::Vector3i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector4d const&>(alg->Zero(typeid(Eigen::Vector4d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector4f const&>(alg->Zero(typeid(Eigen::Vector4f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector4i const&>(alg->Zero(typeid(Eigen::Vector4i))).norm(), 0.0);

  const Eigen::VectorXd vecd = boost::any_cast<Eigen::VectorXd const&>(alg->Zero(typeid(Eigen::VectorXd), 13));
  const Eigen::VectorXf vecf = boost::any_cast<Eigen::VectorXf const&>(alg->Zero(typeid(Eigen::VectorXf), 13));
  const Eigen::VectorXi veci = boost::any_cast<Eigen::VectorXi const&>(alg->Zero(typeid(Eigen::VectorXi), 13));
  EXPECT_DOUBLE_EQ(vecd.norm(), 0.0);
  EXPECT_EQ(vecd.size(), 13);
  EXPECT_DOUBLE_EQ(vecf.norm(), 0.0);
  EXPECT_EQ(vecf.size(), 13);
  EXPECT_DOUBLE_EQ(veci.norm(), 0.0);
  EXPECT_EQ(veci.size(), 13);
}

TEST(EigenVectorAlgebraTests, IsZero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  EXPECT_TRUE(alg->IsZero((Eigen::Vector2d)Eigen::Vector2d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector2d)Eigen::Vector2d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector2f)Eigen::Vector2f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector2f)Eigen::Vector2f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector2i)Eigen::Vector2i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector2i)Eigen::Vector2i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Vector3d)Eigen::Vector3d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector3d)Eigen::Vector3d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector3f)Eigen::Vector3f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector3f)Eigen::Vector3f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector3i)Eigen::Vector3i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector3i)Eigen::Vector3i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Vector4d)Eigen::Vector4d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector4d)Eigen::Vector4d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector4f)Eigen::Vector4f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector4f)Eigen::Vector4f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector4i)Eigen::Vector4i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector4i)Eigen::Vector4i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::VectorXd)Eigen::VectorXd::Zero(21)));
  EXPECT_FALSE(alg->IsZero((Eigen::VectorXd)Eigen::VectorXd::Random(38)));
  EXPECT_TRUE(alg->IsZero((Eigen::VectorXf)Eigen::VectorXf::Zero(19)));
  EXPECT_FALSE(alg->IsZero((Eigen::VectorXf)Eigen::VectorXf::Random(22)));
  EXPECT_TRUE(alg->IsZero((Eigen::VectorXi)Eigen::VectorXi::Zero(36)));
  EXPECT_FALSE(alg->IsZero((Eigen::VectorXi)Eigen::VectorXi::Random(5)));
}

TEST(EigenVectorAlgebraTests, Identity) {
    auto alg = std::shared_ptr<AnyAlgebra>();

    // Eigen::Vector
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d>(alg->Identity(typeid(Eigen::Vector2d))).array()==Eigen::Matrix2d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f>(alg->Identity(typeid(Eigen::Vector2f))).array()==Eigen::Matrix2f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i>(alg->Identity(typeid(Eigen::Vector2i))).array()==Eigen::Matrix2i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d>(alg->Identity(typeid(Eigen::Vector3d))).array()==Eigen::Matrix3d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f>(alg->Identity(typeid(Eigen::Vector3f))).array()==Eigen::Matrix3f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i>(alg->Identity(typeid(Eigen::Vector3i))).array()==Eigen::Matrix3i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d>(alg->Identity(typeid(Eigen::Vector4d))).array()==Eigen::Matrix4d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f>(alg->Identity(typeid(Eigen::Vector4f))).array()==Eigen::Matrix4f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i>(alg->Identity(typeid(Eigen::Vector4i))).array()==Eigen::Matrix4i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd>(alg->Identity(typeid(Eigen::VectorXd), 4, 9)).array()==Eigen::MatrixXd::Identity(4, 9).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf>(alg->Identity(typeid(Eigen::VectorXf), 43, 21)).array()==Eigen::MatrixXf::Identity(43, 21).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi>(alg->Identity(typeid(Eigen::VectorXi), 23, 23)).array()==Eigen::MatrixXi::Identity(23, 23).array()).all());
}

TEST(EigenVectorAlgebraTests, Norm) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test Eigen::Vector norm
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vecXd = Eigen::VectorXd::Random(8);
  EXPECT_DOUBLE_EQ(alg->Norm(vec2d), vec2d.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vec3d), vec3d.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vec4d), vec4d.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vecXd), vecXd.norm());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vecXf = Eigen::VectorXf::Random(8);
  EXPECT_DOUBLE_EQ(alg->Norm(vec2f), vec2f.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vec3f), vec3f.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vec4f), vec4f.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vecXf), vecXf.norm());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vecXi = Eigen::VectorXi::Random(8);
  EXPECT_DOUBLE_EQ(alg->Norm(vec2i), vec2i.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vec3i), vec3i.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vec4i), vec4i.norm());
  EXPECT_DOUBLE_EQ(alg->Norm(vecXi), vecXi.norm());
}

TEST(EigenVectorAlgebraTests, InnerProduct) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  const Eigen::Vector2d vec2d(2.0, 3.0);
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  const Eigen::Vector3d vec3d(2.0, 3.0, 4.0);
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  const Eigen::Vector4d vec4d(2.0, 3.0, 4.0, 5.0);
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_EQ(alg->InnerProduct(vec2d, vec2d), vec2d.dot(vec2d));
  EXPECT_EQ(alg->InnerProduct(vec2d, vec2Xd), vec2d.dot(vec2Xd));
  EXPECT_EQ(alg->InnerProduct(vec3d, vec3d), vec3d.dot(vec3d));
  EXPECT_EQ(alg->InnerProduct(vec3d, vec3Xd), vec3d.dot(vec3Xd));
  EXPECT_EQ(alg->InnerProduct(vec4d, vec4d), vec4d.dot(vec4d));
  EXPECT_EQ(alg->InnerProduct(vec4d, vec4Xd), vec4d.dot(vec4Xd));
  EXPECT_EQ(alg->InnerProduct(vec3Xd, vec3Xd), vec3Xd.dot(vec3Xd));

  const Eigen::Vector2f vec2f(2.0, 3.0);
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  const Eigen::Vector3f vec3f(2.0, 3.0, 4.0);
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  const Eigen::Vector4f vec4f(2.0, 3.0, 4.0, 5.0);
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_EQ(alg->InnerProduct(vec2f, vec2f), vec2f.dot(vec2f));
  EXPECT_EQ(alg->InnerProduct(vec2f, vec2Xf), vec2f.dot(vec2Xf));
  EXPECT_EQ(alg->InnerProduct(vec3f, vec3f), vec3f.dot(vec3f));
  EXPECT_EQ(alg->InnerProduct(vec3f, vec3Xf), vec3f.dot(vec3Xf));
  EXPECT_EQ(alg->InnerProduct(vec4f, vec4f), vec4f.dot(vec4f));
  EXPECT_EQ(alg->InnerProduct(vec4f, vec4Xf), vec4f.dot(vec4Xf));
  EXPECT_EQ(alg->InnerProduct(vec3Xf, vec3Xf), vec3Xf.dot(vec3Xf));

  const Eigen::Vector2i vec2i(2, 3);
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  const Eigen::Vector3i vec3i(2, 3, 4);
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  const Eigen::Vector4i vec4i(2, 3, 4, 5);
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_EQ(alg->InnerProduct(vec2i, vec2i), vec2i.dot(vec2i));
  EXPECT_EQ(alg->InnerProduct(vec2i, vec2Xi), vec2i.dot(vec2Xi));
  EXPECT_EQ(alg->InnerProduct(vec3i, vec3i), vec3i.dot(vec3i));
  EXPECT_EQ(alg->InnerProduct(vec3i, vec3Xi), vec3i.dot(vec3Xi));
  EXPECT_EQ(alg->InnerProduct(vec4i, vec4i), vec4i.dot(vec4i));
  EXPECT_EQ(alg->InnerProduct(vec4i, vec4Xi), vec4i.dot(vec4Xi));
  EXPECT_EQ(alg->InnerProduct(vec3Xi, vec3Xi), vec3Xi.dot(vec3Xi));
}

TEST(EigenVectorAlgebraTests, OuterProduct) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  const Eigen::Vector2d vec2d(2.0, 3.0);
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  const Eigen::Vector3d vec3d(2.0, 3.0, 4.0);
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  const Eigen::Vector4d vec4d(2.0, 3.0, 4.0, 5.0);
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix2d>(alg->OuterProduct(vec2d, vec2d))-vec2d*vec2d.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec2d, vec2Xd))-vec2d*vec2Xd.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec2Xd, vec2d))-vec2Xd*vec2d.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix3d>(alg->OuterProduct(vec3d, vec3d))-vec3d*vec3d.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec3d, vec3Xd))-vec3d*vec3Xd.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec3Xd, vec3d))-vec3Xd*vec3d.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix4d>(alg->OuterProduct(vec4d, vec4d))-vec4d*vec4d.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec4d, vec4Xd))-vec4d*vec4Xd.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec4Xd, vec4d))-vec4Xd*vec4d.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXd>(alg->OuterProduct(vec3Xd, vec3Xd))-vec3Xd*vec3Xd.transpose()).norm(), 0.0);

  const Eigen::Vector2f vec2f(2.0, 3.0);
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  const Eigen::Vector3f vec3f(2.0, 3.0, 4.0);
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  const Eigen::Vector4f vec4f(2.0, 3.0, 4.0, 5.0);
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix2f>(alg->OuterProduct(vec2f, vec2f))-vec2f*vec2f.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec2f, vec2Xf))-vec2f*vec2Xf.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec2Xf, vec2f))-vec2Xf*vec2f.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix3f>(alg->OuterProduct(vec3f, vec3f))-vec3f*vec3f.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec3f, vec3Xf))-vec3f*vec3Xf.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec3Xf, vec3f))-vec3Xf*vec3f.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix4f>(alg->OuterProduct(vec4f, vec4f))-vec4f*vec4f.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec4f, vec4Xf))-vec4f*vec4Xf.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec4Xf, vec4f))-vec4Xf*vec4f.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXf>(alg->OuterProduct(vec3Xf, vec3Xf))-vec3Xf*vec3Xf.transpose()).norm(), 0.0);

  const Eigen::Vector2i vec2i(2, 3);
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  const Eigen::Vector3i vec3i(2, 3, 4);
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  const Eigen::Vector4i vec4i(2, 3, 4, 5);
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix2i>(alg->OuterProduct(vec2i, vec2i))-vec2i*vec2i.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec2i, vec2Xi))-vec2i*vec2Xi.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec2Xi, vec2i))-vec2Xi*vec2i.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix3i>(alg->OuterProduct(vec3i, vec3i))-vec3i*vec3i.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec3i, vec3Xi))-vec3i*vec3Xi.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec3Xi, vec3i))-vec3Xi*vec3i.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::Matrix4i>(alg->OuterProduct(vec4i, vec4i))-vec4i*vec4i.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec4i, vec4Xi))-vec4i*vec4Xi.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec4Xi, vec4i))-vec4Xi*vec4i.transpose()).norm(), 0.0);
  EXPECT_DOUBLE_EQ((boost::any_cast<Eigen::MatrixXi>(alg->OuterProduct(vec3Xi, vec3Xi))-vec3Xi*vec3Xi.transpose()).norm(), 0.0);
}

TEST(EigenVectorAlgebraTests, LogDeterminate) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random().cwiseAbs();
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random().cwiseAbs();
  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random().cwiseAbs();
  const Eigen::VectorXd vecXd = Eigen::VectorXd::Random(13).cwiseAbs();
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vec2d)), vec2d.array().log().sum());
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vec3d)), vec3d.array().log().sum());
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vec4d)), vec4d.array().log().sum());
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vecXd)), vecXd.array().log().sum());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random().cwiseAbs();
  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random().cwiseAbs();
  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random().cwiseAbs();
  const Eigen::VectorXf vecXf = Eigen::VectorXf::Random(13).cwiseAbs();
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vec2f)), vec2f.array().log().sum());
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vec3f)), vec3f.array().log().sum());
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vec4f)), vec4f.array().log().sum());
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->LogDeterminate(vecXf)), vecXf.array().log().sum());
}

TEST(EigenVectorAlgebraTests, SquareRoot) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random().cwiseAbs();
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random().cwiseAbs();
  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random().cwiseAbs();
  const Eigen::VectorXd vecXd = Eigen::VectorXd::Random(13).cwiseAbs();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->SquareRoot(vec2d)).array()==(Eigen::Array2d)vec2d.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->SquareRoot(vec3d)).array()==(Eigen::Array3d)vec3d.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->SquareRoot(vec4d)).array()==(Eigen::Array4d)vec4d.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->SquareRoot(vecXd)).array()==(Eigen::ArrayXd)vecXd.cwiseSqrt()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random().cwiseAbs();
  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random().cwiseAbs();
  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random().cwiseAbs();
  const Eigen::VectorXf vecXf = Eigen::VectorXf::Random(13).cwiseAbs();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->SquareRoot(vec2f)).array()==(Eigen::Array2f)vec2f.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->SquareRoot(vec3f)).array()==(Eigen::Array3f)vec3f.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->SquareRoot(vec4f)).array()==(Eigen::Array4f)vec4f.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->SquareRoot(vecXf)).array()==(Eigen::ArrayXf)vecXf.cwiseSqrt()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random().cwiseAbs();
  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random().cwiseAbs();
  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random().cwiseAbs();
  const Eigen::VectorXi vecXi = Eigen::VectorXi::Random(13).cwiseAbs();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->SquareRoot(vec2i)).array()==(Eigen::Array2i)vec2i.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->SquareRoot(vec3i)).array()==(Eigen::Array3i)vec3i.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->SquareRoot(vec4i)).array()==(Eigen::Array4i)vec4i.cwiseSqrt()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->SquareRoot(vecXi)).array()==(Eigen::ArrayXi)vecXi.cwiseSqrt()).all());
}

TEST(EigenVectorAlgebraTests, Add) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Add(vec2d, vec2d)).array()==(vec2d+vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Add(vec2d, vec2Xd)).array()==(vec2d+vec2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vec2Xd, vec2d)).array()==(vec2d+vec2Xd).array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Add(vec2f, vec2f)).array()==(vec2f+vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Add(vec2f, vec2Xf)).array()==(vec2f+vec2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vec2Xf, vec2f)).array()==(vec2f+vec2Xf).array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Add(vec2i, vec2i)).array()==(vec2i+vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Add(vec2i, vec2Xi)).array()==(vec2i+vec2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(vec2Xi, vec2i)).array()==(vec2i+vec2Xi).array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Add(vec3d, vec3d)).array()==(vec3d+vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Add(vec3d, vec3Xd)).array()==(vec3d+vec3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vec3Xd, vec3d)).array()==(vec3d+vec3Xd).array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Add(vec3f, vec3f)).array()==(vec3f+vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Add(vec3f, vec3Xf)).array()==(vec3f+vec3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vec3Xf, vec3f)).array()==(vec3f+vec3Xf).array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Add(vec3i, vec3i)).array()==(vec3i+vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Add(vec3i, vec3Xi)).array()==(vec3i+vec3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(vec3Xi, vec3i)).array()==(vec3i+vec3Xi).array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Add(vec4d, vec4d)).array()==(vec4d+vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Add(vec4d, vec4Xd)).array()==(vec4d+vec4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vec4Xd, vec4d)).array()==(vec4d+vec4Xd).array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Add(vec4f, vec4f)).array()==(vec4f+vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Add(vec4f, vec4Xf)).array()==(vec4f+vec4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vec4Xf, vec4f)).array()==(vec4f+vec4Xf).array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Add(vec4i, vec4i)).array()==(vec4i+vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Add(vec4i, vec4Xi)).array()==(vec4i+vec4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(vec4Xi, vec4i)).array()==(vec4i+vec4Xi).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(8);
  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vecd, vecd)).array()==(vecd+vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vecf, vecf)).array()==(vecf+vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(veci, veci)).array()==(veci+veci).array()).all());
}

TEST(EigenVectorAlgebraTests, Subtract) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Subtract(vec2d, vec2d)).array()==(vec2d-vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Subtract(vec2d, vec2Xd)).array()==(vec2d-vec2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vec2Xd, vec2d)).array()==(vec2Xd-vec2d).array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Subtract(vec2f, vec2f)).array()==(vec2f-vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Subtract(vec2f, vec2Xf)).array()==(vec2f-vec2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vec2Xf, vec2f)).array()==(vec2Xf-vec2f).array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Subtract(vec2i, vec2i)).array()==(vec2i-vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Subtract(vec2i, vec2Xi)).array()==(vec2i-vec2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(vec2Xi, vec2i)).array()==(vec2Xi-vec2i).array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Subtract(vec3d, vec3d)).array()==(vec3d-vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Subtract(vec3d, vec3Xd)).array()==(vec3d-vec3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vec3Xd, vec3d)).array()==(vec3Xd-vec3d).array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Subtract(vec3f, vec3f)).array()==(vec3f-vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Subtract(vec3f, vec3Xf)).array()==(vec3f-vec3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vec3Xf, vec3f)).array()==(vec3Xf-vec3f).array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Subtract(vec3i, vec3i)).array()==(vec3i-vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Subtract(vec3i, vec3Xi)).array()==(vec3i-vec3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(vec3Xi, vec3i)).array()==(vec3Xi-vec3i).array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Subtract(vec4d, vec4d)).array()==(vec4d-vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Subtract(vec4d, vec4Xd)).array()==(vec4d-vec4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vec4Xd, vec4d)).array()==(vec4Xd-vec4d).array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Subtract(vec4f, vec4f)).array()==(vec4f-vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Subtract(vec4f, vec4Xf)).array()==(vec4f-vec4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vec4Xf, vec4f)).array()==(vec4Xf-vec4f).array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Subtract(vec4i, vec4i)).array()==(vec4i-vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Subtract(vec4i, vec4Xi)).array()==(vec4i-vec4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(vec4Xi, vec4i)).array()==(vec4Xi-vec4i).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(8);
  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vecd, vecd)).array()==(vecd-vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vecf, vecf)).array()==(vecf-vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(veci, veci)).array()==(veci-veci).array()).all());
}

TEST(EigenVectorAlgebraTests, Multiply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Multiply(xd, vec2d)).array()==(xd*vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Multiply(vec2d, xd)).array()==(xd*vec2d).array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Multiply(xf, vec2f)).array()==(xf*vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Multiply(vec2f, xf)).array()==(xf*vec2f).array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(xi, vec2i)).array()==(xi*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(vec2i, xi)).array()==(xi*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(xui, vec2i)).array()==(xui*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(vec2i, xui)).array()==(xui*vec2i).array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Multiply(xd, vec3d)).array()==(xd*vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Multiply(vec3d, xd)).array()==(xd*vec3d).array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Multiply(xf, vec3f)).array()==(xf*vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Multiply(vec3f, xf)).array()==(xf*vec3f).array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(xi, vec3i)).array()==(xi*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(vec3i, xi)).array()==(xi*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(xui, vec3i)).array()==(xui*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(vec3i, xui)).array()==(xui*vec3i).array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Multiply(xd, vec4d)).array()==(xd*vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Multiply(vec4d, xd)).array()==(xd*vec4d).array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Multiply(xf, vec4f)).array()==(xf*vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Multiply(vec4f, xf)).array()==(xf*vec4f).array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(xi, vec4i)).array()==(xi*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(vec4i, xi)).array()==(xi*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(xui, vec4i)).array()==(xui*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(vec4i, xui)).array()==(xui*vec4i).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(12);
  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Multiply(xd, vecd)).array()==(xd*vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Multiply(vecd, xd)).array()==(xd*vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Multiply(xf, vecf)).array()==(xf*vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Multiply(vecf, xf)).array()==(xf*vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(xi, veci)).array()==(xi*veci).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(veci, xi)).array()==(xi*veci).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(xui, veci)).array()==(xui*veci).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(veci, xui)).array()==(xui*veci).array()).all());
}

TEST(EigenVectorAlgebraTests, Apply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  const Eigen::Vector2d vm_2d2d = vec2d.asDiagonal()*vec2d;
  const Eigen::Vector2d vm_2d2Xd = vec2d.asDiagonal()*vec2Xd;
  const Eigen::Vector2d vm_2Xd2d = vec2Xd.asDiagonal()*vec2d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->Apply(vec2d, vec2d)).array()==vm_2d2d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vec2d, vec2Xd)).array()==vm_2d2Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->Apply(vec2Xd, vec2d)).array()==vm_2Xd2d.array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  const Eigen::Vector2f vm_2f2f = vec2f.asDiagonal()*vec2f;
  const Eigen::Vector2f vm_2f2Xf = vec2f.asDiagonal()*vec2Xf;
  const Eigen::Vector2f vm_2Xf2f = vec2Xf.asDiagonal()*vec2f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->Apply(vec2f, vec2f)).array()==vm_2f2f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vec2f, vec2Xf)).array()==vm_2f2Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->Apply(vec2Xf, vec2f)).array()==vm_2Xf2f.array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  const Eigen::Vector2i vm_2i2i = vec2i.asDiagonal()*vec2i;
  const Eigen::Vector2i vm_2i2Xi = vec2i.asDiagonal()*vec2Xi;
  const Eigen::Vector2i vm_2Xi2i = vec2Xi.asDiagonal()*vec2i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->Apply(vec2i, vec2i)).array()==vm_2i2i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(vec2i, vec2Xi)).array()==vm_2i2Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->Apply(vec2Xi, vec2i)).array()==vm_2Xi2i.array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  const Eigen::Vector3d vm_3d3d = vec3d.asDiagonal()*vec3d;
  const Eigen::Vector3d vm_3d3Xd = vec3d.asDiagonal()*vec3Xd;
  const Eigen::Vector3d vm_3Xd3d = vec3Xd.asDiagonal()*vec3d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->Apply(vec3d, vec3d)).array()==vm_3d3d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vec3d, vec3Xd)).array()==vm_3d3Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->Apply(vec3Xd, vec3d)).array()==vm_3Xd3d.array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  const Eigen::Vector3f vm_3f3f = vec3f.asDiagonal()*vec3f;
  const Eigen::Vector3f vm_3f3Xf = vec3f.asDiagonal()*vec3Xf;
  const Eigen::Vector3f vm_3Xf3f = vec3Xf.asDiagonal()*vec3f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->Apply(vec3f, vec3f)).array()==vm_3f3f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vec3f, vec3Xf)).array()==vm_3f3Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->Apply(vec3Xf, vec3f)).array()==vm_3Xf3f.array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  const Eigen::Vector3i vm_3i3i = vec3i.asDiagonal()*vec3i;
  const Eigen::Vector3i vm_3i3Xi = vec3i.asDiagonal()*vec3Xi;
  const Eigen::Vector3i vm_3Xi3i = vec3Xi.asDiagonal()*vec3i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->Apply(vec3i, vec3i)).array()==vm_3i3i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(vec3i, vec3Xi)).array()==vm_3i3Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->Apply(vec3Xi, vec3i)).array()==vm_3Xi3i.array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  const Eigen::Vector4d vm_4d4d = vec4d.asDiagonal()*vec4d;
  const Eigen::Vector4d vm_4d4Xd = vec4d.asDiagonal()*vec4Xd;
  const Eigen::Vector4d vm_4Xd4d = vec4Xd.asDiagonal()*vec4d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->Apply(vec4d, vec4d)).array()==vm_4d4d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vec4d, vec4Xd)).array()==vm_4d4Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->Apply(vec4Xd, vec4d)).array()==vm_4Xd4d.array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  const Eigen::Vector4f vm_4f4f = vec4f.asDiagonal()*vec4f;
  const Eigen::Vector4f vm_4f4Xf = vec4f.asDiagonal()*vec4Xf;
  const Eigen::Vector4f vm_4Xf4f = vec4Xf.asDiagonal()*vec4f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->Apply(vec4f, vec4f)).array()==vm_4f4f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vec4f, vec4Xf)).array()==vm_4f4Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->Apply(vec4Xf, vec4f)).array()==vm_4Xf4f.array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  const Eigen::Vector4i vm_4i4i = vec4i.asDiagonal()*vec4i;
  const Eigen::Vector4i vm_4i4Xi = vec4i.asDiagonal()*vec4Xi;
  const Eigen::Vector4i vm_4Xi4i = vec4Xi.asDiagonal()*vec4i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->Apply(vec4i, vec4i)).array()==vm_4i4i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(vec4i, vec4Xi)).array()==vm_4i4Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->Apply(vec4Xi, vec4i)).array()==vm_4Xi4i.array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(13);
  const Eigen::VectorXd vm_dd = vecd.asDiagonal()*vecd;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vecd, vecd)).array()==vm_dd.array()).all());

  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(13);
  const Eigen::VectorXf vm_ff = vecf.asDiagonal()*vecf;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vecf, vecf)).array()==vm_ff.array()).all());

  const Eigen::VectorXi veci = Eigen::VectorXi::Random(13);
  const Eigen::VectorXi vm_ii = veci.asDiagonal()*veci;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(veci, veci)).array()==vm_ii.array()).all());
}

TEST(EigenVectorAlgebraTests, ApplyInverse) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  const Eigen::Vector2d vm_2d2d = (1.0/vec2d.array()).matrix().asDiagonal()*vec2d;
  const Eigen::Vector2d vm_2d2Xd = (1.0/vec2d.array()).matrix().asDiagonal()*vec2Xd;
  const Eigen::Vector2d vm_2Xd2d = (1.0/vec2Xd.array()).matrix().asDiagonal()*vec2d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->ApplyInverse(vec2d, vec2d)).array()==vm_2d2d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vec2d, vec2Xd)).array()==vm_2d2Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->ApplyInverse(vec2Xd, vec2d)).array()==vm_2Xd2d.array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  const Eigen::Vector2f vm_2f2f = (1.0/vec2f.array()).matrix().asDiagonal()*vec2f;
  const Eigen::Vector2f vm_2f2Xf = (1.0/vec2f.array()).matrix().asDiagonal()*vec2Xf;
  const Eigen::Vector2f vm_2Xf2f = (1.0/vec2Xf.array()).matrix().asDiagonal()*vec2f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->ApplyInverse(vec2f, vec2f)).array()==vm_2f2f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vec2f, vec2Xf)).array()==vm_2f2Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->ApplyInverse(vec2Xf, vec2f)).array()==vm_2Xf2f.array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  const Eigen::Vector2i vm_2i2i = (1/vec2i.array()).matrix().asDiagonal()*vec2i;
  const Eigen::Vector2i vm_2i2Xi = (1/vec2i.array()).matrix().asDiagonal()*vec2Xi;
  const Eigen::Vector2i vm_2Xi2i = (1/vec2Xi.array()).matrix().asDiagonal()*vec2i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->ApplyInverse(vec2i, vec2i)).array()==vm_2i2i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(vec2i, vec2Xi)).array()==vm_2i2Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->ApplyInverse(vec2Xi, vec2i)).array()==vm_2Xi2i.array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  const Eigen::Vector3d vm_3d3d = (1.0/vec3d.array()).matrix().asDiagonal()*vec3d;
  const Eigen::Vector3d vm_3d3Xd = (1.0/vec3d.array()).matrix().asDiagonal()*vec3Xd;
  const Eigen::Vector3d vm_3Xd3d = (1.0/vec3Xd.array()).matrix().asDiagonal()*vec3d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->ApplyInverse(vec3d, vec3d)).array()==vm_3d3d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vec3d, vec3Xd)).array()==vm_3d3Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->ApplyInverse(vec3Xd, vec3d)).array()==vm_3Xd3d.array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  const Eigen::Vector3f vm_3f3f = (1.0/vec3f.array()).matrix().asDiagonal()*vec3f;
  const Eigen::Vector3f vm_3f3Xf = (1.0/vec3f.array()).matrix().asDiagonal()*vec3Xf;
  const Eigen::Vector3f vm_3Xf3f = (1.0/vec3Xf.array()).matrix().asDiagonal()*vec3f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->ApplyInverse(vec3f, vec3f)).array()==vm_3f3f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vec3f, vec3Xf)).array()==vm_3f3Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->ApplyInverse(vec3Xf, vec3f)).array()==vm_3Xf3f.array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  const Eigen::Vector3i vm_3i3i = (1/vec3i.array()).matrix().asDiagonal()*vec3i;
  const Eigen::Vector3i vm_3i3Xi = (1/vec3i.array()).matrix().asDiagonal()*vec3Xi;
  const Eigen::Vector3i vm_3Xi3i = (1/vec3Xi.array()).matrix().asDiagonal()*vec3i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->ApplyInverse(vec3i, vec3i)).array()==vm_3i3i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(vec3i, vec3Xi)).array()==vm_3i3Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->ApplyInverse(vec3Xi, vec3i)).array()==vm_3Xi3i.array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  const Eigen::Vector4d vm_4d4d = (1.0/vec4d.array()).matrix().asDiagonal()*vec4d;
  const Eigen::Vector4d vm_4d4Xd = (1.0/vec4d.array()).matrix().asDiagonal()*vec4Xd;
  const Eigen::Vector4d vm_4Xd4d = (1.0/vec4Xd.array()).matrix().asDiagonal()*vec4d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->ApplyInverse(vec4d, vec4d)).array()==vm_4d4d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vec4d, vec4Xd)).array()==vm_4d4Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->ApplyInverse(vec4Xd, vec4d)).array()==vm_4Xd4d.array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  const Eigen::Vector4f vm_4f4f = (1.0/vec4f.array()).matrix().asDiagonal()*vec4f;
  const Eigen::Vector4f vm_4f4Xf = (1.0/vec4f.array()).matrix().asDiagonal()*vec4Xf;
  const Eigen::Vector4f vm_4Xf4f = (1.0/vec4Xf.array()).matrix().asDiagonal()*vec4f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->ApplyInverse(vec4f, vec4f)).array()==vm_4f4f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vec4f, vec4Xf)).array()==vm_4f4Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->ApplyInverse(vec4Xf, vec4f)).array()==vm_4Xf4f.array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  const Eigen::Vector4i vm_4i4i = (1/vec4i.array()).matrix().asDiagonal()*vec4i;
  const Eigen::Vector4i vm_4i4Xi = (1/vec4i.array()).matrix().asDiagonal()*vec4Xi;
  const Eigen::Vector4i vm_4Xi4i = (1/vec4Xi.array()).matrix().asDiagonal()*vec4i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->ApplyInverse(vec4i, vec4i)).array()==vm_4i4i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(vec4i, vec4Xi)).array()==vm_4i4Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->ApplyInverse(vec4Xi, vec4i)).array()==vm_4Xi4i.array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(13);
  const Eigen::VectorXd vm_dd = (1.0/vecd.array()).matrix().asDiagonal()*vecd;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vecd, vecd)).array()==vm_dd.array()).all());

  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(13);
  const Eigen::VectorXf vm_ff = (1.0/vecf.array()).matrix().asDiagonal()*vecf;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vecf, vecf)).array()==vm_ff.array()).all());

  const Eigen::VectorXi veci = Eigen::VectorXi::Random(13);
  const Eigen::VectorXi vm_ii = (1/veci.array()).matrix().asDiagonal()*veci;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(veci, veci)).array()==vm_ii.array()).all());
}
