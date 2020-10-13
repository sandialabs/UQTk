#include "gtest/gtest.h"

#include <cstdio>

#include "MUQ/Utilities/HDF5/HDF5File.h"
#include "EigenTestUtils.h"

using namespace muq::Utilities;

class HDF5FileTest : public::testing::Test {
protected:

  virtual void SetUp() override {
    filename = "test.h5";

    // create a new file
    hdf5file = std::make_shared<HDF5File>(filename);
  }
  
  virtual void TearDown() override {
    // close the file. Note: this is not strictly necessary since the file is automatically closed when HDF5File is destroyed.
    hdf5file->Close();

    // remove the file
    std::remove(filename.c_str());
  }

  /// The hdf5 file (muq wrapper)
  std::shared_ptr<HDF5File> hdf5file;

  /// the temporary filename
  std::string filename;
};

TEST_F(HDF5FileTest, ReadWrite) {
    
  // create a vector and matrix
  Eigen::VectorXd testVec = Eigen::VectorXd::Random(5);
  Eigen::MatrixXi testMat = Eigen::MatrixXi::Ones(3, 2);

  // write the vector and matrix to file
  hdf5file->WriteMatrix("/vector", testVec);
  hdf5file->WriteMatrix<int>("/matrix", testMat);

  // the vector and matrix should exist
  EXPECT_TRUE(hdf5file->DoesDataSetExist("/vector"));
  EXPECT_TRUE(hdf5file->DoesDataSetExist("/matrix"));

  // get the dimensions of the vector and matrix
  Eigen::VectorXi testVecSize = hdf5file->GetDataSetSize("/vector");
  Eigen::VectorXi testMatSize = hdf5file->GetDataSetSize("/matrix");

  // make sure they are the correct size
  EXPECT_EQ(5, testVecSize(0));
  EXPECT_EQ(1, testVecSize(1));
  EXPECT_EQ(3, testMatSize(0));
  EXPECT_EQ(2, testMatSize(1));

  // read the vector and matrix 
  Eigen::VectorXd truthVec = hdf5file->ReadMatrix("/vector");
  Eigen::MatrixXi truthMat = hdf5file->ReadMatrix<int>("/matrix");

  // Read onlyt part of the vector
  Eigen::VectorXd partVec = hdf5file->ReadPartialMatrix("/vector", 0, 0, 2, 1);
  EXPECT_EQ(2, partVec.rows());
  EXPECT_EQ(1, partVec.cols());
  EXPECT_DOUBLE_EQ(testVec(0), partVec(0));
  EXPECT_DOUBLE_EQ(testVec(1), partVec(1));
  
  // they should match
  EXPECT_PRED_FORMAT3(MatrixApproxEqual, truthVec, testVec, 1e-10);
  EXPECT_PRED_FORMAT2(MatrixEqual, truthMat, testMat);

  // create a different vector and matrix (change sizes)
  testVec = Eigen::VectorXd::Random(3);
  testMat = Eigen::MatrixXi::Ones(5, 8);

  // write the vector and matrix to file
  hdf5file->WriteMatrix("/vector", testVec);
  hdf5file->WriteMatrix<int>("/matrix", testMat);

  // the vector and matrix should exist
  EXPECT_TRUE(hdf5file->DoesDataSetExist("/vector"));
  EXPECT_TRUE(hdf5file->DoesDataSetExist("/matrix"));

  // get the dimensions of the vector and matrix
  testVecSize = hdf5file->GetDataSetSize("/vector");
  testMatSize = hdf5file->GetDataSetSize("/matrix");

  // make sure they are the correct size
  EXPECT_EQ(3, testVecSize(0));
  EXPECT_EQ(1, testVecSize(1));
  EXPECT_EQ(5, testMatSize(0));
  EXPECT_EQ(8, testMatSize(1));

  // read the vector and matrix 
  truthVec = hdf5file->ReadMatrix("/vector");
  truthMat = hdf5file->ReadMatrix<int>("/matrix");

  // they should match
  EXPECT_PRED_FORMAT3(MatrixApproxEqual, truthVec, testVec, 1e-10);
  EXPECT_PRED_FORMAT2(MatrixEqual, truthMat, testMat);
}

TEST_F(HDF5FileTest, ChildList)
{
  // create some groups
  hdf5file->CreateGroup("/a/crazy/path");
  hdf5file->CreateGroup("/here/is/a");
  hdf5file->CreateGroup("/here/is/b/");

  // create a dataset
  const Eigen::VectorXd testVec = Eigen::VectorXd::Random(3);
  hdf5file->WriteMatrix("/here/is/a/vector", testVec);

  std::vector<std::string> children = hdf5file->GetChildren();
  EXPECT_EQ(2,children.size());
  EXPECT_EQ("a",children.at(0));
  EXPECT_EQ("here",children.at(1));

  children = hdf5file->GetChildren("/here");
  EXPECT_EQ(1,children.size());
  EXPECT_EQ("is",children.at(0));

  children = hdf5file->GetChildren("/a/");
  EXPECT_EQ(1,children.size());
  EXPECT_EQ("crazy", children.at(0));
  
  children = hdf5file->GetChildren("/here/is");
  EXPECT_EQ(2,children.size());
  EXPECT_EQ("a",children.at(0));
  EXPECT_EQ("b",children.at(1)); 
}

TEST_F(HDF5FileTest, MultiLevelGroup) {
  // some data
  const Eigen::VectorXd testVec = Eigen::VectorXd::Random(3);
  const Eigen::MatrixXd testMat = Eigen::MatrixXd::Random(3, 2);

  //some tests on an empty file
  EXPECT_FALSE(hdf5file->DoesDataSetExist("/here/is/a/vector"));
  EXPECT_FALSE(hdf5file->DoesGroupExist("/here/"));
  EXPECT_FALSE(hdf5file->DoesGroupExist("/here/is"));
  EXPECT_FALSE(hdf5file->DoesGroupExist("/here/is/a"));
  EXPECT_FALSE(hdf5file->DoesGroupExist("/here/is/b"));

  // create some groups
  hdf5file->CreateGroup("/a/crazy/path");
  hdf5file->CreateGroup("/here/is/a");
  hdf5file->CreateGroup("/here/is/b/");

  EXPECT_TRUE(hdf5file->IsGroup("/a/crazy/path"));
  EXPECT_FALSE(hdf5file->IsGroup("/a/not/so/crazy/path"));
  
  // make sure the groups were created
  EXPECT_TRUE(hdf5file->DoesGroupExist("/a/crazy/path/"));
  EXPECT_TRUE(hdf5file->DoesGroupExist("/here/"));
  EXPECT_TRUE(hdf5file->DoesGroupExist("/here/is"));
  EXPECT_TRUE(hdf5file->DoesGroupExist("/here/is/a/"));
  EXPECT_TRUE(hdf5file->DoesGroupExist("/here/is/b"));

  // write the vector/matrix to file
  hdf5file->WriteMatrix("/here/is/a/vector", testVec);
  hdf5file->WriteMatrix("/here/is/a/matrix", testMat);

  EXPECT_TRUE(hdf5file->IsDataSet("/here/is/a/vector"));
  EXPECT_TRUE(hdf5file->IsDataSet("/here/is/a/matrix"));
  EXPECT_FALSE(hdf5file->IsGroup("/here/is/a/vector"));
  EXPECT_FALSE(hdf5file->IsGroup("/here/is/a/vector"));
  
  // make sure it is there
  ASSERT_TRUE(hdf5file->DoesDataSetExist("/here/is/a/vector"));
  ASSERT_TRUE(hdf5file->DoesDataSetExist("/here/is/a/matrix"));

  // read it back
  const Eigen::VectorXd truthVec = hdf5file->ReadMatrix("/here/is/a/vector");
  const Eigen::MatrixXd truthMat = hdf5file->ReadMatrix("/here/is/a/matrix");

  // it should match
  EXPECT_PRED_FORMAT3(MatrixApproxEqual, truthVec, testVec, 1e-10);
  EXPECT_PRED_FORMAT3(MatrixApproxEqual, truthMat, testMat, 1e-10);
}

TEST_F(HDF5FileTest, Attributes) {
  // a group name
  std::string groupName0 = "/a/crazy//path";
  std::string groupName1 = "/another/crazy//path";

  // create the group
  hdf5file->CreateGroup(groupName0);

  // make sure it was created (or not)
  EXPECT_TRUE(hdf5file->DoesGroupExist(groupName0));
  EXPECT_FALSE(hdf5file->DoesGroupExist(groupName1));

  // write some attributes to file
  hdf5file->WriteScalarAttribute<double>(groupName0, "derivedAtt2", 9);
  hdf5file->WriteScalarAttribute<double>(groupName0, "doubleAtt", 1.234);
  hdf5file->WriteScalarAttribute<float>(groupName1, "floatAtt", 1.234);
  hdf5file->WriteScalarAttribute<int>(groupName1, "intAtt", -11);
  hdf5file->WriteScalarAttribute<unsigned>(groupName1, "unsignedAtt", 11);
  hdf5file->WriteStringAttribute(groupName0, "stringAtt", "myTestString");

  // make sure it was created
  EXPECT_TRUE(hdf5file->DoesGroupExist(groupName0));
  EXPECT_TRUE(hdf5file->DoesGroupExist(groupName1));

  // close the file
  hdf5file->Close();

  // open the file
  hdf5file->Open(filename);

  EXPECT_EQ(9, hdf5file->GetScalarAttribute<double>(groupName0, "derivedAtt2"));
  EXPECT_DOUBLE_EQ(1.234, hdf5file->GetScalarAttribute<double>(groupName0, "doubleAtt"));
  EXPECT_FLOAT_EQ(1.234, hdf5file->GetScalarAttribute<float>(groupName1, "floatAtt"));
  EXPECT_EQ(-11, hdf5file->GetScalarAttribute<int>(groupName1, "intAtt"));
  EXPECT_EQ(11, hdf5file->GetScalarAttribute<unsigned>(groupName1, "unsignedAtt"));
  EXPECT_EQ("myTestString", hdf5file->GetStringAttribute(groupName0, "stringAtt"));
}

TEST_F(HDF5FileTest, WritePartial1) {
  // create the full matrix of nan's
  const Eigen::MatrixXi nanMat = Eigen::MatrixXi::Constant(5, 6, -1);
  Eigen::MatrixXi fullMat = Eigen::MatrixXi::Constant(5, 6, -1);
  
  // write the matrix to the file to create the full dataset 
  hdf5file->WriteMatrix("/path/to/matrix", nanMat);

  // create the patrial matrices 
  const Eigen::VectorXi topleft = Eigen::VectorXi::Ones(3);
  const Eigen::MatrixXi topright = 2 * Eigen::MatrixXi::Ones(3, 5);
  const Eigen::MatrixXi bottomleft = 3 * Eigen::MatrixXi::Ones(2, 2);
  const Eigen::MatrixXi bottomright = 4 * Eigen::MatrixXi::Ones(2, 4);

  // write the partial matrices to file
  hdf5file->WritePartialMatrix("/path/to/matrix", topleft, 0, 0);
  fullMat.block<3, 1>(0, 0) = topleft;

  EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));

  // write the partial matrices to file
  hdf5file->WritePartialMatrix("/path/to/matrix", topright, 0, 1);
  fullMat.block<3, 5>(0, 1) = topright;

  EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));

  // write the partial matrices to file
  hdf5file->WritePartialMatrix("/path/to/matrix", bottomleft, 3, 0);
  fullMat.block<2, 2>(3, 0) = bottomleft;

  EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));

  // write the partial matrices to file
  hdf5file->WritePartialMatrix("/path/to/matrix", bottomright, 3, 2);
  fullMat.block<2, 4>(3, 2) = bottomright;

  EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));
}

TEST_F(HDF5FileTest, WritePartial2) {

    // create the full matrix
    hdf5file->CreateDataset<int>("/path/to/matrix", 5, 6);
    
    Eigen::MatrixXi fullMat = Eigen::MatrixXi::Zero(5,6);
    hdf5file->WriteMatrix("/path/to/matrix", fullMat);
    
    // create the patrial matrices
    const Eigen::VectorXi topleft = Eigen::VectorXi::Ones(3);
    const Eigen::MatrixXi topright = 2 * Eigen::MatrixXi::Ones(3, 5);
    const Eigen::MatrixXi bottomleft = 3 * Eigen::MatrixXi::Ones(2, 2);
    const Eigen::MatrixXi bottomright = 4 * Eigen::MatrixXi::Ones(2, 4);
    
    // write the partial matrices to file
    hdf5file->WritePartialMatrix("/path/to/matrix", topleft, 0, 0);
    fullMat.block<3, 1>(0, 0) = topleft;
    
    EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));
    
    // write the partial matrices to file
    hdf5file->WritePartialMatrix("/path/to/matrix", topright, 0, 1);
    fullMat.block<3, 5>(0, 1) = topright;
    
    EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));
    
    // write the partial matrices to file
    hdf5file->WritePartialMatrix("/path/to/matrix", bottomleft, 3, 0);
    fullMat.block<2, 2>(3, 0) = bottomleft;
    
    EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));
    
    // write the partial matrices to file
    hdf5file->WritePartialMatrix("/path/to/matrix", bottomright, 3, 2);
    fullMat.block<2, 4>(3, 2) = bottomright;
    
    EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));
    
    // copy the second column into the first and make sure it worked
    hdf5file->WritePartialMatrix<int,Eigen::Dynamic,1>("/path/to/matrix", fullMat.col(1), 0,0);
    fullMat.col(0) = fullMat.col(1);
    EXPECT_PRED_FORMAT2(MatrixEqual, fullMat, hdf5file->ReadMatrix<int>("/path/to/matrix"));
    
}

TEST_F(HDF5FileTest, MergeFiles) {
    std::string tempName = "otherFile.h5";

    // create a new file
    auto otherFile = std::make_shared<HDF5File>(tempName);

    // test long vector and large matrix to read and write to file 
    const Eigen::VectorXd testVec = Eigen::VectorXd::Random(10);
    const Eigen::MatrixXd testMat = Eigen::MatrixXd::Random(13, 10);
    
    // write the vector and matrix to seperate files
    hdf5file->WriteMatrix("/vector", testVec);
    otherFile->WriteMatrix("/matrix", testMat);
    
    hdf5file->MergeFile(otherFile);

    // close the other file
    otherFile->Close();

    // remove the file
    std::remove(tempName.c_str());

    // read the vector and matrix from file
    const Eigen::VectorXd truthVec = hdf5file->ReadMatrix("/vector");
    const Eigen::MatrixXd truthMat = hdf5file->ReadMatrix("/matrix");
    
    // it should match
    EXPECT_PRED_FORMAT3(MatrixApproxEqual, truthVec, testVec, 1e-10);
    EXPECT_PRED_FORMAT3(MatrixApproxEqual, truthMat, testMat, 1e-10);
}

