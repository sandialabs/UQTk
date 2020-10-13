#include <gtest/gtest.h>

#include "MUQ/Utilities/HDF5/H5Object.h"




class DatasetTest : public::testing::Test {
protected:

    virtual void SetUp() override {
    	filename = "test.h5";
    	hdf5file = std::make_shared<muq::Utilities::HDF5File>(filename);

    	// create some groups
    	hdf5file->CreateGroup("/a/crazy/path");
    	hdf5file->CreateGroup("/here/is/a");
    	hdf5file->CreateGroup("/here/is/b/");

    	// create a dataset
    	const Eigen::VectorXd testVec = Eigen::VectorXd::Random(3);
    	hdf5file->WriteMatrix("/here/is/a/vector", testVec);
    }

    virtual void TearDown() override {
    	hdf5file->Close();
    	std::remove(filename.c_str());
    }

    std::shared_ptr<muq::Utilities::HDF5File> hdf5file;
    std::string filename;
};

TEST_F(DatasetTest, Groups)
{
    muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

    auto g = f.CreateGroup("/something");

    g["/first"]  = 1.0;
    g["/second"] = 2.0;

    double val = f["/something/first"](0);
    EXPECT_DOUBLE_EQ(1.0, val);

    val = g["/second"](0);
    EXPECT_DOUBLE_EQ(2.0, val);
}


TEST_F(DatasetTest, DoubleSetGet)
{
    muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

    f["/something/first"] = 1.0;
    f["/something/second"] = 2.0;

    const int numThird = 10;
    f["/something/else/third"] = Eigen::VectorXd::LinSpaced(numThird,0,1);

    double firstVal = f["/something/first"](0);
    double secondVal = f["/something/second"](0);

    Eigen::VectorXd third = f["/something/else/third"];

    EXPECT_DOUBLE_EQ(1.0, firstVal);
    EXPECT_DOUBLE_EQ(2.0, secondVal);

    for(int i=0; i<numThird; ++i)
	EXPECT_DOUBLE_EQ(double(i)/double(numThird-1), third(i));
}

TEST_F(DatasetTest, CopyDataset)
{
  muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

  f["/something/first"] = 1.0;

  f["/something/second"] = f["/something/first"];

  EXPECT_EQ(1, f["/something/second"].rows());
  EXPECT_EQ(1, f["/something/second"].cols());
  EXPECT_DOUBLE_EQ(1.0, f["/something/second"](0));

}

TEST_F(DatasetTest, CopyDatasetWithAttribute)
{
  muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

  f["/something/first"] = 1.0;
  f["/something/first"].attrs["meta1"] = "Units";
  f["/something/first"].attrs["meta2"] = 10.0;

  f["/something/second"] = f["/something/first"];

  EXPECT_EQ(1, f["/something/second"].rows());
  EXPECT_EQ(1, f["/something/second"].cols());
  EXPECT_DOUBLE_EQ(1.0, f["/something/second"](0));

  std::string meta1 = f["/something/second"].attrs["meta1"];
  EXPECT_EQ(std::string("Units"), meta1);
  EXPECT_DOUBLE_EQ(10.0, f["/something/second"].attrs["meta2"]);
}

TEST_F(DatasetTest, CopyGroup)
{
  muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

  f["/something/first"] = 1.0;
  f["/something/first"].attrs["meta1"] = "Units";
  f["/something/first"].attrs["meta2"] = 10.0;

  f["/something/second"] = f["/something/first"];

  f["/somethingelse"] = f["/something"];

  EXPECT_EQ(1, f["/somethingelse/first"].rows());
  EXPECT_EQ(1, f["/somethingelse/first"].cols());
  EXPECT_DOUBLE_EQ(1.0, f["/somethingelse/first"](0));

  std::string meta1 = f["/somethingelse/first"].attrs["meta1"];
  EXPECT_EQ(std::string("Units"), meta1);
  EXPECT_DOUBLE_EQ(10.0, f["/somethingelse/first"].attrs["meta2"]);

  EXPECT_EQ(1, f["/somethingelse/second"].rows());
  EXPECT_EQ(1, f["/somethingelse/second"].cols());
  EXPECT_DOUBLE_EQ(1.0, f["/somethingelse/second"](0));

  meta1 = std::string(f["/somethingelse/second"].attrs["meta1"]);
  EXPECT_EQ(std::string("Units"), meta1);
  EXPECT_DOUBLE_EQ(10.0, f["/somethingelse/second"].attrs["meta2"]);
}

TEST_F(DatasetTest, AttributeSetGet)
{
    muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

    f["/something/first"] = 1.0;
    f["/something/second"] = 2.0;

    {
	auto& attrs = f["/something"].attrs;
	attrs["Meta1"] = 10.0;
	attrs["Meta2"] = "This is an attribute.";

	auto& attrs2 = f["/something/first"].attrs;
	attrs2["OtherVec"] = Eigen::VectorXd::Ones(10).eval();

	// Make sure contents have been written to file
	f.Flush();
    }

    {
	double val1 = f["/something"].attrs["Meta1"];
	EXPECT_DOUBLE_EQ(10.0, val1);

	std::string val2 = f["/something"].attrs["Meta2"];
	EXPECT_EQ("This is an attribute.", val2);


	Eigen::VectorXd temp = f["/something/first"].attrs["OtherVec"];
	EXPECT_EQ(10, temp.size());
	for(int i=0; i<10; ++i)
	    EXPECT_DOUBLE_EQ(1.0, temp(i));
    }
}


TEST_F(DatasetTest, BlockReadOps)
{

    Eigen::MatrixXd A(3,3);
    A << 2, 1, 0,
	 1, 2, 1,
	 0, 2, 3;

    Eigen::VectorXd b(3);
    b << 1,2,3;

    {
	muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

	f["/A"] = A;
	f["/b"] = b;

	EXPECT_EQ(3, f["/A"].rows());
	EXPECT_EQ(3, f["/A"].cols());
	EXPECT_EQ(9, f["/A"].size());

	EXPECT_EQ(3, f["/b"].rows());
	EXPECT_EQ(1, f["/b"].cols());
	EXPECT_EQ(3, f["/b"].size());

	f.Flush();
    }

    {

	muq::Utilities::H5Object f = muq::Utilities::OpenFile("test.h5");

	Eigen::MatrixXd A2 = f["/A"].col(2);
	EXPECT_EQ(3,A2.size());
	EXPECT_EQ(0, A2(0));
	EXPECT_EQ(1, A2(1));
	EXPECT_EQ(3, A2(2));

	A2 = f["/A"].row(1);
	EXPECT_EQ(3,A2.size());
	EXPECT_EQ(1, A2(0));
	EXPECT_EQ(2, A2(1));
	EXPECT_EQ(1, A2(2));

	A2 = f["/A"].block(1,1,2,2);
	EXPECT_EQ(2, A2.rows());
	EXPECT_EQ(2, A2.cols());

	EXPECT_EQ(2, A2(0,0));
	EXPECT_EQ(1, A2(0,1));
	EXPECT_EQ(2, A2(1,0));
	EXPECT_EQ(3, A2(1,1));


	Eigen::VectorXd b2 = f["/b"].segment(0,2);

	EXPECT_EQ(2, b2.size());
	EXPECT_EQ(1, b2(0));
	EXPECT_EQ(2, b2(1));
    }

}

TEST_F(DatasetTest, BlockWriteOps)
{

    muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

    Eigen::MatrixXd A(3,3);
    A << 2, 1, 0,
	       1, 2, 1,
	       0, 2, 3;

    f["/A"] = A;

    Eigen::VectorXd newCol(3);
    newCol << -1, -2, -3;

    f["/A"].col(1) = newCol;

    Eigen::MatrixXd A2 = f["/A"];
    EXPECT_DOUBLE_EQ(2, A2(0,0));
    EXPECT_DOUBLE_EQ(-1, A2(0,1));
    EXPECT_DOUBLE_EQ(-2, A2(1,1));
    EXPECT_DOUBLE_EQ(-3, A2(2,1));
}

TEST_F(DatasetTest, IntSetGet)
{
    muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

    f["/something/first"] = int(1);

    const int numThird = 10;
    f["/something/else/third"] = 10*Eigen::VectorXi::Ones(numThird);

    int firstVal = f["/something/first"](0);

    Eigen::VectorXi third = f["/something/else/third"];

    EXPECT_DOUBLE_EQ(1, firstVal);

    for(int i=0; i<numThird; ++i)
	EXPECT_DOUBLE_EQ(10, third(i));


}

TEST_F(DatasetTest, CloseOpen)
{
    {
	muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

	f["/something/first"] = 1.0;
	f["/something/second"] = 2.0;

	f.file->Close();
    }

    {
	muq::Utilities::H5Object f = muq::Utilities::OpenFile("test.h5");

	double result = f["/something/first"](0);
	EXPECT_DOUBLE_EQ(1.0, result);

        result = f["/something/second"](0);
	EXPECT_DOUBLE_EQ(2.0, result);
    }

}

TEST_F(DatasetTest, GroupSizeFail)
{
  muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

  EXPECT_THROW(f.rows(), std::runtime_error);
  EXPECT_THROW(f.cols(), std::runtime_error);
  EXPECT_THROW(f.size(), std::runtime_error);
}

TEST_F(DatasetTest, CreateDataset)
{
  muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

  const int numRows = 10;
  const int numCols = 20;

  auto ds = f.CreateDataset<double>("/My/New/Dataset", numRows, numCols);

  EXPECT_EQ(numRows,ds.rows());
  EXPECT_EQ(numCols,ds.cols());

  Eigen::MatrixXd temp = f["/My/New/Dataset"];
  for(int j=0; j<numCols; ++j){
    for(int i=0; i<numRows; ++i)
      EXPECT_DOUBLE_EQ(0.0,temp(i,j));
  }

  ds.col(0) = Eigen::VectorXd::Ones(numRows).eval();

  temp = f["/My/New/Dataset"];
  for(int i=0; i<numRows; ++i)
    EXPECT_DOUBLE_EQ(1.0,temp(i,0));

  for(int j=1; j<numCols; ++j){
    for(int i=0; i<numRows; ++i)
      EXPECT_DOUBLE_EQ(0.0,temp(i,j));
  }
}

TEST_F(DatasetTest, FromH5File)
{

    muq::Utilities::H5Object f = AddChildren(hdf5file, "/");

    f["/here/is/a/vector"] = 2.0;
    f["/here/is/a/vector"] = Eigen::VectorXd::Ones(10);

    //std::cout << f["/here/is/a/vector"](0) << std::endl;
}
