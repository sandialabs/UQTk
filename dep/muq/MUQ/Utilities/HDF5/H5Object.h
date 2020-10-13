#ifndef H5OBJECT_H
#define H5OBJECT_H

#include "MUQ/Utilities/HDF5/HDF5File.h"
#include "MUQ/Utilities/HDF5/BlockDataset.h"
#include "MUQ/Utilities/HDF5/Attributes.h"

#include <Eigen/Core>
#include <string>
#include <map>
#include <iostream>
#include <memory>
#include <type_traits>

namespace muq
{

namespace Utilities
{


class H5Object;

/// Recursively add children to create an HDF5 file hierarchy
H5Object AddChildren(std::shared_ptr<HDF5File>        file,
       std::string               const& groupName);

/// Open an HDF5 file and return the root object
/**
   @brief Primary way to open an HDF5file.  It returns the root group.
   @details
@code
// Open the HDF5 file for reading and writing
auto f = muq::Utilities::OpenFile("Data.h5");

// Create a new group
f.CreateGroup("/NewGroup");

// Add an attribute to the group
f["/NewGroup"].attrs["Some Metadata"] = "Created with MH5!";

// Store a vector in the group
f["/NewGroup/Ones"] = Eigen::VectorXd::Ones(10);

// Store another vector in a different group.  The group is automatically generated
f["/AnotherGroup/Zeros"] = Eigen::MatrixXd::Zero(5,5);
@endcode

 */
H5Object OpenFile(std::string const& filename);


class H5Object
{
    friend H5Object muq::Utilities::AddChildren(std::shared_ptr<HDF5File> file, std::string const& path);

public:

    H5Object(){};

    H5Object(std::shared_ptr<HDF5File>        file_,
       std::string               const& path_,
       bool                             isDataset_) : file(file_),
                                                            attrs(file_, path_),
 	                                                    path(path_),
                                                      isDataset(isDataset_){};

    typedef std::function<void(boost::any const&, H5Object& )> AnyWriterType;
    typedef std::unordered_map<std::type_index, AnyWriterType> AnyWriterMapType;

    static std::shared_ptr<AnyWriterMapType> GetAnyWriterMap();

    H5Object& operator=(boost::any const& val);


    // Use this templated function for arithmetic types
    template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
    H5Object& operator=(ScalarType val)
    {
      assert(path.length()>0);
      if(isDataset)
      {
        Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> temp(1,1);
        temp(0,0) = val;
        file->WriteMatrix(path, temp);
      }
      else
      {
        assert(false);
      }

        return *this;
    };

    // Use this templated function for non-arithmetic types
    template<typename Derived>
    H5Object& operator=(Eigen::DenseBase<Derived> const& val)
    {
      return (*this)=val.eval();
    };

    template<typename ScalarType, int fixedRows, int fixedCols>
    H5Object& operator=(Eigen::Matrix<ScalarType, fixedRows, fixedCols> const& val)
    {
      assert(path.length()>0);
      if(isDataset)
      {
        file->WriteMatrix(path, val);
      }
      else
      {
        assert(false);
      }

        return *this;
    };

    template<typename scalarType, int rows, int cols>
    operator Eigen::Matrix<scalarType,rows,cols>()
    {
      return eval<scalarType,rows,cols>();
    }

    /**
        @brief Performs a deep copy of the datasets and groups in otherObj into this object.
    */
    H5Object& operator=(H5Object const& otherObj);

    template<typename scalarType=double, int rows=Eigen::Dynamic, int cols=Eigen::Dynamic>
    Eigen::Matrix<scalarType,rows,cols> eval()
    {
      if(isDataset)
      {
        return file->ReadMatrix<scalarType,rows,cols>(path).template cast<scalarType>();
      }
      else
      {
        assert(false);
      }
    }
    //////////////////////////////////////////////////////
    // Create groups and datsets
    //////////////////////////////////////////////////////
    H5Object& CreatePlaceholder(std::string const& grpName);
    H5Object& CreateGroup(std::string const& grpName);

    /** Creates a dataset with a particular size and type. */
    template<typename ScalarType>
    H5Object& CreateDataset(std::string const& setName,
                          unsigned int rows,
                          unsigned int cols=0)
    {
      H5Object& temp = CreatePlaceholder(setName);
      temp.file->CreateDataset<ScalarType>(temp.path, rows, cols);
      return temp;
    }

    //////////////////////////////////////////////////////
    // Accessors
    //////////////////////////////////////////////////////
    H5Object& operator[](std::string const& targetPath);

    //const H5Object& operator[](std::string const& path) const;

    BlockDataset block(unsigned startRow, unsigned startCol, unsigned numRows, unsigned numCols) const;

    BlockDataset topLeftCorner(unsigned numRows, unsigned numCols) const;
    BlockDataset bottomLeftCorner(unsigned numRows, unsigned numCols) const;

    BlockDataset topRightCorner(unsigned numRows, unsigned numCols) const;
    BlockDataset bottomRightCorner(unsigned numRows, unsigned numCols) const;

    BlockDataset topRows(unsigned numRows) const;
    BlockDataset bottomRows(unsigned numRows) const;

    BlockDataset leftCols(unsigned numCols) const;
    BlockDataset rightCols(unsigned numCols) const;

    BlockDataset col(unsigned col) const;

    BlockDataset row(unsigned row) const;

    BlockDataset segment(unsigned startInd, unsigned numInds) const;
    BlockDataset head(unsigned numInds) const;
    BlockDataset tail(unsigned numInds) const;

    unsigned rows() const;

    unsigned cols() const;

    unsigned size() const;

    double operator()(int i) const;

    double operator()(int i, int j) const;

    //////////////////////////////////////////////////////
    // Utility Functions
    //////////////////////////////////////////////////////

    void Flush();

    void Print(std::string prefix = "") const;

    std::shared_ptr<HDF5File> file;

    AttributeList attrs;

private:

    /** @brief Creates an exact copy.
        @details Equivalent to the default assignment operator.  Does not copy datasets,
                 only the path, children, and isDataset values from another H5Object.
    */
    void ExactCopy(H5Object const& otherObj);


    // Copy the other objects content into the current dataset
    void DeepCopy(H5Object const& otherObj);

    std::string path;

    std::map<std::string, H5Object> children;

    bool isDataset;

};


#ifndef REGISTER_HDF5OBJECT_ANYTYPE
#define REGISTER_HDF5OBJECT_ANYTYPE(REGNAME, NAME) static auto regHDF ##REGNAME \
        = muq::Utilities::H5Object::GetAnyWriterMap()->insert(std::make_pair(std::type_index(typeid(NAME)), muq::Utilities::AnyWriter<NAME>() ));

#endif



} // namespace Utilities
} // namespace muq


#endif
