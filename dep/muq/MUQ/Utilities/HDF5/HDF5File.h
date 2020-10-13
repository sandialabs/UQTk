#ifndef HDF5FILE_H_
#define HDF5FILE_H_

#include <Eigen/Core>

#include "MUQ/Utilities/HDF5/HDF5Types.h"
#include "MUQ/Utilities/HDF5/PathTools.h"

#include <vector>
#include <memory>

#include <hdf5.h>
#include <hdf5_hl.h>


namespace muq
{

namespace Utilities
{

    /// A wrapper around HDF5 (and PHDF5) to handle file I/O
    class HDF5File : public std::enable_shared_from_this<HDF5File> {
    public:

	/// Open or create the file
	/**
	If the input file exists, it is opened; if the does not exist, it is created.
	@param[in] filename_ The name (including the path) of the file
	@param[in] write The write node that does all of the write-to-file (defaults to 0 --- only important if MPI is enabled)
	@param[in] multipleFiles True: each rank opens filename seperately and write to its own file (each rank should have a unique filename); False (default): rank 0 opens filename and the other send data to it to write to file.
	*/
	HDF5File(std::string const& filename_);

	/// If HDF5File is destroyed, the file should be closed
	virtual ~HDF5File();

	/// Opens or creates the file
	/**
	   If the input file exists, it is opened; if the does not exist, it is created.
	   @param[in] filename_ The name (including the path) of the file
	   \return The file id of the opened or created file
	*/
	void Open(std::string const& filename_);

	/// Close the file
	/**
	   If the file is not open (HDF5File::fileID is <= 0) then this function does nothing, otherwise the file is closed.
	*/
	void Close();

        /// Copy the contents of one dataset into another
        /**
           Copies the dataset at location srcName in the file srcFile into the location destName of this file.  This is a deep copy.
         */
        void Copy(std::string const& destName, std::shared_ptr<HDF5File> srcFile, std::string const& srcName);

	/// Write data set to file
	/**
	   Write a new data set to file.  If the dataset (and corresponding group) does not exist it will be created.  If the dataset does exist it is deleted and overwritten.
	   <h3> In serial (MPI not enabled): </h3>
	   The dataset is written to the user-provided data set name.
	   <h3> In parallel (MPI enabled): </h3>
	   In parallel, multiple datasets are created and (potentially different but the same size) data is written to each one.  The names must be unique
	   @param[in] name The name of the data set
	   @param[in] dataset The matrix or vector to save to file (the dataset)
	   @param[in] singleProc True if only a single processor (asumed to be the write parameter) is calling this function (defaults to false, does not matter if MPI is not enabled)
	*/
	template<typename scalarType, int fixedRows, int fixedCols>
	void WriteMatrix(std::string const& name,
			 Eigen::Matrix<scalarType, fixedRows, fixedCols> const& dataset) {

	    // data set name must begin with '/'
	    if( name.at(0)!='/' ) {
		std::cerr << std::endl << "ERROR: Paths in the HDF5 file must start with a forward slash (/)" << std::endl << "\tHDF5File::WriteMatrix(std::string const&, Eigen::Matrix<scalarType, fixedRows, fixedCols> const&)" << std::endl << std::endl;
		assert(name.at(0) == '/');
	    }

	    // make sure the file is open
	    assert(fileID>0);

	    // set the maximum size of the data set to be unlimited
	    const hsize_t maxdim[2] = {H5S_UNLIMITED, H5S_UNLIMITED};

	    const hsize_t dimsf[2] = {(hsize_t)dataset.rows(), (hsize_t)dataset.cols()};

	    // create the dataspace for the dataset.
	    const hid_t filespace = H5Screate_simple(2, dimsf, maxdim);
	    assert(filespace>0);

	    // the dataset we are writing to
	    hid_t dataID;
	    if( DoesDataSetExist(name) ) {

		// open the data
		dataID = H5Dopen(fileID, name.c_str(), H5P_DEFAULT);

		// get the dataspace dimension size
		hsize_t* dims = (hsize_t*)malloc(2*sizeof(hsize_t));
		H5Sget_simple_extent_dims(H5Dget_space(dataID), dims, nullptr);

		if( (dimsf[0]!=dims[0]) || (dimsf[1]!=dims[1]) ) { // if the data size changes ...
		    // Extend the dataset
		    H5Dset_extent(dataID, dimsf);
		}

	    } else {

		// get the parent path
		std::string parentPath = GetParentPath(name);

		// make sure the parent exists (and create it if it does not)
		CreateGroup(parentPath);

		// modify data set creation properties --- enable chunking
		const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);

		// set chunk size to be (1,1) -- may be a bad idea of the data set is extremely large
		hsize_t chunk[2] = {100,100};
		if(fixedRows>0)
		{
		    chunk[0] = std::min<hsize_t>(chunk[0],fixedRows);
		}
		if(fixedCols>0)
		{
		    chunk[1] = std::min<hsize_t>(chunk[1],fixedRows);
		}
		H5Pset_chunk(prop, 2, chunk);

		// create a dataset for each process
		dataID = H5Dcreate(fileID, name.c_str(), HDF5_Type<scalarType>::GetFlag(), filespace, H5P_DEFAULT, prop, H5P_DEFAULT);

		// close the creation properties
		H5Pclose(prop);
	    }

	    // constant reference (avoid copying) to the data transpose because of column versus row major conventions
	    const Eigen::Matrix<scalarType, fixedRows, fixedCols>& dataset_ = dataset.transpose();

	    // write the data to file
	    H5Dwrite(dataID, HDF5_Type<scalarType>::GetFlag(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_.data());

	    // close the filespace
	    H5Sclose(filespace);

	    // close the dataset
	    H5Dclose(dataID);

	}

	/// Write part of a data set to file
	/**
	   Write part of a data set to file.  The data set must already exist and be large enough to hold the data given to this function.
	   @param[in] name The name of the data set
	   @param[in] data The matrix or vector to save to file (part of the dataset)
	   @param[in] row The row to start writing the data
	   @param[in] col The column to start writting the data
	   @param[in] singleProc True if only a single processor (asumed to be the write parameter) is calling this function (defaults to false, does not matter if MPI is not enabled)
	*/
	template<typename scalarType, int fixedRows, int fixedCols>
        void WritePartialMatrix(std::string const& name,
				Eigen::Matrix<scalarType, fixedRows, fixedCols> const& data,
				unsigned int const row, unsigned int const col)
	{
	    // data set name must begin with '/'
	    if( name.at(0)!='/' ) {
		std::cerr << std::endl << "ERROR: Paths in the HDF5 file must start with a forward slash (/)" << std::endl << "\tHDF5File::WriteMatrix(std::string const&, Eigen::Matrix<scalarType, fixedRows, fixedCols> const&)" << std::endl << std::endl;
		assert(name.at(0) == '/');
	    }

	    // make sure the file is open
	    assert(fileID>0);

	    // how many elements to move in each dimension
	    const hsize_t stride[2] = {1, 1};

	    // how many blocks to move in each dimension
	    const hsize_t count[2] = {1, 1};

	    // make sure the data set exists
	    if( !DoesDataSetExist(name) ) {
		std::cerr << std::endl << "ERROR: Dataset " << name << " does not exsts." << std::endl << std::endl;
		assert(DoesDataSetExist(name));
	    }

	    // open the data set
	    const hid_t dataID = H5Dopen(fileID, name.c_str(), H5P_DEFAULT);

	    // get the space in the file
	    const hid_t fullspace = H5Dget_space(dataID);

	    // the location in the dataset where we start writing
	    const hsize_t start[2] = {row, col};

	    // the size of each block
	    const hsize_t block[2] = {(hsize_t)data.rows(), (hsize_t)data.cols()};

	    // get the hyperslab (the subspace needed for the partial write)
	    H5Sselect_hyperslab(fullspace, H5S_SELECT_SET, start, stride, count, block);

	    // get the dataspace for the partial write
	    const hsize_t limitedBlock = H5Screate_simple(2, block, block);

	    // constant reference (prevent copying) to transpose data
	    const Eigen::Matrix<scalarType, fixedRows, fixedCols>& data_ = data.transpose();

	    // write to file
	    auto write = H5Dwrite(dataID,
				  HDF5_Type<scalarType>::GetFlag(),
				  limitedBlock,
				  fullspace,
				  H5P_DEFAULT,
				  data_.data());

	    // close the spaces
	    H5Sclose(limitedBlock);
	    H5Sclose(fullspace);
	    H5Dclose(dataID);
	}

	/// Read a dataset from the file
	/**
	   This function will return an empty matrix if the data set does not exist.

	   Unless the matrix is full of doubles, the scalarType template parameter should be set to the expected type.  Also, fixed size matrices can be returned by setting the fixedRows and/or fixedCols template parameters to something other than Eigen::Dynamic.
	   @param[in] name The name of the dataset.
	*/
	template<typename scalarType = double, int fixedRows = Eigen::Dynamic, int fixedCols = Eigen::Dynamic>
        Eigen::Matrix<scalarType, fixedRows, fixedCols> ReadMatrix(std::string const& name) const {

	    // make sure the file is open
	    assert(fileID>0);

	    if( !DoesDataSetExist(name) ) { // if the data does not exist
		      // ... return an empty matrix.
          return Eigen::Matrix<scalarType, fixedRows, fixedCols>();
	    }

	    // get the dataset ID
	    const hid_t dataset = H5Dopen2(fileID, name.c_str(), H5P_DEFAULT);

	    // get the space associated with this data set
	    const hid_t filespace = H5Dget_space(dataset);

      const int ndims = H5Sget_simple_extent_ndims(filespace);

	    // get the size of the dataset
	    std::vector<hsize_t> dimsr(ndims);
	    H5Sget_simple_extent_dims(filespace, dimsr.data(), nullptr);

	    // declare the matrix that we will read into
	    Eigen::Matrix<scalarType, Eigen::Dynamic, Eigen::Dynamic> pt;

	    if( ndims==2 ) { // if it is a matrix ...
		      pt.resize(dimsr[1], dimsr[0]);
	    } else { // if it is a vector ...
          if(fixedRows==1){
              pt = Eigen::Matrix<scalarType, 1, fixedCols>(dimsr[0]);
          }else{
              pt = Eigen::Matrix<scalarType, fixedRows, 1>(dimsr[0]);
          }

	    }

	    // get the property list
	    const hid_t prop = H5Dget_create_plist(dataset);

	    // get the memory for this data set
	    const hid_t memspace = H5Screate_simple(ndims, dimsr.data(), nullptr);

	    // read the data
      //std::cout << "Using type " << HDF5_Type<scalarType>::GetFlag() << " vs " << HDF5_Type<double>::GetFlag() << " vs " << H5T_NATIVE_DOUBLE << std::endl;
	    H5Dread(dataset, HDF5_Type<scalarType>::GetFlag(), memspace, filespace, H5P_DEFAULT, pt.data());

	    // close the creation properties
	    H5Pclose(prop);

	    // close the filespace
	    H5Sclose(filespace);

	    // close the memory space
	    H5Sclose(memspace);

	    // close the dataset
	    H5Dclose(dataset);

	    if(ndims==2) {
		    // hdf5 and eigen have different memory convention so we read in the matrix transpose
        if(fixedRows==1){
          return pt;
        }else{
          return pt.transpose();
        }
	    }

	    // return the data
	    return pt;
	}


	template<typename scalarType=double, int fixedRows = Eigen::Dynamic, int fixedCols = Eigen::Dynamic>
	Eigen::Matrix<scalarType, fixedRows, fixedCols> ReadPartialMatrix(std::string const& name,
									  int rowStart,
									  int colStart,
									  int numRows,
									  int numCols) const
	{


	    // make sure the file is open
	    assert(fileID>0);

	    if( !DoesDataSetExist(name) ) { // if the data does not exist
	    	// ... return an empty matrix.
	    	return Eigen::Matrix<scalarType, fixedRows, fixedCols>();
	    }

	    // get the dataset ID
	    const hid_t dataset = H5Dopen2(fileID, name.c_str(), H5P_DEFAULT);

	    // get the space associated with this data set
	    const hid_t dataspace = H5Dget_space(dataset);

	    // get the size of the dataset
	    hsize_t dimsr[2];
	    H5Sget_simple_extent_dims(dataspace, dimsr, nullptr);

	    // declare the matrix that we will read into
	    Eigen::Matrix<scalarType, fixedCols, fixedRows> pt;

	    pt.resize(numCols,numRows); //<- will flip to proper orientation below
	    if( dimsr[1]>0 ) { // if it is a matrix ...
	    	assert(numRows<=dimsr[0]);
	    	assert(numCols<=dimsr[1]);
	    } else { // if it is a vector ...
	    	assert(numRows<dimsr[0]);
	    }

	    const hsize_t start[2]  = {(hsize_t) rowStart, (hsize_t) colStart};
	    const hsize_t count[2]  = {(hsize_t) numRows, (hsize_t) numCols};

	    herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
	    // get the memory for this data set
	    dimsr[0] = numRows;
	    dimsr[1] = numCols;
	    hid_t memspace = H5Screate_simple(2, dimsr, nullptr);

	    // read the data
	    H5Dread(dataset, HDF5_Type<scalarType>::GetFlag(), memspace, dataspace, H5P_DEFAULT, pt.data());

	    // close the filespace
	    H5Sclose(dataspace);

	    // close the memory space
	    H5Sclose(memspace);

	    // close the dataset
	    H5Dclose(dataset);

	    if( dimsr[1]>0 ) {
	    	// hdf5 and eigen have different memory convention so we read in the matrix transpose
	    	return pt.transpose();
	    }

	    // return the data
	    return pt;

	}

	/// Check to see if a group exists
	/**
	   @param[in] name The group name (including the path)
	   \return True if the group exists, false if not
	*/
	bool DoesGroupExist(std::string const& name) const;

	/// Check to see if a data set exists
	/**
	   @param[in] name The dataset name (including the path)
	   \return True if the dataset exists, false if not
	*/
	bool DoesDataSetExist(std::string const& name) const;

	bool IsDataSet(std::string const& name) const;
	bool IsGroup(std::string const& name) const;

	/// Get the size of a dataset (rows,cols)
	/**
	   An empty size zero vector is returned if the dataset could not be found.
	   @param[in] name The name of the data set (include path)
	*/
	Eigen::VectorXi GetDataSetSize(std::string const name) const;

	/// Create a new group in the file
	/**
	   @param[in] name The name of the group
	*/
	void CreateGroup(std::string const& name);



	/// Create a new empty dataset in the file
	/**
	   @param[in] name The name of the dataset to create.  If needed, parent groups will also be created.
	   @param[in] rows The number of rows in the dataset
	   @param[in] cols The number of columns in the dataset.
	   @param[in] singleProc True if only a single processor (asumed to be the write parameter) is calling this function (defaults to false, does not matter if MPI is not enabled)
	*/
	template<typename scalarType>
        void CreateDataset(std::string const& name, int rows, int cols)
	{

            // make sure the file is open
            assert(fileID>0);

            // set the maximum size of the data set to be unlimited
            const hsize_t maxdim[2] = {H5S_UNLIMITED, H5S_UNLIMITED};

            // the sizes of the data set
            const hsize_t dimsf[2] = {(hsize_t)rows, (hsize_t)cols};

            // create the dataspace for the dataset.
            const hid_t filespace = H5Screate_simple(2, dimsf, maxdim);
            assert(filespace>0);

            // the dataset we are writing to
            hid_t dataID;

            if( DoesDataSetExist(name) ) {

                // open the data
                dataID = H5Dopen(fileID, name.c_str(), H5P_DEFAULT);

                // get the dataspace dimension size
                hsize_t dims[2];
                H5Sget_simple_extent_dims(H5Dget_space(dataID), dims, nullptr);

                if( (dimsf[0]!=dims[0]) || (dimsf[1]!=dims[1]) ) { // if the data size changes ...
                    // Extend the dataset
                    H5Dset_extent(dataID, dimsf);
                }

            } else {

                // get the parent path
		std::string parentPath = GetParentPath(name);

                // make sure the parent exists (and create it if it does not)
                CreateGroup(parentPath);

                // modify data set creation properties --- enable chunking
                const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);

                // set chunk size to be (1,1) -- may be a bad idea of the data set is extremely large
                hsize_t chunk[2] = {std::min<hsize_t>(100,rows), std::min<hsize_t>(100,cols)};
                H5Pset_chunk(prop, 2, chunk);

                // create a dataset for each process
                dataID = H5Dcreate(fileID, name.c_str(), HDF5_Type<scalarType>::GetFlag(), filespace, H5P_DEFAULT, prop, H5P_DEFAULT);
                assert(dataID>0);
                
                // close the creation properties
                H5Pclose(prop);
            }

            // close the filespace
            H5Sclose(filespace);

            // close the dataset
            H5Dclose(dataID);
        }

	/// Add a vector attricute to a dataset or group
	/**
	   This function adds a vector valued numeric attribute (such as metadata) to an existing dataset of group.  If the given dataset or group does not exist this function creates it (as a group, not a data set).
	   @param[in] datasetName The name of the dataset or group (starting with a "/" and including the path) to attach the attribute to.
	   @param[in] attributeName The name of the attribute.
	   @param[in] attribute The value of the attribute.
	   @param[in] singleProc True if only a single processor (asumed to be the write parameter) is calling this function (defaults to false, does not matter if MPI is not enabled)
    */
	template<typename scalarType, int fixedRows>
	void WriteVectorAttribute(std::string const& datasetName,
				  std::string const& attributeName,
				  Eigen::Matrix<scalarType, fixedRows, 1> const& attribute)
	{

	    // make sure the file is open
	    assert(fileID>0);

	    // Create the group if necessary
	    if( !DoesDataSetExist(datasetName) || !DoesGroupExist(datasetName) )
		CreateGroup(datasetName);

	    // set the attribute
	    const herr_t set_result = H5LTset_attribute(fileID,
							datasetName.c_str(),
							attributeName.c_str(),
							HDF5_Type<scalarType>::GetFlag(),
							(void*) attribute.data(),
							attribute.rows());
	    assert(set_result>=0);
	}

	/// Add a scalar attribute to a dataset or group
	/**
	   This function adds a scalar valued numeric attribute (such as metadata) to an existing dataset of group.  If the given dataset or group does not exist this function creates it (as a group, not a data set).
	   @param[in] datasetName The name of the dataset or group (starting with a "/" and including the path) that you want to attach this attribute to.
	   @param[in] attributeName The name of the attribute.
	   @param[in] attribute The value of the attribute.
	   @param[in] singleProc True if only a single processor (asumed to be the write parameter) is calling this function (defaults to false, does not matter if MPI is not enabled)
	*/
	template<typename attType = double>
	void WriteScalarAttribute(std::string const& datasetName,
				  std::string const& attributeName,
				  attType const& attribute)
	{
	    // make the scalar a vector
	    Eigen::Matrix<attType, 1, 1> vec(1);
	    vec << attribute;

	    // call the vector attribute version
	    WriteVectorAttribute<attType, 1>(datasetName, attributeName, vec);
	}

	/// Write a string attribute to a dataset or group
	/**
	   This function adds a scalar valued numeric attribute (such as metadata) to an existing dataset of group.  If the given dataset or group does not exist this function creates it (as a group, not a data set).
	   @param[in] datasetName The name of the dataset or group (starting with a "/" and including the path) that you want to attach this attribute to.
	   @param[in] attributeName The name of the attribute.
	   @param[in] attribute The value of the attribute.
	   @param[in] singleProc True if only a single processor (asumed to be the write parameter) is calling this function (defaults to false, does not matter if MPI is not enabled)
	*/
	void WriteStringAttribute(std::string const& datasetName,
				  std::string const& attributeName,
				  std::string const& attribute);

	/// Read a scalar attribute from the HDF5 file
	/**
	   \param[in] datasetName Name of the dataset or group that the attribute is attached to (starting with a "/" and including the path).  If the path does not exist, an assert will fail.
	   \param[in] attributeName The name of the attribute to read.  If the attribute does not exist, the HDF5 library may throw an error.
	*/
	template<typename scalarType>
	scalarType GetScalarAttribute(std::string const& datasetName,
				      std::string const& attributeName) const {

	    // read as a vector attribute
	    Eigen::Matrix<scalarType, Eigen::Dynamic, 1> att = GetVectorAttribute<scalarType>(datasetName, attributeName);
	    assert(att.size() == 1);

	    // return the sclar
	    return att(0);
	}

	/// Read a vector attribute from the HDF5 file
	/**
	   \param[in] datasetName Name of the dataset or group that the attribute is attached to (starting with a "/" and including the path).  If the path does not exist, an assert will fail.
	   \param[in] attributeName The name of the attribute to read.  If the attribute does not exist, the HDF5 library may throw an error.
	*/
	template<typename scalarType>
	Eigen::Matrix<scalarType, Eigen::Dynamic, 1> GetVectorAttribute(std::string const& datasetName,
									std::string const& attributeName) const
	{
	    // make sure the file is open
	    assert(fileID>0);

	    // make sure the dataset exists
	    assert(DoesDataSetExist(datasetName) || DoesGroupExist(datasetName));

	    // get the dimensions of the attribute
	    int rank[1];
	    hsize_t dims[2];
	    auto get_ndim = H5LTget_attribute_ndims(fileID, datasetName.c_str(), attributeName.c_str(), rank);
	    assert(get_ndim >= 0);
	    assert(rank[0]==1);

	    H5T_class_t typeClass;
	    size_t typeSize;
	    auto status = H5LTget_attribute_info(fileID, datasetName.c_str(), attributeName.c_str(), dims, &typeClass, &typeSize);
	    assert(status >= 0);

	    // create the vector to return
	    Eigen::Matrix<scalarType, Eigen::Dynamic, 1> data(dims[0]);

	    // get the attribute

	    auto get_att = H5LTget_attribute(fileID, datasetName.c_str(), attributeName.c_str(), HDF5_Type<scalarType>::GetFlag(), data.data());
	    assert(get_att >= 0);

	    // translate it into and Eigen type
	    //Translater<std::vector<scalarType>, Eigen::Matrix<scalarType, Eigen::Dynamic, 1> > trans(data);
	    return data;
	}

	/// Read a string attribute from the HDF5 file
	/**
	   \param[in] datasetName Name of the dataset or group that the attribute is attached to (starting with a "/" and including the path).  If the path does not exist, an assert will fail.
	   \param[in] attributeName The name of the attribute to read.  If the attribute does not exist, the HDF5 library may throw an error.
	*/
	std::string GetStringAttribute(std::string const& datasetName, std::string const& attributeName) const;


	/// Get a list of immediate children of a group
	std::vector<std::string> GetChildren(std::string base = "/") const;

	/// Flush any data in the HDF5 buffer to the file.
	void FlushFile();

	/// Merge another file into this file
	/**
	   Both files must already be open.
	   @param[in] otherFile The file to merge into this one.
	*/
	void MergeFile(std::shared_ptr<HDF5File> const& otherFile);

	/// The HDF5 file ID
	/**
	   The default value -1 indicates the file is not open.
	*/
	hid_t fileID = -1;

	/// The name of the file
	/**
	   Defaults to an empty file
	*/
	std::string filename = "";

    private:

	bool DoesFileExist(const std::string& name) const;

    };

} // namespace Utilities
} // namespace muq

#endif
