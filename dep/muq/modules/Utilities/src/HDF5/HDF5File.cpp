#include "MUQ/Utilities/HDF5/HDF5File.h"

#include<iostream>
#include<fstream>

using namespace muq::Utilities;

HDF5File::HDF5File(std::string const& filename_){

  // make sure the file is not open
  assert(fileID<0);

  // create (or open) the file
  Open(filename_);
}

HDF5File::~HDF5File() {
  // close the file
  Close();

  // make sure the file is closed
  assert(fileID<=0);
}

bool HDF5File::DoesFileExist(const std::string& name) const {
  std::ifstream f(name.c_str());
  return f.good();
}

void HDF5File::Open(std::string const& filename_) {

    if( fileID>=0 ) { // if a file is already open ...
      // ... close it.
      Close();
    }

    // save the file name;
    filename = filename_;

    // Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

    if( DoesFileExist(filename) ){ // if the file exists ...
      // ... open it.
      fileID = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
    } else { // if the file does not exist ...
      // ... create it.
      fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    }

    // close the property list
    H5Pclose(plist_id);

    // make sure the file is open
    assert(fileID>=0);
}

void HDF5File::Close() {

    if( fileID<0 ) { // if the file is already closed ...
      // ... do nothing
      return;
    }

    // flush the file
    FlushFile();

    // close the file
    H5Fclose(fileID);

    // set the file ID to something invalid
    fileID = -1;
    filename = "";
}

void HDF5File::Copy(std::string const& dstName, std::shared_ptr<HDF5File> srcFile, std::string const& srcName)
{

  // make sure both files are open
  assert(fileID>0);
  assert(srcFile->fileID>0);

  herr_t err;
  err = H5Ocopy(srcFile->fileID, srcName.c_str(), fileID, dstName.c_str(), H5P_DEFAULT, H5P_DEFAULT);

  if(err<0)
  {
    std::cerr << "WARNING: HDF5 could not copy " << srcName << " to " << dstName << std::endl;
  };

}

bool HDF5File::DoesGroupExist(std::string const& name) const {

  // if the group is the root, return true
  if( (name.compare("/")==0) || (name.compare("")==0) || (name.compare("/.")==0) ) {
    return true;
  }

  // make sure the file is open
  assert(fileID>0);

  // get the group path and the path to it's parent
  std::string parentPath = GetParentPath(name);

  // recursivly check if the parent exists and make sure the current group exists
  return DoesGroupExist(parentPath) && (H5Lexists(fileID, name.c_str(), H5P_DEFAULT)>0);
}

bool HDF5File::DoesDataSetExist(std::string const& name) const {

  // make sure the file is open
  assert(fileID>0);

  // get the group path and the path to it's parent
  std::string parentPath = GetParentPath(name);

  // recursivly check if the parent group exists and make sure the current data set exists
  return DoesGroupExist(parentPath) && (H5Lexists(fileID, name.c_str(), H5P_DEFAULT) > 0);
}

Eigen::VectorXi HDF5File::GetDataSetSize(std::string const name) const {

  // make sure the file is open
  assert(fileID>0);

  if( !DoesDataSetExist(name) ) { // if the data set does not exist ...
    // return an empty vector.
    return Eigen::VectorXi();
  }

  // make sure the file is open
  assert(fileID>0);

  // open the data
  hid_t dataset = H5Dopen2(fileID, name.c_str(), H5P_DEFAULT);

  // get the id for the dataspace of the dataset
  hid_t space_id = H5Dget_space(dataset);

  // get the dimensionality of the dataspace
  int rank = H5Sget_simple_extent_ndims(space_id);

  // get the dataspace dimension size and the max. size
  hsize_t* dims = (hsize_t*)malloc(rank*sizeof(hsize_t));
  hsize_t* max_dims = (hsize_t*)malloc(rank*sizeof(hsize_t));
  H5Sget_simple_extent_dims(space_id, dims, max_dims);

  // close the dataspace and the dataset
  H5Sclose(space_id);
  H5Dclose(dataset);

  // convert the dimensionality into an Eigen::VectorXi
  Eigen::VectorXi output(rank);
  for( int i=0; i<rank; ++i ) {
    output(i) = dims[i];
  }

  // free the memory
  free(dims);
  free(max_dims);

  // return the dimensionality
  return output;
}

void HDF5File::CreateGroup(std::string const& name) {

    // make sure the file is open
    assert(fileID>0);

    if( (DoesGroupExist(name))||(name.compare("")==0)||(name.compare("/")==0) ) { return; }

    // get the group path and the path to it's parent
    std::string parentPath = GetParentPath(name);

    // make sure the parent exists by recursively creating it
    if(!DoesGroupExist(parentPath))
	CreateGroup(parentPath);

    // create the group
    hid_t newgroup = H5Gcreate2(fileID, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // close the group
    H5Gclose(newgroup);
}

void HDF5File::WriteStringAttribute(std::string const& datasetName,
				    std::string const& attributeName,
				    std::string const& attribute)
{
    // make sure the file is open
    assert(fileID>0);

    // Create the group or dataset in necessary
    if( !DoesDataSetExist(datasetName) || !DoesGroupExist(datasetName) )
	CreateGroup(datasetName);

    // write the attribute
    H5LTset_attribute_string(fileID, datasetName.c_str(), attributeName.c_str(), attribute.c_str());
}

std::string HDF5File::GetStringAttribute(std::string const& datasetName, std::string const& attributeName) const {
  /*#if MUQ_MPI==1
  std::unique_ptr<mpi::communicator> worldComm(new mpi::communicator);

  assert(worldComm->rank()==write);
  #endif*/

  // make sure the file is open
  assert(fileID>0);

  // make sure the dataset exists
  assert(DoesDataSetExist(datasetName) || DoesGroupExist(datasetName));

  // get the string attribute
  char tempStr[256];
  H5LTget_attribute_string(fileID, datasetName.c_str(), attributeName.c_str(), tempStr);

  // return it as a strng
  return std::string(tempStr);
}

void HDF5File::FlushFile() {
  if( fileID>0 ) { // if the file is open ...
    // flush it.
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
  }
}

struct DataFileInfo {
    DataFileInfo(std::shared_ptr<HDF5File> const& hdf5file) : hdf5file(hdf5file) {}

    const std::shared_ptr<HDF5File> hdf5file;
};

herr_t CopyObjectToGlobalFile(hid_t o_id, const char *name, const H5O_info_t *info, void *op_data) {
    std::string nameBuffer(name);
    std::string fullGroupName = "/" + nameBuffer;

    // get the file we are copying into
    DataFileInfo* fileInfo = static_cast<DataFileInfo*>(op_data);

    if( info->type==H5O_TYPE_DATASET ) {  // data sets
	if( !fileInfo->hdf5file->DoesDataSetExist(fullGroupName) ) { // if the data set does not exist ...
	    // ... copy it over
	    H5Ocopy(o_id, name, fileInfo->hdf5file->fileID, fullGroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT);
	}
    } else if( info->type == H5O_TYPE_GROUP ) { // groups
	if( !fileInfo->hdf5file->DoesGroupExist(fullGroupName) ) { // if the group does not exist ...
	    // ... copy it over.
	    H5Ocopy(o_id, name, fileInfo->hdf5file->fileID, fullGroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT);
	}
    }

    return 0;
}

void HDF5File::MergeFile(std::shared_ptr<HDF5File> const& otherFile) {

    // make sure the other file is open
    assert(otherFile->fileID>0);

    // make sure this file is open
    assert(fileID>0);

    // open the root group in the other file
    const std::string rootGroupName = "/";
    const hid_t otherRootGroup = H5Gopen2(otherFile->fileID, rootGroupName.c_str(), H5P_DEFAULT);

    auto dataInfo = std::make_shared<DataFileInfo>(shared_from_this());

    // copy the file
    const herr_t status = H5Ovisit(otherRootGroup, H5_INDEX_NAME, H5_ITER_NATIVE, &CopyObjectToGlobalFile, static_cast<void*>(dataInfo.get()));

    assert(status >= 0);

    // close the other file's root group
    H5Gclose(otherRootGroup);
}


bool HDF5File::IsDataSet(std::string const& name) const
{

    if(!DoesDataSetExist(name))
	return false;

    herr_t status;
    H5O_info_t info;

    status = H5Oget_info_by_name(fileID, name.c_str(), &info, H5P_DEFAULT);

    if(status<0)
	return false;

    return info.type == H5O_TYPE_DATASET;
}

bool HDF5File::IsGroup(std::string const& name) const
{

    if(!DoesGroupExist(name))
	return false;

    herr_t status;
    H5O_info_t info;

    status = H5Oget_info_by_name(fileID, name.c_str(), &info, H5P_DEFAULT);

    if(status<0)
	return false;

    return info.type == H5O_TYPE_GROUP;
}

std::vector<std::string> HDF5File::GetChildren(std::string base) const
{
    // Make sure the HDF5 file is open
    assert(fileID>0);

    if(IsDataSet(base))
	return std::vector<std::string>();

    // Make sure the group exists
    assert(DoesGroupExist(base));

    // open the group
    hid_t gid = H5Gopen2(fileID, base.c_str(), H5P_DEFAULT);

    char name[1024];
    ssize_t len;
    hsize_t nobj;

    // get the number of objects in this group
    herr_t status = H5Gget_num_objs(gid, &nobj);

    // Intialize the vector of strings
    std::vector<std::string> output(nobj);

    // Fill in the output vector
    for(int i = 0; i < nobj; i++)
    {
	len = H5Gget_objname_by_idx(gid, (hsize_t)i, name, (size_t)1024);
	output.at(i) = std::string(name,name + len);
    }

    return output;

};
