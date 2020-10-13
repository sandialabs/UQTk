# make sure that the HDF5HL library is available
set(CMAKE_REQUIRED_LIBRARIES ${HDF5HL_LIBRARIES} ${HDF5_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${HDF5HL_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <hdf5.h>
  #include <hdf5_hl.h>
  #include <stdlib.h>

  #define ATTR_SIZE 5

  int main( void )
  {
   hid_t   file_id;
   hid_t   dset_id;
   hid_t   space_id;
   hsize_t dims[1] = { ATTR_SIZE };
   int     data[ATTR_SIZE] = {1,2,3,4,5};
   herr_t  status;
   int     i;

   char str[128];
   file_id = H5Fcreate(str, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   space_id = H5Screate_simple(1, dims, NULL);

   dset_id = H5Dcreate2(file_id, str, H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dclose(dset_id);
   status = H5Sclose(space_id);

   status = H5LTset_attribute_int(file_id, str, str, data, ATTR_SIZE);
   status = H5LTget_attribute_int(file_id, str, str, data);

   status = H5Fclose(file_id);

   return 0;
  }
  
  "
  HDF5HL_COMPILES)

	
	if(NOT HDF5HL_COMPILES)
		set(HDF5HL_TEST_FAIL 1)
	else()
		set(HDF5HL_TEST_FAIL 0)
	endif()