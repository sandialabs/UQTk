# make sure that the HDF5 library is available
set(CMAKE_REQUIRED_LIBRARIES ${HDF5_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <hdf5.h>
  #define NX     5                      /* dataset dimensions */
  #define NY     6
  #define RANK   2

  int
  main (void)
  {
      hid_t       file, dataset;         /* file and dataset handles */
      hid_t       datatype, dataspace;   /* handles */
      hsize_t     dimsf[2];              /* dataset dimensions */
      herr_t      status;                             
      int         data[NX][NY];          /* data to write */
      int         i, j;

	   char NAME[128];
      file = H5Fcreate(NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
	  dimsf[0] = NX;
      dimsf[1] = NY;
      dataspace = H5Screate_simple(RANK, dimsf, NULL); 
      datatype = H5Tcopy(H5T_NATIVE_INT);
	  
      status = H5Tset_order(datatype, H5T_ORDER_LE);
      //dataset = H5Dcreate(file, NAME, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

      H5Sclose(dataspace);
      H5Tclose(datatype);
      H5Dclose(dataset);
      H5Fclose(file);
 
      return 0;
  }     

  
  "
  HDF5_COMPILES)
  
if(MUQ_USE_OPENMPI)
  
  set(CMAKE_REQUIRED_LIBRARIES ${HDF5_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIR})
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <hdf5.h>
    #define NX     5                      /* dataset dimensions */
    #define NY     6
    #define RANK   2

    int
    main (void)
    {
      MPI_Comm comm     = MPI_COMM_SELF;
      MPI_Info info     = MPI_INFO_NULL;
      hid_t    plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, comm, info);
      return 0;
    }     
    "
    HDF5_MPI_COMPILES)
else()
  set(HDF5_MPI_COMPILES 1)
endif()
  
if(NOT HDF5_COMPILES OR NOT HDF5_MPI_COMPILES)
        set(HDF5_TEST_FAIL 1)
else()
        set(HDF5_TEST_FAIL 0)
endif()

