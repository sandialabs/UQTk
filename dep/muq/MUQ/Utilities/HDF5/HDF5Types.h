#ifndef HDF5TYPES_H_
#define HDF5TYPES_H_

#include <iostream>

#include <hdf5.h>
#include <hdf5_hl.h>

namespace muq
{

namespace Utilities
{
    herr_t H5LTset_attribute(hid_t loc_id, const char *obj_name, const char *attr_name, hid_t mem_type_id, void *buffer, size_t size);

#if MUQ_MPI==1
    template<typename scalarType>
      struct MPI_Type {
	static MPI_Datatype GetFlag() { return MPI_DATATYPE_NULL; }
      };

    template<>
      struct MPI_Type<double> {
      static MPI_Datatype GetFlag() { return MPI_DOUBLE; }
    };

    template<>
      struct MPI_Type<char> {
      static MPI_Datatype GetFlag() { return MPI_CHAR; }
    };

    template<>
      struct MPI_Type<short> {
      static MPI_Datatype GetFlag() { return MPI_SHORT; }
    };

    template<>
      struct MPI_Type<int> {
      static MPI_Datatype GetFlag() { return MPI_INT; }
    };

    template<>
      struct MPI_Type<long> {
      static MPI_Datatype GetFlag() { return MPI_LONG; }
    };

    template<>
      struct MPI_Type<float> {
      static MPI_Datatype GetFlag() { return MPI_FLOAT; }
    };

    template<>
      struct MPI_Type<unsigned short> {
      static MPI_Datatype GetFlag() { return MPI_UNSIGNED_SHORT; }
    };

    template<>
      struct MPI_Type<unsigned int> {
      static MPI_Datatype GetFlag() { return MPI_UNSIGNED; }
    };

    template<>
      struct MPI_Type<unsigned long> {
      static MPI_Datatype GetFlag() { return MPI_UNSIGNED_LONG; }
    };
#endif

    template<typename scalarType>
      struct HDF5_Type{
	static hid_t GetFlag() { return -1; }
      };

    template<>
    struct HDF5_Type<double>{
      static hid_t GetFlag(){return H5T_NATIVE_DOUBLE;};
    };
    template<>
    struct HDF5_Type<long double>{
      static hid_t GetFlag(){return H5T_NATIVE_LDOUBLE;};
    };
    template<>
    struct HDF5_Type<int>{
      static hid_t GetFlag(){return H5T_NATIVE_INT;};
    };
    template<>
    struct HDF5_Type<long>{
      static hid_t GetFlag(){return H5T_NATIVE_LONG;};
    };
    template<>
    struct HDF5_Type<unsigned long>{
      static hid_t GetFlag(){return H5T_NATIVE_ULONG;};
    };
    template<>
    struct HDF5_Type<unsigned>{
      static hid_t GetFlag(){return H5T_NATIVE_UINT;};
    };

    template<>
    struct HDF5_Type<float>{
      static hid_t GetFlag(){return H5T_NATIVE_FLOAT;};
    };

    template<>
    struct HDF5_Type<unsigned short>{
      static hid_t GetFlag(){return H5T_NATIVE_USHORT;};
    };
    template<>
    struct HDF5_Type<short>{
      static hid_t GetFlag(){return H5T_NATIVE_SHORT;};
    };

    template<>
    struct HDF5_Type<char>{
      static hid_t GetFlag(){return H5T_NATIVE_CHAR;};
    };

    template<>
    struct HDF5_Type<bool>{

      static hid_t GetFlag();
    };

} // namespace Utilities
} // namespace muq

#endif // HDF5TYPES_H_
