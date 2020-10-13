#include "MUQ/Utilities/HDF5/HDF5Types.h"
#include <iostream>

using namespace muq::Utilities;

hid_t HDF5_Type<bool>::GetFlag()
{
  if (sizeof(bool) == sizeof(int)) {
    return H5T_NATIVE_INT;
  } else if (sizeof(bool) == sizeof(short)) {
    return H5T_NATIVE_SHORT;
  } else if (sizeof(bool) == sizeof(char)) {
    return H5T_NATIVE_CHAR;
  } else {
    return H5T_NATIVE_CHAR; //safe default, don't read more than one byte
  }
}

herr_t muq::Utilities::H5LTset_attribute(hid_t       loc_id,
                                         const char *obj_name,
                                         const char *attr_name,
                                         hid_t       mem_type_id,
                                         void       *buffer,
                                         size_t      size)
{
  if (mem_type_id == H5T_NATIVE_DOUBLE) {
    return H5LTset_attribute_double(loc_id, obj_name, attr_name, (double *)buffer, size);
  } else if (mem_type_id == H5T_NATIVE_FLOAT) {
    return H5LTset_attribute_float(loc_id, obj_name, attr_name, (float *)buffer, size);
  } else if (mem_type_id == H5T_NATIVE_INT) {
    return H5LTset_attribute_int(loc_id, obj_name, attr_name, (int *)buffer, size);
  } else if (mem_type_id == H5T_NATIVE_UINT) {
    return H5LTset_attribute_uint(loc_id, obj_name, attr_name, (unsigned *)buffer, size);
  } else if (mem_type_id == H5T_NATIVE_SHORT) {
    return H5LTset_attribute_short(loc_id, obj_name, attr_name, (short *)buffer, size);
  } else if (mem_type_id == H5T_NATIVE_CHAR) {
    return H5LTset_attribute_char(loc_id, obj_name, attr_name, (char *)buffer, size);
  } else {
    return -1;
  }
  return -1;
}
