set(STANMATH_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/stanmath/src/STANMATH")


include(ExternalProject)
if(NOT STANMATH_EXTERNAL_SOURCE)
	set(STANMATH_EXTERNAL_SOURCE https://github.com/stan-dev/math/archive/release/v2.18.0.zip)
endif()

set(STANMATH_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external/)

ExternalProject_Add(
 STANMATH
 PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/stanmath
 URL ${STANMATH_EXTERNAL_SOURCE}
 CONFIGURE_COMMAND ""
 BUILD_COMMAND ""
 INSTALL_COMMAND cp -r ${STANMATH_BUILD_DIR}/stan/ ${STANMATH_INSTALL_DIR}include/stan/
)




set_property( TARGET STANMATH PROPERTY FOLDER "Externals")

set(STANMATH_INCLUDE_DIRS ${STANMATH_INSTALL_DIR}/include/)
