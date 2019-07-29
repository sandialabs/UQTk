/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */

/** \file ftndefs.h
 * header file with definitions and macros to use when calling fortran
 * functions from C/C++ and vice versa
 * inspired by http://www.thp.univie.ac.at/~jthorn/c2f.html
 */

/* this header file is idempotent */
#ifndef FTNDEFS_H_SEEN
#define FTNDEFS_H_SEEN
/**********************************************************************/

/*
 * C/C++ compatibility:
 */

/*
 * Use this in prototypes like this:  extern FTN_FUNC void foo(...)
 *
 * At present, this is set up to tell a C++ compiler that  foo()  uses
 * a C-compatible calling convention.
 */
#ifdef __cplusplus
  #define FTN_FUNC	"C"
#else
  #define FTN_FUNC	/* empty */
#endif

/*
 * Macro to convert function names to Fortran format
 */

/*
 * Names of Fortran routines are often altered by the compiler/loader.  The
 * following macro should be used to call a Fortran routine from C code, i.e.
 *	call sgefa(...)			-- Fortran code
 *	FTN_NAME(sgefa)(...);		-- C code to do the same thing
 * Make sure the fortran names are always lower case in the C code.
 * Specify "wsu" for "with single underscore" if the Fortran compiler appends a single
 * underscore to function names. Use "wdu" for "with double underscore".
 */

#if defined(__wsu)
  #define FTN_NAME(n_)   n_ ## _
#elif defined(__wdu)
  #define FTN_NAME(n_)   n_ ## __
#else
  #error "FTN_NAME macros not defined for this system"
#endif

/**********************************************************************/
#endif	/* FNTDEFS_H_SEEN */
