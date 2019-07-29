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
/** \file gen_defs.h
 * Miscellaneous definitions used in several utilities. 
 */

#include "assert.h"

#ifndef GENDEFS_H
#define GENDEFS_H

#define BOOLEAN unsigned char

#ifndef TRUE
#  define TRUE  ((BOOLEAN)1)
#  define FALSE ((BOOLEAN)0)
#endif

#define MAX(a,b)        (((a) > (b)) ? (a) : (b))
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))
#define SIGN(a,b) ( (b) >= 0.0 ? fabs(a) : -fabs(a) )
#define SWAP(a,b) { double temp=(a); (a)=(b); (b)=temp; }
#define CHECKEQ(a,b) { if ( (a) != (b) ) { cout << "Equality check failed: " << (a) << " not equal to " << (b) << endl; assert( (a)==(b) );} }

extern char *optarg;




#endif // GENDEFS_H
