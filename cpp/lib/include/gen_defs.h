/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
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
