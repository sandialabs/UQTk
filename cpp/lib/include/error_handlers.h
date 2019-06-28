/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
#ifndef ERROR_HANDLERS_H_SEEN
#define ERROR_HANDLERS_H_SEEN

#include <string>
#include <iostream>

using namespace std;

/// \class Tantrum
/// \brief Error handler: allows program to display an error message
/// and exit
/// \todo Build in more smarts for tracing error etc.
struct Tantrum {
        /// \brief String literal to store error message in
	const char* p;
        /// \brief Constructor that takes a string literal as error message
	/// to display before core dumping
	Tantrum(const char*q) {p=q; cout << p << endl;}
        /// \brief Constructor that takes a string as error message
	/// to display before core dumping
	Tantrum(std::string q) {p=q.c_str(); cout << p << endl;}
};



#endif /* ERROR_HANDLERS_H_SEEN */
