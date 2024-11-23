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
