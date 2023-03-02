/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#ifndef PCMAPS_H
#define PCMAPS_H

/** \file pcmaps.h
 * \brief Header for suite of functions to help map one kind of a PC variable to another.
 * \todo Perhaps use more robust tools, like dcdflib.
 * \todo Need more testing of these tools.
 */

/// \brief Implements a map y=f(x), where f is a function mapping one PC domain (pcIn with parameters in1,in2) 
/// to another (pcOut with parameters out1,out2) 
/// \note Besides standard PC types, it also incorporates 
/// a) 'TG' : truncated-gaussian variable, 
/// b) 'RB' : Roe-Baker PDF from Roe, G. H., & Baker, M. B. (2007). Why is climate sensitivity so unpredictable?. Science, 318(5850), 629-632.
/// c) 'pdf': Given arbitrary cumulative distribution function values in cdf.dat it maps corresponding r.v. to PC variables as well
/// \param[in] x     : Input scalar x
/// \param[in] pcIn  : PC type for input x (options are LU, HG, LG, JB, SW, TG, RB, pdf)
/// \param[in] in1   : Parameter #1 for input PC (if relevant, i.e. for LG, JB, SW, TG, RB) 
/// \param[in] in2   : Parameter #2 for input PC (if relevant, i.e. for JB, SW, RB) 
/// \param[in] pcOut : PC type for output y (options are LU, HG, LG, JB, SW, TG, RB, pdf)
/// \param[in] out1  : Parameter #1 for output PC (if relevant, i.e. for LG, JB, SW, TG, RB) 
/// \param[in] out2  : Parameter #2 for output PC (if relevant, i.e. for JB, SW, RB) 
/// \return    y     : Output scalar y
/// \note The user is free to choose any value of x that is in the PC domain defined by pcIn
/// \note The map f() is a map frpm pcIn germ to a pcOut germ
/// \note For example, y=PCtoPC(x,'HG',0,0,'LU',0,0) maps \f$(-\infty,\infty)\f$ to \f$[-1,1]\f$ as a map from standard normal r.v. to uniform r.v.
double PCtoPC(double x, const std::string pcIn, double in1, double in2, const std::string pcOut, double out1, double out2);

/// \brief Implements PCtoPC() map entrywise from array xx to yy
void PCtoPC(Array2D<double>& xx, const std::string pcIn, double in1, double in2, Array2D<double>& yy, const std::string pcOut, double out1, double out2);

/// \brief Bisection method for root-finding, modified to invert PCtoPC maps
double rtbis_mod(double func(double,const std::string,double,double,const std::string,double,double), const double x1, const double x2, const double xacc,double x, const std::string pcIn, double in1, double in2, const std::string pcOut, double out1, double out2);

/// \brief Auxiliary linear interpolation function
void linint( Array2D<double> &xydata, const double x, double &y, int col ) ;

/// \brief Auxiliary linear interpolation function
/// \note Currently not used, as there is an overloaded linint() function that is more general
void linint( Array2D<double> &xydata, const double x, double &y ) ;


//---------------------------------------------------------------------------------------
#endif // PCMAPS_H
