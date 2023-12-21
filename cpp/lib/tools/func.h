/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file func.h
/// \brief Header for implementation of functions of form \f$y=f(\lambda;x)\f$
/// \note Functions of form \f$y=f(\lambda;x)\f$ for 
///       \f$x\in\mathbf{R}^s\f$ and \f$\lambda\in\mathbf{R}^d\f$ at
///       \f$r\f$ values of model parameters \f$\lambda\f$ and \f$n\f$ values of design parameters \f$x\f$
/// \param p Model parameters \f$\lambda\f$ as a matrix \f$r\times d\f$
/// \param x Design parameters \f$x\f$ as a matrix \f$n\times s\f$
/// \param *funcinfo Potentially function-specific information
/// \return y  Output as a matrix \f$r\times n\f$

#ifndef FUNC_H_SEEN
#define FUNC_H_SEEN

#include "Array1D.h"
#include "Array2D.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

using namespace std; // needed for python string conversion



/// \brief \f$y=f(\lambda;x)=\lambda x\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^1\f$
Array2D<double> Func_Prop(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=\lambda_1 x+\lambda_2x^2\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^2\f$
Array2D<double> Func_PropQuad(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=e^{\lambda_1 +\lambda_2x}\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^2\f$
Array2D<double> Func_Exp(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=e^{\lambda_1 +\lambda_2x+\lambda_3x^2}\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^3\f$
Array2D<double> Func_ExpQuad(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=\lambda\f$ for \f$x\in\mathbf{R}^s\f$ and \f$\lambda\in\mathbf{R}^1\f$
Array2D<double> Func_Const(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=\lambda_1+\lambda_2x\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^2\f$
Array2D<double> Func_Linear(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)\f$ a black-box function with a script bb.x which takes p.dat and x.dat and returns output in y.dat
Array2D<double> Func_BB(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief Heat_transfer1: a custom model designed for a tutorial case of a heat conduction problem
/// \brief \f$y=f(\lambda;x)=\frac{x d_w}{A_w \lambda}+T_0\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^1\f$
/// \note hardwired parameters: \f$d_w=0.1, A_w=0.04, T_0=273\f$
Array2D<double> Func_HT1(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief Heat_transfer2: a custom model designed for a tutorial case of a heat conduction problem
/// \brief \f$y=f(\lambda;x)=\frac{x Q}{A_w \lambda_1}+\lambda_2\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^2\f$
/// \note hardwired parameters: \f$A_w=0.04, Q=20.0\f$
Array2D<double> Func_HT2(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=\lambda_1+\lambda_2 x+\lambda_3 x^2+ \lambda_4 (x+1)^{3.5}\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^4\f$
Array2D<double> Func_FracPower(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x)=\lambda_2 e^{\lambda_1 x} - 2\f$ for \f$x\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^2\f$
Array2D<double> Func_ExpSketch(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief \f$y=f(\lambda;x_i)=\lambda_i\f$ for \f$i=1,...,d\f$, \f$x_i\in\mathbf{R}^1\f$ and \f$\lambda\in\mathbf{R}^d\f$
Array2D<double> Func_Inputs(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief Legendre PC expansion with \f$\lambda\f$'s as coefficients
/// \brief \f$y=f(\lambda;x)=\sum_{\alpha\in{\cal S}} \lambda_\alpha \Psi_\alpha(x)\f$ for \f$x\in\mathbf{R}^s\f$ and \f$\lambda\in\mathbf{R}^{|{\cal S}|}\f$
/// \note hardwired parameter: multiindex set \f${\cal S}\f$ is given in a file mindexx.dat
Array2D<double> Func_PCl(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief Legendre PC expansion with respect to \f$z=(\lambda,x)\f$ 
/// \brief \f$y=f(\lambda;x)=\sum_{\alpha\in{\cal S}} c_\alpha \Psi_\alpha(\lambda,x)\f$ for \f$x\in\mathbf{R}^s\f$ and \f$\lambda\in\mathbf{R}^d\f$
/// \note hardwired parameters: multiindex set \f${\cal S}\f$ is given in a file mindexpx.dat, coefficients \f$c_\alpha\f$ given in a file pccfpx.dat
Array2D<double> Func_PCx(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief Legendre PC expansion for each value of \f$x\f$ 
/// \brief \f$y=f(\lambda;x^{(i)})=\sum_{\alpha\in{\cal S}} c_{\alpha,i} \Psi_\alpha(\lambda)\f$ for \f$x\in\mathbf{R}^s\f$ and \f$\lambda\in\mathbf{R}^d\f$
/// \note hardwired parameters: common multiindex set for all PCs \f${\cal S}\f$ is given in a file mindexp.dat, coefficients \f$c_{\alpha,i}\f$ are given in a file pccf_all.dat
Array2D<double> Func_PC(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);

/// \brief Legendre PC expansion for each value of \f$x\f$ 
/// \brief \f$y=f(\lambda;x^{(i)})=\sum_{\alpha\in{\cal S}_i} c_{\alpha,i} \Psi_\alpha(\lambda)\f$ for \f$x\in\mathbf{R}^s\f$ and \f$\lambda\in\mathbf{R}^d\f$
/// \note hardwired parameters: multiindex sets for all PCs \f${\cal S}\f$ are given in files mindexp.i.dat, coefficients \f$c_{\alpha,i}\f$ are given in files pccfp.i.dat
Array2D<double> Func_PCs(Array2D<double>& p, Array2D<double>& x, Array2D<double>& fixindnom, void* funcinfo);


/// \brief Augments a parameter matrix with 'fixed' columns given indices and nominal values of those
Array2D<double> augment(Array2D<double>& p, Array2D<double>& fixindnom);


#endif /* FUNC_H_SEEN */
