//=====================================================================================
//
//                      The UQ Toolkit (UQTk) version 3.1.2
//                          Copyright (2022) NTESS
//                        https://www.sandia.gov/UQToolkit/
//                        https://github.com/sandialabs/UQTk
//
//     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
//     retains certain rights in this software.
//
//     This file is part of The UQ Toolkit (UQTk)
//
//     UQTk is open source software: you can redistribute it and/or modify
//     it under the terms of BSD 3-Clause License
//
//     UQTk is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     BSD 3 Clause License for more details.
//
//     You should have received a copy of the BSD 3 Clause License
//     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
//
//     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
//     Sandia National Laboratories, Livermore, CA, USA
//=====================================================================================
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

#include "Array2D.h"
#include "Array1D.h"
#include "dsfmt_add.h"
#include "KCenterClustering.h"
#include "figtree_internal.h"

#include "probability.h"
#include "combin.h"
#include "minmax.h"
#include "multiindex.h"
#include "pcmaps.h"
#include "gq.h"
#include "rosenblatt.h"
#include "func.h"

namespace py=pybind11;

PYBIND11_MODULE(_tools, m) {
  // Combin.h
  m.def("choose",&choose);
  m.def("factorial",&factorial);
  m.def("logfactorial",&logfactorial);
  m.def("chooseComb",&chooseComb);
  m.def("get_perm",static_cast<void (*)(int, int*,int)>(&get_perm));
  m.def("get_perm",static_cast<void (*)(Array1D<int>&, int)>(&get_perm));
  m.def("gammai",&gammai);
  m.def("beta",&beta);
  m.def("betai",&betai);
  m.def("digama",&digama);
  m.def("clust",&clust);
  m.def("clust_best",&clust_best);
  m.def("findNumCl",&findNumCl);

  // rosenblatt.h
  m.def("invRos",static_cast<void (*)(Array1D<double>&, Array2D<double>&, Array1D<double>&, Array1D<double>&)>(&invRos));
  m.def("invRos",static_cast<void (*)(Array1D<double>&, Array2D<double>&, Array1D<double>&, double)>(&invRos));
  m.def("invRos",static_cast<void (*)(Array1D<double>&, Array2D<double>&, Array1D<double>&)>(&invRos));
  m.def("invRos",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&)>(&invRos));
  m.def("get_opt_KDEbdwth",&get_opt_KDEbdwth);
  m.def("Rosen",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, Array1D<double>&)>(&Rosen));
  m.def("Rosen",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, double)>(&Rosen));
  m.def("Rosen",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&)>(&Rosen));

  // multiindex.h
  m.def("head_ext_",&heap_ext_);
  m.def("computeNPCTerms",&computeNPCTerms);
  m.def("computeMultiIndex",static_cast<int (*)(int,int, Array2D<int> &)>(&computeMultiIndex));
  m.def("computeMultiIndex",static_cast<int (*)(int,int, Array2D<int> &,string)>(&computeMultiIndex));
  m.def("computeMultiIndexT",&computeMultiIndexT);
  m.def("computeMultiIndexTP",&computeMultiIndexTP);
  m.def("computeNPCTermsHDMR",&computeNPCTermsHDMR);
  m.def("computeMultiIndexHDMR",&computeMultiIndexHDMR);
  m.def("decodeMindex",&decodeMindex);
  m.def("upOrder",&upOrder);
  m.def("is_admis",&is_admis);
  m.def("getOrders",&getOrders);
  m.def("get_invmindex",&get_invmindex);
  m.def("get_invmindex_ord",&get_invmindex_ord);

  //probability.h
  m.def("erff",static_cast<double (*)(const double)>(&erff));
  m.def("inverf",&inverf);
  m.def("invnormcdf",&invnormcdf);
  m.def("normcdf",&normcdf);
  m.def("normcdfc",&normcdfc);
  m.def("generate_uniform",static_cast<void (*)(double*,int, int, int)>(&generate_uniform));
  m.def("generate_uniform",static_cast<void (*)(Array2D<double>&,int)>(&generate_uniform));
  m.def("generate_uniform",static_cast<void (*)(double *, int, int, dsfmt_t *)>(&generate_uniform));
  m.def("generate_uniform",static_cast<void (*)(Array2D<double> &, dsfmt_t *)>(&generate_uniform));
  m.def("generate_uniform_lhs",static_cast<void (*)(double*,int, int, int)>(&generate_uniform_lhs));
  m.def("generate_uniform_lhs",static_cast<void (*)(Array2D<double>&,int)>(&generate_uniform_lhs));
  m.def("generate_uniform_lhs",static_cast<void (*)(double *, int, int, dsfmt_t *)>(&generate_uniform_lhs));
  m.def("generate_uniform_lhs",static_cast<void (*)(Array2D<double> &, dsfmt_t *)>(&generate_uniform_lhs));
  m.def("generate_normal",&generate_normal);
  m.def("generate_normal_lhs",&generate_normal_lhs);
  m.def("get_median",&get_median);
  m.def("get_mean",static_cast<double (*)(const Array1D<double>&)>(&get_mean));
  m.def("get_mean",static_cast<double (*)(const Array2D<double>&)>(&get_mean));
  m.def("get_std",&get_std);
  m.def("get_var",&get_var);
  m.def("getMean_Variance",&getMean_Variance);
  m.def("getMean",static_cast<void (*)(Array2D<double>&, Array1D<double>&)>(&getMean));
  m.def("getMean",static_cast<void (*)(Array2D<double>&, Array1D<double>&,char *)>(&getMean));
  m.def("rperm",&rperm);
  m.def("getPdf_figtree",&getPdf_figtree);
  m.def("getPdf_cl",&getPdf_cl);
  m.def("covariance",&covariance);
  m.def("ihsU",static_cast<void (*)(Array2D<double> &, int, dsfmt_t *)>(&ihsU));
  m.def("ihsU",static_cast<void (*)(int, int, double *, int, dsfmt_t *)>(&ihsU));
  m.def("ihsP",&ihsP);
  m.def("distCorr",&distCorr);

  //pcmaps.h
  m.def("PCtoPC",static_cast<double (*)(double, const std::string, double, double, const std::string, double, double)>(&PCtoPC));
  m.def("PCtoPC",static_cast<void (*)(Array2D<double>&, const std::string, double, double, Array2D<double>&, const std::string, double, double)>(&PCtoPC));
  m.def("rtbis_mod",&rtbis_mod);
  m.def("linint",static_cast<void (*)(Array2D<double> &, const double, double &, int)>(&linint));
  m.def("linint",static_cast<void (*)(Array2D<double> &, const double, double &)>(&linint));

  //func.h
  m.def("Func_Prop",&Func_Prop);
  m.def("Func_PropQuad",&Func_PropQuad);
  m.def("Func_Exp",&Func_Exp);
  m.def("Func_ExpQuad",&Func_ExpQuad);
  m.def("Func_Const",&Func_Const);
  m.def("Func_Linear",&Func_Linear);
  m.def("Func_BB",&Func_BB);
  m.def("Func_HT1",&Func_HT1);
  m.def("Func_HT2",&Func_HT2);
  m.def("Func_FracPower",&Func_FracPower);
  m.def("Func_ExpSketch",&Func_ExpSketch);
  m.def("Func_Inputs",&Func_Inputs);
  m.def("Func_PCl",&Func_PCl);
  m.def("Func_PCx",&Func_PCx);
  m.def("Func_PC",&Func_PC);
  m.def("Func_PCs",&Func_PCs);
  m.def("augment",&augment);

  //minmax.h
  m.def("getDomain",&getDomain);
  m.def("maxVal",static_cast<double (*)(Array1D<double>&)>(&maxVal));
  m.def("maxVal",static_cast<int (*)(const Array1D<int>&)>(&maxVal));
  m.def("maxVal",static_cast<double (*)(const Array2D<double>&)>(&maxVal));
  m.def("maxVal",static_cast<int (*)(const Array2D<int>&)>(&maxVal));
  m.def("minVal",static_cast<double (*)(const Array1D<double>&)>(&minVal));
  m.def("minVal",static_cast<int (*)(const Array1D<int>&)>(&minVal));
  m.def("minVal",static_cast<double (*)(const Array2D<double>&)>(&minVal));
  m.def("minVal",static_cast<int (*)(const Array2D<int>&)>(&minVal));
  m.def("maxIndex",static_cast<int (*)(Array1D<double>&)>(&maxIndex));
  m.def("maxIndex",static_cast<int (*)(Array1D<int>&)>(&maxIndex));
  m.def("minIndex",static_cast<int (*)(Array1D<double>&)>(&minIndex));
  m.def("minIndex",static_cast<int (*)(Array1D<int>&)>(&minIndex));
  m.def("maxIndexC_2D",&maxIndexC_2D);
  m.def("maxIndexC_2D",&minIndexC_2D);

  //gq.h
  m.def("gq",static_cast<void (*)(const int, const double, const double, Array1D<double> &, Array1D<double> &)>(&gq));
  m.def("gq",static_cast<void (*)(const int, const int, const double, const double, double *, double *)>(&gq));
  m.def("gq_gen",&gq_gen);
  m.def("vandermonde_gq",&vandermonde_gq);
  m.def("gchb",&gchb);



}
