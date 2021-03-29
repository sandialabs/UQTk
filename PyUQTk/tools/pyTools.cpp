//=====================================================================================
//
//                      The UQ Toolkit (UQTk) version @UQTKVERSION@
//                          Copyright (@UQTKYEAR@) NTESS
//                        https://www.sandia.gov/UQToolkit/
//                        https://github.com/sandialabs/UQTk
//
//     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

PYBIND11_MODULE(tools, m) {
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

  m.def("invRos",static_cast<void (*)(Array1D<double>&, Array2D<double>&, Array1D<double>&, Array1D<double>&)>(&invRos));
  m.def("invRos",static_cast<void (*)(Array1D<double>&, Array2D<double>&, Array1D<double>&, double)>(&invRos));
  m.def("invRos",static_cast<void (*)(Array1D<double>&, Array2D<double>&, Array1D<double>&)>(&invRos));
  m.def("invRos",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&)>(&invRos));
  m.def("get_opt_KDEbdwth",&get_opt_KDEbdwth);
  m.def("Rosen",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, Array1D<double>&)>(&Rosen));
  m.def("Rosen",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, double)>(&Rosen));
  m.def("Rosen",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&)>(&Rosen));

}
