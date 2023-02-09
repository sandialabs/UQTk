/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

// Include std header files
#include <math.h>
#include <algorithm>

// Include UQTk header files
#include "PCSet.h"
#include "model.h"

double eval_surrogate_model(Array2D<double> &parameters, int d, int n, Info &info)
{
    // Evaluate the PCE
    Array1D<double> out(1);
    info.pces(d)(n)->EvalPCAtCustPoints(out, parameters, info.pccfs(d)(n));

    // Return value
    return out(0);
}

void eval_grad_surrogate_model(Array2D<double> &grad, Array2D<double> &parameters, int d, int n, Info &info)
{
    info.pces(d)(n)->dPhi(parameters, info.mindices(d)(n), grad, info.pccfs(d)(n));
}