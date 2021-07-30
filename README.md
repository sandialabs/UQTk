# Sandia National Labs
# Uncertainty Quantification Toolkit (UQTk) version 3.1.1

#### Bert Debusschere, Cosmin Safta, Katherine Johnston, Kenny Chowdhary, Khachik Sargsyan, Luke Boll, Mohammad Khalil, Prashant Rai, Tiernan Casey, Xiaoshu Zeng

## Overview
The UQ Toolkit (UQTk) is a collection of libraries and tools for the
quantification of uncertainty in numerical model predictions. Version
3.1.1 offers Polynomial Chaos Expansions to represent random variables,
intrusive and non-intrusive methods for propagating uncertainties through
computational models, tools for sensitivity analysis, methods for sparse
surrogate construction, and Bayesian inference tools for inferring parameters
and model uncertainties from experimental data.

## Documentation
For documentation on how to install and use UQTk, please refer to the manual,
which is included as a PDF file in the directory doc/UQTk_manual.pdf. For more
detailed documentation on the actual source code for development purposes,
please refer to the doxygen documentation in the directory doc/doxy or
[online at https://www.sandia.gov/UQToolkit/doc/](https://www.sandia.gov/UQToolkit/doc/). If you are
familiar with UQTk and would like just a high level overview of where to find
everything and how to install, see the sections on Directory Structure and
Installation below.

## Directory Structure
On a high level UQTk is organized as follows:
* config: Example CMake configuration scripts
* cpp/lib: Core C++ libraries
* cpp/app: Standalone apps that make UQTk functionality available to the command line
* cpp/tests: CMake Unit Tests
* dep: Third party libraries that UQTk depends on
* examples: short tutorial style examples that illustrate key UQTk capabilities
* PyUQTk: Python wrappers for the core C++ libraries as well as additional Python tools

In many key directories, README files have been included to further lay out the
contents of their subdirectories.

## Capabilities
Below is a list of key UQTk capabilities, along with examples that illustrate those
capabilities:
* Intrusive Forward UQ: examples/ops (C++), examples/surf_rxn/SurfRxnISP.cpp (C++)
* Non-Intrusive Forward UQ: examples/surf_rxn/SurfRxnNISP.cpp (C++), examples/fwd_prop (Python),
  examples/window (Python), examples/uqpc (Command Line/Python)
* Non-Intrusive Surrogate Construction: examples/uqpc (Command Line/Python)
* Bayesian Compressive Sensing (BCS): examples/pce_bcs (C++)
* Global Sensitivity Analysis: examples/uqpc (Command Line/Python), examples/pce_bcs (C++), examples/sensMC (Command Line)
* Bayesian Inference: examples/line_infer (C++), examples/iuq (Command Line/Python), examples/polynomial (Python)
* Bayesian model selection: examples/polynomial (Python)
* Transitional Markov chain Monte Carlo (TMCMC): examples/tmcmc_bimodal (C++/Python/Command Line)
* Karhunen-Loève decompositions: examples/kle_ex1 (C++)
* Data Free Inference (Inference based on summary statistics): examples/dfi (C++)
* Forward Propagation with Basis Adaptation: examples/d_spring_series (Python)
* Numerical Integration (Quadrature): examples/num_integ (Python)

For more details on these capabilities, please refer to the UQTk manual in PDF format.

## Installation
To install UQTk, first create a build directory outside of the UQTk repository.
From within the build directory, configure the distribution via CMake. See example
CMake configuration scripts in the directory config
Then build via ``make``, and test with ``ctest``. Install with ``make install``
For example:
```
% mkdir build
% cd build
% ../UQTk/config/config-gcc-Python.sh
% make -j 8
% ctest
% make install
```

For more details, please refer to the UQTk manual in PDF format.


## How to Cite
To cite UQTk, please use the following publications:

```
@ARTICLE{DebusscherePCE:2004,
  author   =  {B.J. Debusschere and H.N. Najm and P.P. P\'ebay and O.M. Knio
               and R.G. Ghanem and O.P. {Le Ma{\^\i}tre}},
  title    =  {Numerical challenges in the use of polynomial chaos representations
               for stochastic processes},
  journal  =  {{SIAM} Journal on Scientific Computing},
  year     =  {2004},
  volume   =  {26},
  pages    =  {698-719},
  number   =  {2}
  url      =  {http://dx.doi.org/10.1137/S1064827503427741}
}

@InCollection{DebusschereUQTk:2017,
  author    = {B. Debusschere and K. Sargsyan and C. Safta and K. Chowdhary},
  title     = {The Uncertainty Quantification Toolkit (UQTk)},
  booktitle = {Handbook of Uncertainty Quantification},
  editor    = {R. Ghanem and D. Higdon and H. Owhadi},
  year      = {2017},
  pages     = {1807--1827},
  publisher = {Springer}
  url       = {http://www.springer.com/us/book/9783319123844}
}
```

## Contact Us
For more information, visit the [UQTk website at https://www.sandia.gov/UQToolkit/](https://www.sandia.gov/UQToolkit/) or
contact the [UQTk Developers at <uqtk-developers@software.sandia.gov>](mailto:uqtk-developers@software.sandia.gov)
Sandia National Laboratories, Livermore, CA, USA.
