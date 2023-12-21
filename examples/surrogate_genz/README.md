## Surrogate Construction Tutorials

The Jupyter notebookes in this folder illustrate several ways to construct
surrogate models for Genz functions

* surrogate_genz-BCS.ipynb:        Constructs a surrogate for a Genz function using BCS, calculating the NRMSE
* surrogate_genz-Galerkin.ipynb:   Constructs a surrogate for a Genz function using Galerkin projection (with full and sparse quadrature), calculating the NRMSE
* surrogate_genz-Regression.ipynb: Constructs a surrogate for a Genz function using regression, calculating the NRMSE
* surrogate_genz.ipynb:            Constructs a surrogate for a Genz function using Galerkin Projection, regression, and BCS, calculating the NRMSE for each

Additionally, the file surrogate_genz-BCS_detailed.ipynb gives more detailed illustrations of some of the BCS functionality available through the Python interface.
