# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
	APPROXIMATION_REGRESSION
	"Polynomial regression (global and local)"
	muqApproximation
	"UTILITIES_MULTIINDEX;MODELING_FLANN;MODELING_CORE;APPROXIMATION_POLYNOMIALS;OPTIMIZATION_CORE"
	"NLOPT"
	""
	Regression/Regression.cpp
	Regression/LocalRegression.cpp
)

CreateCompileGroup(
	APPROXIMATION_POLYNOMIALS
	"Polynomials and multiindices"
	muqApproximation
	"MODELING_CORE;UTILITIES_MULTIINDEX;APPROXIMATION_QUADRATURE"
	""
	""
	Polynomials/IndexedScalarBasis.cpp
	Polynomials/OrthogonalPolynomial.cpp
	Polynomials/Monomial.cpp
	Polynomials/PhysicistHermite.cpp
	Polynomials/ProbabilistHermite.cpp
	Polynomials/Legendre.cpp
	Polynomials/Laguerre.cpp
	Polynomials/Jacobi.cpp
	Polynomials/HermiteFunction.cpp
  Polynomials/BasisExpansion.cpp
  Polynomials/MonotoneExpansion.cpp
)

CreateCompileGroup(
	APPROXIMATION_QUADRATURE
	"Quadrature rules"
	muqApproximation
	"APPROXIMATION_POLYNOMIALS;UTILITIES_MULTIINDEX"
	""
	""
	Quadrature/Quadrature.cpp
  Quadrature/GaussQuadrature.cpp
	Quadrature/FullTensorQuadrature.cpp
	Quadrature/SmolyakQuadrature.cpp
	Quadrature/ClenshawCurtisQuadrature.cpp
	Quadrature/GaussPattersonQuadrature.cpp
	Quadrature/ExponentialGrowthQuadrature.cpp
)

CreateCompileGroup(
	APPROXIMATION_POLYNOMIALCHAOS
	"Polynomial chaos tools"
	muqApproximation
	"MODELING_CORE;APPROXIMATION_QUADRATURE;APPROXIMATION_POLYNOMIALS"
	""
	""
	PolynomialChaos/PolynomialChaosExpansion.cpp
	PolynomialChaos/PCEFactory.cpp
	PolynomialChaos/SmolyakEstimator.cpp
	PolynomialChaos/AdaptiveSmolyakPCE.cpp
	Quadrature/AdaptiveSmolyakQuadrature.cpp
)

CreateCompileGroup(
    APPROXIMATION_GP_Kernels
    "Gaussian process covariance kernels."
    muqApproximation
    "MODELING_CORE"
    "EIGEN3;BOOST;STANMATH"
    ""
    GaussianProcesses/KernelBase.cpp
    GaussianProcesses/MaternKernel.cpp
    GaussianProcesses/ConcatenateKernel.cpp
		GaussianProcesses/ProductKernel.cpp
    GaussianProcesses/SumKernel.cpp
    GaussianProcesses/PeriodicKernel.cpp
    GaussianProcesses/KarhunenLoeveExpansion.cpp
    #GaussianProcesses/SeparableKarhunenLoeve.cpp
    GaussianProcesses/ObservationInformation.cpp

)

CreateCompileGroup(
    APPROXIMATION_GP
    "Tools for constructing Gaussian processes."
    muqApproximation
    "MODELING_CORE;APPROXIMATION_GP_Kernels;INFERENCE_FILTERING"
    "EIGEN3;BOOST"
    "NLOPT"

    GaussianProcesses/StateSpaceGP.cpp
    GaussianProcesses/GaussianProcess.cpp
)
