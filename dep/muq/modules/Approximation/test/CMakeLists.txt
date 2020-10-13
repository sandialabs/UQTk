CreateTestGroup(
  APPROXIMATION_GP_Kernels

  GaussianProcesses/KernelTests.cpp
  #GaussianProcesses/KarhunenLoeveTest.cpp
)

CreateTestGroup(
  APPROXIMATION_GP

  GaussianProcesses/GaussianProcessTests.cpp
  GaussianProcesses/StateSpaceTests.cpp
  GaussianProcesses/KarhunenLoeveTest.cpp
)

CreateTestGroup(
  APPROXIMATION_POLYNOMIALS

  Polynomials/BasisExpansionTest.cpp
  Polynomials/MonotoneExpansionTest.cpp
  Polynomials/PolynomialTests.cpp
)

CreateTestGroup(
  APPROXIMATION_POLYNOMIALCHAOS

  PolynomialChaos/PCEFactoryTest.cpp
  PolynomialChaos/PolynomialChaosExpansionTest.cpp
  PolynomialChaos/AdaptiveSmolyakPCETest.cpp
  Quadrature/AdaptiveSmolyakQuadratureTest.cpp
)

CreateTestGroup(
  APPROXIMATION_QUADRATURE

  Quadrature/GaussQuadratureTests.cpp
  Quadrature/FullTensorQuadratureTest.cpp
  Quadrature/ClenshawCurtisTests.cpp
)

CreateTestGroup(
	APPROXIMATION_REGRESSION

	RegressionTests/RegressionTests.cpp
	RegressionTests/LocalRegressionTests.cpp
)

CreateParallelTestGroup(
	APPROXIMATION_REGRESSION

	RegressionTests/LocalRegressionParallelTests.cpp
)
