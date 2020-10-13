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
    MODELING_CORE
    "The basic building blocks of the Modeling library."
    muqModeling
    "MODELING_LINEARALGEBRA;UTILITIES_CORE"
    "EIGEN3;BOOST"
    ""
    WorkPiece.cpp

    NodeNameFinder.cpp

    WorkGraph.cpp
    WorkGraphNode.cpp
    WorkGraphEdge.cpp

    PyModPiece.cpp

    ModPiece.cpp
    ModGraphPiece.cpp
    WorkGraphPiece.cpp
    LinearSDE.cpp

    ConstantPiece.cpp
    ConstantVector.cpp
    IdentityPiece.cpp
    ProductPiece.cpp
    SplitVector.cpp
    ScaleVector.cpp
    CombineVectors.cpp
    SumPiece.cpp

    MultiLogisticLikelihood.cpp

    ReplicateOperator.cpp

    LinearSDE.cpp

)

CreateCompileGroup(
  MODELING_STAN
  "Modeling tools that depend on STAN."
  muqModeling
  "MODELING_CORE"
  "EIGEN3;BOOST;STANMATH"
  ""
)

CreateCompileGroup(
  MODELING_FLANN
  "Modeling tools that also depend on FLANN."
  muqModeling
  "MODELING_CORE"
  "EIGEN3;BOOST;NANOFLANN"
  ""
  Flann/FlannCache.cpp
)

CreateCompileGroup(
    MODELING_LINEARALGEBRA
    "Extensions and generalizations of Eigen for linear algebra."
    muqModeling
    "MODELING_CORE"
    "EIGEN3"
    ""
    LinearAlgebra/LinearOperator.cpp
    LinearAlgebra/CompanionMatrix.cpp
    LinearAlgebra/BlockDiagonalOperator.cpp
    LinearAlgebra/BlockRowOperator.cpp
    LinearAlgebra/KroneckerProductOperator.cpp
    LinearAlgebra/SumOperator.cpp
    LinearAlgebra/ConcatenateOperator.cpp
    LinearAlgebra/DiagonalOperator.cpp
    LinearAlgebra/ProductOperator.cpp
    LinearAlgebra/AffineOperator.cpp
    
    #LinearAlgebra/AnyAlgebra.cpp
    #LinearAlgebra/ScalarAlgebra.cpp
    #LinearAlgebra/EigenVectorAlgebra.cpp
    #LinearAlgebra/EigenMatrixAlgebra.cpp
    #LinearAlgebra/SundialsAlgebra.cpp
)

CreateCompileGroup(
    MODELING_DISTRIBUTIONS
    "Implement probability distributions, which allow us to evaluate log-densities and/or sample the distirubtion."
    muqModeling
    ""
    "MODELING_CORE;UTILITIES_CORE"
    ""

    Distributions/Distribution.cpp
    Distributions/PyDistribution.cpp
    Distributions/UniformBox.cpp
    Distributions/Gaussian.cpp
    Distributions/GaussianBase.cpp
    Distributions/InverseGamma.cpp
    Distributions/DensityProduct.cpp
    Distributions/Density.cpp
    Distributions/RandomVariable.cpp
)

CreateCompileGroup(
    MODELING_SUNDIALS_MODELS
    "ODE integrators using Sundials"
    muqModeling
    "MODELING_CORE"
    "EIGEN3;BOOST;SUNDIALS"
    ""
    ODEBase.cpp
    ODEData.cpp
    ODE.cpp
    #RootfindingIVP.cpp
)
