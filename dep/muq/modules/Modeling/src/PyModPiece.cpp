#include "MUQ/Modeling/PyModPiece.h"

using namespace muq::Modeling;

PyModPiece::PyModPiece(Eigen::VectorXi const& inputSizes,
                       Eigen::VectorXi const& outputSizes)
                       : ModPiece(inputSizes, outputSizes) {}

void PyModPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  EvaluateImpl(ToStdVec(input));
}

void PyModPiece::GradientImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input,
                              Eigen::VectorXd             const& sensitivity) {
  GradientImpl(outputDimWrt, inputDimWrt, ToStdVec(input), sensitivity);
}

void PyModPiece::GradientImpl(unsigned int                 const  outputDimWrt,
                          unsigned int                 const  inputDimWrt,
                          std::vector<Eigen::VectorXd> const& input,
                          Eigen::VectorXd              const& sensitivity)
{
    gradient = GradientByFD(outputDimWrt, inputDimWrt, input,sensitivity);
}

void PyModPiece::JacobianImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input) {
  JacobianImpl(outputDimWrt, inputDimWrt, ToStdVec(input));
}

void PyModPiece::JacobianImpl(unsigned int                 const  outputDimWrt,
                          unsigned int                 const  inputDimWrt,
                          std::vector<Eigen::VectorXd> const& input)
{
  jacobian = JacobianByFD(outputDimWrt, inputDimWrt, input);
}


void PyModPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& vec) {
  ApplyJacobianImpl(outputDimWrt, inputDimWrt, ToStdVec(input), vec);
}

void PyModPiece::ApplyJacobianImpl(unsigned int                 const  outputDimWrt,
                                   unsigned int                 const  inputDimWrt,
                                   std::vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd              const& vec)
{
  jacobianAction = ApplyJacobianByFD(outputDimWrt, inputDimWrt, input, vec);
}
