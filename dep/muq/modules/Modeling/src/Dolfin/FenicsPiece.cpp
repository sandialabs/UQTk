#include "MUQ/Modeling/Dolfin/FenicsPiece.h"

using namespace muq::Modeling;

/**
FenicsPiece::FenicsPiece(pybind11::object const& problemIn, pybind11::object const& outputFieldIn, std::vector<pybind11::object> const& inputs)
{
    ExtractInputs(inputs);
    
    // First, make sure the input object is actually an instance of LinearVariationalProblem
    CheckProblemType(problemIn, "LinearVariationalProblem");
    CheckProblemType(outputFieldIn, "Function");
    
    problem = SwigExtract(problemIn);
    outputField = SwigExtract(outputFieldIn);
};
*/
