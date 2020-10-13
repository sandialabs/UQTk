#include "AllClassWrappers.h"

#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/AffineOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"
//#include "MUQ/Modeling/LinearAlgebra/SparseLinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockRowOperator.h"
#include "MUQ/Modeling/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Modeling/LinearAlgebra/ConcatenateOperator.h"
#include "MUQ/Modeling/LinearAlgebra/DiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/KroneckerProductOperator.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/ProductOperator.h"
#include "MUQ/Modeling/LinearAlgebra/SumOperator.h"
#include "MUQ/Modeling/LinearAlgebra/ZeroOperator.h"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Modeling::PythonBindings;
using namespace muq::Modeling;
namespace py = pybind11;


void muq::Modeling::PythonBindings::LinearOperatorWrapper(py::module &m)
{

  py::class_<LinearOperator, ModPiece, WorkPiece, std::shared_ptr<LinearOperator>> lo(m, "LinearOperator");
  lo
    .def("rows", &LinearOperator::rows)
    .def("cols", &LinearOperator::cols)
    .def("GetMatrix", &LinearOperator::GetMatrix);

  py::class_<AffineOperator, ModPiece, WorkPiece, std::shared_ptr<AffineOperator>>(m, "AffineOperator")
    .def(py::init([](Eigen::MatrixXd const& Ain, Eigen::VectorXd const& bIn) {
        return std::make_shared<AffineOperator>(Ain, bIn);;
     }))
    .def(py::init<std::shared_ptr<LinearOperator>,Eigen::VectorXd>())
    .def("Linear", &AffineOperator::Linear)
    .def("Offset", &AffineOperator::Offset)
    .def("rows", &AffineOperator::rows)
    .def("cols", &AffineOperator::cols);

  py::class_<EigenLinearOperator<Eigen::MatrixXd>, LinearOperator, ModPiece, WorkPiece, std::shared_ptr<EigenLinearOperator<Eigen::MatrixXd>>> elo(m, "DenseLinearOperator");
  elo
    .def(py::init<Eigen::MatrixXd>())
    .def("Apply", &EigenLinearOperator<Eigen::MatrixXd>::Apply)
    .def("ApplyTranspose", &EigenLinearOperator<Eigen::MatrixXd>::ApplyTranspose);

  py::class_<IdentityOperator, LinearOperator, ModPiece, WorkPiece, std::shared_ptr<IdentityOperator>> io(m, "IdentityOperator");
  io
    .def(py::init<unsigned int>())
    .def("Apply", &IdentityOperator::Apply)
    .def("ApplyTranspose", &IdentityOperator::ApplyTranspose);

  /*
  py::class_<SparseLinearOperator, LinearOperator, std::shared_ptr<SparseLinearOperator>> slo(m, "SparseLinearOperator");
  slo
    .def(py::init<Eigen::SparseMatrix const&>())
    .def("Apply", (Eigen::MatrixXd (SparseLinearOperator::*)(Eigen::Ref<Eigen::MatrixXd> const&)) &SparseLinearOperator::Apply)
    .def("Apply", (void (SparseLinearOperator::*)(Eigen::Ref<Eigen::MatrixXd> const&, Eigen::Ref<Eigen::MatrixXd>)) &SparseLinearOperator::Apply)
    .def("ApplyTranspose", (Eigen::MatrixXd (SparseLinearOperator::*)(Eigen::Ref<Eigen::MatrixXd> const&)) &SparseLinearOperator::ApplyTranspose)
    .def("ApplyTranspose", (void (SparseLinearOperator::*)(Eigen::Ref<Eigen::MatrixXd> const&, Eigen::Ref<Eigen::MatrixXd>)) &SparseLinearOperator::ApplyTranspose);
  */

  py::class_<BlockDiagonalOperator, LinearOperator, std::shared_ptr<BlockDiagonalOperator>> bdOp(m, "BlockDiagonalOperator");
  bdOp
    .def(py::init<std::vector<std::shared_ptr<LinearOperator>> const&>())
    .def("Apply", &BlockDiagonalOperator::Apply)
    .def("ApplyTranspose", &BlockDiagonalOperator::ApplyTranspose)
    .def("GetMatrix", &BlockDiagonalOperator::GetMatrix)
    .def("GetBlock", &BlockDiagonalOperator::GetBlock)
    .def("GetBlocks", &BlockDiagonalOperator::GetBlocks);

  py::class_<BlockRowOperator, LinearOperator, std::shared_ptr<BlockRowOperator>> brOp(m, "BlockRowOperator");
  brOp
    .def(py::init<std::vector<std::shared_ptr<LinearOperator>> const&>())
    .def("Apply", &BlockRowOperator::Apply)
    .def("ApplyTranspose", &BlockRowOperator::ApplyTranspose)
    .def("GetMatrix", &BlockRowOperator::GetMatrix)
    .def("GetBlock", &BlockRowOperator::GetBlock)
    .def("GetBlocks", &BlockRowOperator::GetBlocks);

  py::class_<CompanionMatrix, LinearOperator, std::shared_ptr<CompanionMatrix>> compMat(m, "CompanionMatrix");
  compMat
    .def(py::init<Eigen::VectorXd const&>())
    .def("Apply", &CompanionMatrix::Apply)
    .def("ApplyTranspose", &CompanionMatrix::ApplyTranspose)
    .def("GetMatrix", &CompanionMatrix::GetMatrix);

  py::class_<ConcatenateOperator, LinearOperator, std::shared_ptr<ConcatenateOperator>> concatOp(m, "ConcatenateOperator");
  concatOp
    .def(py::init<std::vector<std::shared_ptr<LinearOperator>> const&, const int>())
    .def("Apply", &ConcatenateOperator::Apply)
    .def("ApplyTranspose", &ConcatenateOperator::ApplyTranspose)
    .def("GetMatrix", &ConcatenateOperator::GetMatrix)
    .def("VStack", &ConcatenateOperator::VStack)
    .def("HStack", &ConcatenateOperator::HStack);

  py::class_<DiagonalOperator, LinearOperator, std::shared_ptr<DiagonalOperator>> diagOp(m, "DiagonalOperator");
  diagOp
    .def(py::init<Eigen::VectorXd const&>())
    .def("Apply", &DiagonalOperator::Apply)
    .def("ApplyTranspose", &DiagonalOperator::ApplyTranspose)
    .def("GetMatrix", &DiagonalOperator::GetMatrix);

  py::class_<KroneckerProductOperator, LinearOperator, std::shared_ptr<KroneckerProductOperator>> kpOp(m, "KroneckerProductOperator");
  kpOp
    .def(py::init<std::shared_ptr<LinearOperator>, std::shared_ptr<LinearOperator>>())
    .def("Apply", &KroneckerProductOperator::Apply)
    .def("ApplyTranspose", &KroneckerProductOperator::ApplyTranspose);

  py::class_<ProductOperator, LinearOperator, std::shared_ptr<ProductOperator>> prOp(m, "ProductOperator");
  prOp
    .def(py::init<std::shared_ptr<LinearOperator>, std::shared_ptr<LinearOperator>>())
    .def("Apply", &ProductOperator::Apply)
    .def("ApplyTranspose", &ProductOperator::ApplyTranspose)
    .def("GetMatrix", &ProductOperator::GetMatrix);

  py::class_<SumOperator, LinearOperator, std::shared_ptr<SumOperator>> sumOp(m, "SumOperator");
  sumOp
    .def(py::init<std::shared_ptr<LinearOperator>, std::shared_ptr<LinearOperator>>())
    .def("Apply", &SumOperator::Apply)
    .def("ApplyTranspose", &SumOperator::ApplyTranspose)
    .def("GetMatrix", &SumOperator::GetMatrix);

  py::class_<ZeroOperator, LinearOperator, std::shared_ptr<ZeroOperator>> zeroOp(m, "ZeroOperator");
  zeroOp
    .def(py::init<int, int>())
    .def("Apply", &ZeroOperator::Apply)
    .def("ApplyTranspose", &ZeroOperator::ApplyTranspose);

};
