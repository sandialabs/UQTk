#include "AllClassWrappers.h"

#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"

#include <pybind11/pybind11.h>


using namespace muq::Modeling;
namespace py = pybind11;

void muq::Modeling::PythonBindings::CwiseUnaryOperatorsWrapper(py::module &m) {
  py::class_<ExpOperator, ModPiece, std::shared_ptr<ExpOperator>>(m, "ExpOperator").def(py::init<unsigned int>());
  py::class_<CosOperator, ModPiece, std::shared_ptr<CosOperator>>(m, "CosOperator").def(py::init<unsigned int>());
  py::class_<SinOperator, ModPiece, std::shared_ptr<SinOperator>>(m, "SinOperator").def(py::init<unsigned int>());
  py::class_<AbsOperator, ModPiece, std::shared_ptr<AbsOperator>>(m, "AbsOperator").def(py::init<unsigned int>());
  py::class_<AcosOperator, ModPiece, std::shared_ptr<AcosOperator>>(m, "AcosOperator").def(py::init<unsigned int>());
  py::class_<AsinOperator, ModPiece, std::shared_ptr<AsinOperator>>(m, "AsinOperator").def(py::init<unsigned int>());
  py::class_<AtanOperator, ModPiece, std::shared_ptr<AtanOperator>>(m, "AtanOperator").def(py::init<unsigned int>());
  py::class_<AtanhOperator, ModPiece, std::shared_ptr<AtanhOperator>>(m, "AtanhOperator").def(py::init<unsigned int>());
  py::class_<CbrtOperator, ModPiece, std::shared_ptr<CbrtOperator>>(m, "CbrtOperator").def(py::init<unsigned int>());
  py::class_<CeilOperator, ModPiece, std::shared_ptr<CeilOperator>>(m, "CeilOperator").def(py::init<unsigned int>());
  py::class_<CoshOperator, ModPiece, std::shared_ptr<CoshOperator>>(m, "CoshOperator").def(py::init<unsigned int>());
  py::class_<DigammaOperator, ModPiece, std::shared_ptr<DigammaOperator>>(m, "DigammaOperator").def(py::init<unsigned int>());
  py::class_<ErfOperator, ModPiece, std::shared_ptr<ErfOperator>>(m, "ErfOperator").def(py::init<unsigned int>());
  py::class_<ErfcOperator, ModPiece, std::shared_ptr<ErfcOperator>>(m, "ErfcOperator").def(py::init<unsigned int>());
  py::class_<FloorOperator, ModPiece, std::shared_ptr<FloorOperator>>(m, "FloorOperator").def(py::init<unsigned int>());
  py::class_<InvLogitOperator, ModPiece, std::shared_ptr<InvLogitOperator>>(m, "InvLogitOperator").def(py::init<unsigned int>());
  py::class_<InvPhiOperator, ModPiece, std::shared_ptr<InvPhiOperator>>(m, "InvPhiOperator").def(py::init<unsigned int>());
  py::class_<InvSqrtOperator, ModPiece, std::shared_ptr<InvSqrtOperator>>(m, "InvSqrtOperator").def(py::init<unsigned int>());
  py::class_<InvSquareOperator, ModPiece, std::shared_ptr<InvSquareOperator>>(m, "InvSquareOperator").def(py::init<unsigned int>());
  py::class_<InvOperator, ModPiece, std::shared_ptr<InvOperator>>(m, "InvOperator").def(py::init<unsigned int>());
  py::class_<LogGammaOperator, ModPiece, std::shared_ptr<LogGammaOperator>>(m, "LogGammaOperator").def(py::init<unsigned int>());
  py::class_<LogInvLogitOperator, ModPiece, std::shared_ptr<LogInvLogitOperator>>(m, "LogInvLogitOperator").def(py::init<unsigned int>());
  py::class_<LogOperator, ModPiece, std::shared_ptr<LogOperator>>(m, "LogOperator").def(py::init<unsigned int>());
  py::class_<Log2Operator, ModPiece, std::shared_ptr<Log2Operator>>(m, "Log2Operator").def(py::init<unsigned int>());
  py::class_<Log10Operator, ModPiece, std::shared_ptr<Log10Operator>>(m, "Log10Operator").def(py::init<unsigned int>());
  py::class_<LogitOperator, ModPiece, std::shared_ptr<LogitOperator>>(m, "LogitOperator").def(py::init<unsigned int>());
  py::class_<PhiOperator, ModPiece, std::shared_ptr<PhiOperator>>(m, "PhiOperator").def(py::init<unsigned int>());
  py::class_<RoundOperator, ModPiece, std::shared_ptr<RoundOperator>>(m, "RoundOperator").def(py::init<unsigned int>());
  py::class_<SinhOperator, ModPiece, std::shared_ptr<SinhOperator>>(m, "SinhOperator").def(py::init<unsigned int>());
  py::class_<SqrtOperator, ModPiece, std::shared_ptr<SqrtOperator>>(m, "SqrtOperator").def(py::init<unsigned int>());
  py::class_<SquareOperator, ModPiece, std::shared_ptr<SquareOperator>>(m, "SquareOperator").def(py::init<unsigned int>());
  py::class_<TanOperator, ModPiece, std::shared_ptr<TanOperator>>(m, "TanOperator").def(py::init<unsigned int>());
  py::class_<TanhOperator, ModPiece, std::shared_ptr<TanhOperator>>(m, "TanhOperator").def(py::init<unsigned int>());
  py::class_<TgammaOperator, ModPiece, std::shared_ptr<TgammaOperator>>(m, "TgammaOperator").def(py::init<unsigned int>());
  py::class_<TrigammaOperator, ModPiece, std::shared_ptr<TrigammaOperator>>(m, "TrigammaOperator").def(py::init<unsigned int>());
}
