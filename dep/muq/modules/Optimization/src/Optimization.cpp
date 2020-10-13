#include "MUQ/Optimization/Optimization.h"

#include "MUQ/Utilities/Exceptions.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;

Optimization::Optimization(std::shared_ptr<CostFunction> cost, pt::ptree pt) :
  WorkPiece(cost->InputTypes(), cost->numInputs, std::vector<std::string>({typeid(Eigen::VectorXd).name(), typeid(double).name()})),
  algorithm(NLOptAlgorithm(pt.get<std::string>("Algorithm"))),
  opt(cost, 1),
  ftol_rel(pt.get<double>("Ftol.RelativeTolerance", 1.0e-8)),
  ftol_abs(pt.get<double>("Ftol.AbsoluteTolerance", 1.0e-8)),
  xtol_rel(pt.get<double>("Xtol.RelativeTolerance", 0.0)),
  xtol_abs(pt.get<double>("Xtol.AbsoluteTolerance", 0.0)),
  constraint_tol(pt.get<double>("ConstraintTolerance", 1.0e-8)),
  maxEvals(pt.get<unsigned int>("MaxEvaluations", 100)),
  minimize(pt.get<bool>("Minimize", true))
{}

Optimization::~Optimization() {}

double Optimization::Cost(unsigned int n, const double* x, double* grad, void* f_data) {
  CostHelper* opt = (CostHelper*)f_data;

  Eigen::Map<const Eigen::VectorXd> xmap(x, n);
  const Eigen::VectorXd& xeig = xmap;
  opt->inputs.at(0) = std::cref(xeig);

  if( grad ) {
    Eigen::Map<Eigen::VectorXd> gradmap(grad, n);
    const Eigen::VectorXd& gradeig = opt->cost->Gradient(0, opt->inputs, (Eigen::VectorXd)Eigen::VectorXd::Ones(1));
    gradmap = gradeig;
  }

  return opt->cost->Cost(opt->inputs);
}

nlopt_algorithm Optimization::NLOptAlgorithm(std::string const& alg) const {
  if( alg.compare("DIRECT")==0 ) { return NLOPT_GN_DIRECT; }
  if( alg.compare("DIRECTL")==0 ) { return NLOPT_GN_DIRECT_L; }
  if( alg.compare("CRS")==0 ) { return NLOPT_GN_CRS2_LM; }
  if( alg.compare("MLSL")==0 ) { return NLOPT_G_MLSL_LDS; }
  if( alg.compare("ISRES")==0 ) { return NLOPT_GN_ISRES; }
  if( alg.compare("COBYLA")==0 ) { return NLOPT_LN_COBYLA; }
  if( alg.compare("BOBYQA")==0 ) { return NLOPT_LN_BOBYQA; }
  if( alg.compare("NEWUOA")==0 ) { return NLOPT_LN_NEWUOA_BOUND; }
  if( alg.compare("PRAXIS")==0 ) { return NLOPT_LN_PRAXIS; }
  if( alg.compare("NM")==0 ) { return NLOPT_LN_NELDERMEAD; }
  if( alg.compare("SBPLX")==0 ) { return NLOPT_LN_SBPLX; }
  if( alg.compare("MMA")==0 ) { return NLOPT_LD_MMA; }
  if( alg.compare("SLSQP")==0 ) { return NLOPT_LD_SLSQP; }
  if( alg.compare("LBFGS")==0 ) { return NLOPT_LD_LBFGS; }
  if( alg.compare("PreTN")==0 ) { return NLOPT_LD_TNEWTON_PRECOND_RESTART; }
  if( alg.compare("LMVM")==0 ) { return NLOPT_LD_VAR2; }

  return NLOPT_LN_COBYLA;
}

void Optimization::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  opt.SetInputs(inputs);

  for( auto it=ineqConstraints.begin(); it!=ineqConstraints.end(); ++it ) { it->SetInputs(inputs); }
  for( auto it=eqConstraints.begin(); it!=eqConstraints.end(); ++it ) { it->SetInputs(inputs); }

  // create the optimizer
  auto solver = nlopt_create(algorithm, opt.cost->inputSizes(0));
  if( minimize )
    nlopt_set_min_objective(solver, Optimization::Cost, &opt);
  else
    nlopt_set_max_objective(solver, Optimization::Cost, &opt);

  for( std::vector<CostFunction>::size_type i=0; i<ineqConstraints.size(); ++i ) { nlopt_add_inequality_constraint(solver, Optimization::Cost, &ineqConstraints[i], constraint_tol); }
  for( std::vector<CostFunction>::size_type i=0; i<eqConstraints.size(); ++i ) { nlopt_add_equality_constraint(solver, Optimization::Cost, &eqConstraints[i], constraint_tol); }

  // set the tolerances
  nlopt_set_ftol_rel(solver, ftol_rel);
  nlopt_set_ftol_abs(solver, ftol_abs);
  nlopt_set_xtol_rel(solver, xtol_rel);
  nlopt_set_xtol_abs1(solver, xtol_abs);
  nlopt_set_maxeval(solver, maxEvals);

  const Eigen::VectorXd& xinit = boost::any_cast<Eigen::VectorXd const&>(inputs.at(0));

  outputs.resize(2);
  outputs[0] = xinit;
  Eigen::VectorXd& xopt = boost::any_cast<Eigen::VectorXd&>(outputs.at(0));

  double minf;
  const nlopt_result check = nlopt_optimize(solver, xopt.data(), &minf);
  if( check<0 ) { throw muq::ExternalLibraryError("NLOPT has failed with flag "+std::to_string(check)); }

  outputs[1] = minf;
}

std::pair<Eigen::VectorXd, double> Optimization::Solve(std::vector<Eigen::VectorXd> const& inputs) {
  std::vector<boost::any> anyins;
  anyins.reserve(inputs.size());

  for( auto in : inputs ) {
    anyins.push_back(in);
  }

  return Solve(ToRefVector(anyins));
}

std::pair<Eigen::VectorXd, double> Optimization::Solve(ref_vector<boost::any> const& inputs) {
  Evaluate(inputs);

  return std::pair<Eigen::VectorXd, double>(boost::any_cast<Eigen::VectorXd const&>(outputs[0]), boost::any_cast<double const>(outputs[1]));
}

void Optimization::AddInequalityConstraint(std::shared_ptr<CostFunction> ineq) {
  ineqConstraints.push_back(CostHelper(ineq, numInputs));
  UpdateInputs(ineq->numInputs-1);
}

void Optimization::AddEqualityConstraint(std::shared_ptr<CostFunction> eq) {
  eqConstraints.push_back(CostHelper(eq, numInputs));
  UpdateInputs(eq->numInputs-1);
}

void Optimization::UpdateInputs(unsigned int const numNewIns) {
  for( unsigned int i=numInputs; i<numInputs+numNewIns; ++i ) { inputTypes[i] = typeid(Eigen::VectorXd).name(); }
  numInputs += numNewIns;
}

Optimization::CostHelper::CostHelper(std::shared_ptr<CostFunction> cost, unsigned int const firstin) : cost(cost), firstin(firstin) {}

Optimization::CostHelper::~CostHelper() {}

void Optimization::CostHelper::SetInputs(ref_vector<boost::any> const& ins) {
  inputs.clear();

  // set initial condition
  inputs.push_back(std::cref(boost::any_cast<Eigen::VectorXd const&>(ins[0])));

  // set other inputs
  for( unsigned int i=firstin; i<firstin+cost->numInputs-1; ++i ) {
    std::cout << i << std::endl;
    inputs.push_back(std::cref(boost::any_cast<Eigen::VectorXd const&>(ins[i])));
  }
}
