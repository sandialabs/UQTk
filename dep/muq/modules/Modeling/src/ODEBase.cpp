#include "MUQ/Modeling/ODEBase.h"

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts.
#include <cvodes/cvodes_spgmr.h> // prototypes & constants for CVSPGMR solver
#include <cvodes/cvodes_spbcgs.h> // prototypes & constants for CVSPBCG solver
#include <cvodes/cvodes_sptfqmr.h> // prototypes & constants for SPTFQMR solver
#include <cvodes/cvodes_dense.h> // prototype for CVDense

#include <sundials/sundials_types.h> // definition of type
#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP

#if MUQ_HAS_PARCER==1
#include <nvector/nvector_parallel.h>
#endif

namespace pt = boost::property_tree;
using namespace muq::Modeling;

#if MUQ_HAS_PARCER==1
ODEBase::ODEBase(std::shared_ptr<ModPiece> const& rhs, Eigen::VectorXi const& inputSizes, Eigen::VectorXi const& outputSizes, pt::ptree const& pt, std::shared_ptr<parcer::Communicator> const& comm) :
#else
ODEBase::ODEBase(std::shared_ptr<ModPiece> const& rhs, Eigen::VectorXi const& inputSizes, Eigen::VectorXi const& outputSizes, pt::ptree const& pt) :
#endif
  ModPiece(inputSizes, outputSizes),
  rhs(rhs),
  reltol(pt.get<double>("RelativeTolerance", 1.0e-8)),
  abstol(pt.get<double>("AbsoluteTolerance", 1.0e-8)), maxStepSize(pt.get<double>("MaxStepSize", 1.0)),
  maxNumSteps(pt.get<unsigned int>("MaxNumSteps", 500)),
  autonomous(pt.get<bool>("Autonomous", true)),
  checkPtGap(pt.get<unsigned int>("CheckPointGap", 50))
#if MUQ_HAS_PARCER==1
  , globalSize(pt.get<unsigned int>("GlobalSize", std::numeric_limits<unsigned int>::quiet_NaN())),
  comm(comm)
#endif
   {
    // we must know the number of inputs for the rhs and it must have at least one (the state)
    assert(rhs->numInputs>0);

  // determine the multistep method and the nonlinear solver
  const std::string& multiStepMethod = pt.get<std::string>("MultistepMethod", "BDF");
  assert(multiStepMethod.compare("Adams")==0 || multiStepMethod.compare("BDF")==0); // Adams or BDF
  multiStep = (multiStepMethod.compare("BDF")==0) ? CV_BDF : CV_ADAMS;
  const std::string& nonlinearSolver = pt.get<std::string>("NonlinearSolver", "Newton");
  assert(nonlinearSolver.compare("Iter")==0 || nonlinearSolver.compare("Newton")==0); // Iter or Newton
  solveMethod = (nonlinearSolver.compare("Newton")==0) ? CV_NEWTON : CV_FUNCTIONAL;

  // determine the method of thelinear solver
  const std::string& linearSolver = pt.get<std::string>("LinearSolver", "Dense");
  if( linearSolver.compare("Dense")==0 ) {
    slvr = LinearSolver::Dense;
  } else if( linearSolver.compare("SPGMR")==0 ) {
    slvr = LinearSolver::SPGMR;
  } else if( linearSolver.compare("SPBCG")==0 ) {
    slvr = LinearSolver::SPBCG;
  } else if( linearSolver.compare("SPTFQMR")==0 ) {
    slvr = LinearSolver::SPTFQMR;
  } else {
    std::cerr << "\nInvalid CVODES linear solver type.  Options are Dense, SPGMR, SPBCG, or SPTFQMR\n\n";
    assert(false);
  }
}

ODEBase::~ODEBase() {}

bool ODEBase::CheckFlag(void* flagvalue, std::string const& funcname, unsigned int const opt) const {
  // there are only two options
  assert(opt==0 || opt==1);

  // check if Sundials function returned nullptr pointer - no memory allocated
  if( opt==0 && flagvalue==nullptr ) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname.c_str());

    // return failure
    return false;
  }

  // check if flag<0
  if( opt==1 ) {
    // get int value
    int *errflag = (int *) flagvalue;

    if( *errflag<0 ) { // negative indicates Sundials error
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname.c_str(), *errflag);

      // return failure
      return false;
    }
  }

  // return success
  return true;
}

void ODEBase::ErrorHandler(int error_code, const char *module, const char *function, char *msg, void *user_data) {}

int ODEBase::CreateSolverMemoryB(void* cvode_mem, double const timeFinal, N_Vector const& lambda, N_Vector const& nvGrad, std::shared_ptr<ODEData> data) const {
  int indexB, flag;

  // initialize the adjoint solver
  flag = CVodeAdjInit(cvode_mem, checkPtGap, CV_HERMITE);
  assert(CheckFlag(&flag, "CVodeAdjInit", 1));

  // creat adjoint solver
  flag = CVodeCreateB(cvode_mem, multiStep, solveMethod, &indexB);
  assert(CheckFlag(&flag, "CVodeCreateB", 1));

  // set the pointer to user-defined data
  assert(data);
  flag = CVodeSetUserDataB(cvode_mem, indexB, data.get());
  assert(CheckFlag(&flag, "CVodeSetUserDataB", 1));

  // set the adjoint right hand side
  flag = CVodeInitB(cvode_mem, indexB, AdjointRHS, timeFinal, lambda);
  assert(CheckFlag(&flag, "CVodeInitB", 1));

  // specify the relative and absolute tolerances
  flag = CVodeSStolerancesB(cvode_mem, indexB, reltol, abstol);
  assert(CheckFlag(&flag, "CVodeSStolerancesB", 1));

  // set the maximum time step size
  flag = CVodeSetMaxStepB(cvode_mem, indexB, maxStepSize);
  assert(CheckFlag(&flag, "CVodeSetMaxStepB", 1));

  flag = CVodeSetMaxNumStepsB(cvode_mem, indexB, 10000000);
  assert(CheckFlag(&flag, "CVodesSetMaxNumStepsB", 1));

  // determine which linear solver to use
  if( slvr==LinearSolver::Dense ) { // dense linear solver
    // specify the dense linear solver
#if MUQ_HAS_PARCER==1
    flag = CVDenseB(cvode_mem, indexB, comm? NV_GLOBLENGTH_P(lambda) : NV_LENGTH_S(lambda));
#else
    flag = CVDenseB(cvode_mem, indexB, NV_LENGTH_S(lambda));
#endif
    assert(CheckFlag(&flag, "CVDenseB", 1));

    // set the Jacobian routine to Jac (user-supplied)
    flag = CVDlsSetDenseJacFnB(cvode_mem, indexB, AdjointJacobian);
    assert(CheckFlag(&flag, "CVDlsSetDenseJacFnB", 1));
  } else { // sparse linear solver
    if( slvr==LinearSolver::SPGMR ) {
      flag = CVSpgmrB(cvode_mem, indexB, 0, 0);
      assert(CheckFlag(&flag, "CVSpgmrB", 1));
    } else if( slvr==LinearSolver::SPBCG ) {
      flag = CVSpbcgB(cvode_mem, indexB, 0, 0);
      assert(CheckFlag(&flag, "CVSpbcgB", 1));
    } else if( slvr==LinearSolver::SPTFQMR ) {
      flag = CVSptfqmrB(cvode_mem, indexB, 0, 0);
      assert(CheckFlag(&flag, "CVSptfqmrB", 1));
    } else {
      std::cerr << "\nInvalid CVODES linear solver type.  Options are Dense, SPGMR, SPBCG, or SPTFQMR\n\n";
      assert(false);
    }

    // set the Jacobian-times-vector function
    flag = CVSpilsSetJacTimesVecFnB(cvode_mem, indexB, AdjointJacobianAction);
    assert(CheckFlag(&flag, "CVSpilsSetJacTimesVecFnB", 1));
  }

  flag = CVodeQuadInitB(cvode_mem, indexB, AdjointQuad, nvGrad);
  assert(CheckFlag(&flag, "CVodeQuadInitB", 1));

  flag = CVodeSetQuadErrCon(cvode_mem, true);
  assert(CheckFlag(&flag, "CVodeSetQuadErrCon", 1));

  flag = CVodeQuadSStolerancesB(cvode_mem, indexB, reltol, abstol);
  assert(CheckFlag(&flag, "CVodeQuadSStolerancesB", 1));

  return indexB;
}

void ODEBase::CreateSolverMemory(void* cvode_mem, N_Vector const& state, std::shared_ptr<ODEData> data) const {
  // a flag used for error checking
  int flag;

  // set the pointer to user-defined data
  assert(data);
  flag = CVodeSetUserData(cvode_mem, data.get());
  assert(CheckFlag(&flag, "CVodeSetUserData", 1));

  // tell the solver how to deal with errors
  flag = CVodeSetErrHandlerFn(cvode_mem, ErrorHandler, data.get());
  assert(CheckFlag(&flag, "CVodeSetErrHandlerFun", 1));

  // tell the solver how to evaluate the rhs
  flag = CVodeInit(cvode_mem, EvaluateRHS, 0.0, state);
  assert(CheckFlag(&flag, "CVodeInit", 1));

  // specify the relative and absolute tolerances
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  assert(CheckFlag(&flag, "CVodeSStolerances", 1));

  // set the maximum time step size
  flag = CVodeSetMaxStep(cvode_mem, maxStepSize);
  assert(CheckFlag(&flag, "CVodeSetMaxStep", 1));

  flag = CVodeSetMaxNumSteps(cvode_mem, maxNumSteps);
  assert(CheckFlag(&flag, "CVodesSetMaxNumSteps", 1));

  // determine which linear solver to use
  if( slvr==LinearSolver::Dense ) { // dense linear solver
    // specify the dense linear solver
#if MUQ_HAS_PARCER==1
    flag = CVDense(cvode_mem, comm? NV_GLOBLENGTH_P(state) : NV_LENGTH_S(state));
#else
    flag = CVDense(cvode_mem, NV_LENGTH_S(state));
#endif
    assert(CheckFlag(&flag, "CVDense", 1));

    // set the Jacobian routine to Jac (user-supplied)
    flag = CVDlsSetDenseJacFn(cvode_mem, RHSJacobian);
    assert(CheckFlag(&flag, "CVDlsSetDenseJacFn", 1));
  } else { // sparse linear solver
    if( slvr==LinearSolver::SPGMR ) {
      flag = CVSpgmr(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSpgmr", 1));
    } else if( slvr==LinearSolver::SPBCG ) {
      flag = CVSpbcg(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSpbcg", 1));
    } else if( slvr==LinearSolver::SPTFQMR ) {
      flag = CVSptfqmr(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSptfqmr", 1));
    } else {
      std::cerr << "\nInvalid CVODES linear solver type.  Options are Dense, SPGMR, SPBCG, or SPTFQMR\n\n";
      assert(false);
    }

    // set the Jacobian-times-vector function
    flag = CVSpilsSetJacTimesVecFn(cvode_mem, RHSJacobianAction);
    assert(CheckFlag(&flag, "CVSpilsSetJacTimesVecFn", 1));
  }
}

#define NVectorMap(eigenmap, comm, vec) Eigen::Map<Eigen::VectorXd> eigenmap( \
  comm? NV_DATA_P(vec) : NV_DATA_S(vec), \
  comm? NV_LOCLENGTH_P(vec) : NV_LENGTH_S(vec));

int ODEBase::EvaluateRHS(realtype time, N_Vector state, N_Vector statedot, void *user_data) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

  // evaluate the rhs
#if MUQ_HAS_PARCER==1
  NVectorMap(statedotref, data->comm, statedot);
#else
  Eigen::Map<Eigen::VectorXd> statedotref(NV_DATA_S(statedot), NV_LENGTH_S(statedot));
#endif
  statedotref = data->rhs->Evaluate(data->inputs) [0];

  return 0;
}

int ODEBase::AdjointRHS(realtype time, N_Vector state, N_Vector lambda, N_Vector deriv, void *user_data) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

#if MUQ_HAS_PARCER==1
  NVectorMap(lambdaMap, data->comm, lambda);
#else
  Eigen::Map<Eigen::VectorXd> lambdaMap(NV_DATA_S(lambda), NV_LENGTH_S(lambda));
#endif
  const Eigen::VectorXd& lam = lambdaMap;

#if MUQ_HAS_PARCER==1
  NVectorMap(derivMap, data->comm, deriv);
#else
  Eigen::Map<Eigen::VectorXd> derivMap(NV_DATA_S(deriv), NV_LENGTH_S(deriv));
#endif

#if MUQ_HAS_PARCER==1
  int globalind = 0;
  int localSize = 0;
  if( data->comm ) {
    localSize = NV_LOCLENGTH_P(state);
    for( unsigned int i=0; i<data->comm->GetSize(); ++i ) {
      int temp = localSize;
      data->comm->Bcast(temp, i);
      globalind += i<data->comm->GetRank()? temp : 0;
    }
  } else {
    localSize = NV_LENGTH_S(state);
  }
#else
  int globalind = 0;
  int localSize = NV_LENGTH_S(state);
#endif
  derivMap = -1.0*data->rhs->Gradient(0, data->autonomous? 0 : 1, data->inputs, lam).segment(globalind, localSize);

  return 0;
}

int ODEBase::RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

#if MUQ_HAS_PARCER==1
  NVectorMap(vref, data->comm, v);
#else
  Eigen::Map<Eigen::VectorXd> vref(NV_DATA_S(v), NV_LENGTH_S(v));
#endif
  const Eigen::VectorXd& veig = vref;

  // compute the jacobain wrt the state
#if MUQ_HAS_PARCER==1
  NVectorMap(Jvref, data->comm, Jv);
#else
  Eigen::Map<Eigen::VectorXd> Jvref(NV_DATA_S(Jv), NV_LENGTH_S(Jv));
#endif
  Jvref = data->rhs->ApplyJacobian(0, data->autonomous? 0 : 1, data->inputs, veig);

  return 0;
}

int ODEBase::AdjointJacobianAction(N_Vector target, N_Vector output, realtype time, N_Vector state, N_Vector lambda, N_Vector adjRhs, void *user_data, N_Vector tmp) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

#if MUQ_HAS_PARCER==1
  NVectorMap(targetMap, data->comm, target);
#else
  Eigen::Map<Eigen::VectorXd> targetMap(NV_DATA_S(target), NV_LENGTH_S(target));
#endif
  const Eigen::VectorXd& targ = targetMap;
#if MUQ_HAS_PARCER==1
  NVectorMap(outputMap, data->comm, output);
#else
  Eigen::Map<Eigen::VectorXd> outputMap(NV_DATA_S(output), NV_LENGTH_S(output));
#endif

  outputMap = -1.0*data->rhs->Gradient(0, data->autonomous? 0 : 1, data->inputs, targ);

  return 0;
}

int ODEBase::RHSJacobian(long int N, realtype time, N_Vector state, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

  // evaluate the jacobian
  Eigen::Map<Eigen::MatrixXd> jacref(jac->data, jac->M, jac->N);
  jacref = data->rhs->Jacobian(0, data->autonomous? 0 : 1, data->inputs);

  return 0;
}

int ODEBase::AdjointJacobian(long int N, realtype time, N_Vector state, N_Vector lambda, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

  // evaluate the jacobian
  Eigen::Map<Eigen::MatrixXd> jacref(jac->data, jac->M, jac->N);
  jacref = -1.0*data->rhs->Jacobian(0, data->autonomous? 0 : 1, data->inputs).transpose();

  return 0;
}

int ODEBase::AdjointQuad(realtype time, N_Vector state, N_Vector lambda, N_Vector quadRhs, void *user_data) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, state);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
#endif
  data->UpdateInputs(stateref, time);

#if MUQ_HAS_PARCER==1
  NVectorMap(lambdaMap, data->comm, lambda);
#else
  Eigen::Map<Eigen::VectorXd> lambdaMap(NV_DATA_S(lambda), NV_LENGTH_S(lambda));
#endif
  const Eigen::VectorXd& lam = lambdaMap;

#if MUQ_HAS_PARCER==1
  NVectorMap(nvQuadMap, data->comm, quadRhs);
#else
  Eigen::Map<Eigen::VectorXd> nvQuadMap(NV_DATA_S(quadRhs), NV_LENGTH_S(quadRhs));
#endif

  nvQuadMap = data->rhs->Gradient(0, data->autonomous? data->wrtIn : data->wrtIn+1, data->inputs, lam);

  return 0;
}

void ODEBase::SetUpSensitivity(void *cvode_mem, unsigned int const paramSize, N_Vector *sensState) const {
  // initialze the forward sensitivity solver
  int flag = CVodeSensInit(cvode_mem, paramSize, CV_SIMULTANEOUS, ForwardSensitivityRHS, sensState);
  assert(CheckFlag(&flag, "CVodeSensInit", 1));

  // set sensitivity tolerances
  Eigen::VectorXd absTolVec = Eigen::VectorXd::Constant(paramSize, abstol);
  flag = CVodeSensSStolerances(cvode_mem, reltol, absTolVec.data());
  assert(CheckFlag(&flag, "CVodeSensSStolerances", 1));

  // error control strategy should test the sensitivity variables
  flag = CVodeSetSensErrCon(cvode_mem, true);
  assert(CheckFlag(&flag, "CVodeSetSensErrCon", 1));
}

int ODEBase::ForwardSensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
#if MUQ_HAS_PARCER==1
  NVectorMap(stateref, data->comm, y);
#else
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(y), NV_LENGTH_S(y));
#endif
  data->UpdateInputs(stateref, time);

  for( unsigned int i=0; i<Ns; ++i ) {
#if MUQ_HAS_PARCER==1
    NVectorMap(rhsMap, data->comm, ySdot[i]);
    NVectorMap(sensMap, data->comm, ys[i]);
#else
    Eigen::Map<Eigen::VectorXd> rhsMap(NV_DATA_S(ySdot[i]), stateref.size());
    Eigen::Map<Eigen::VectorXd> sensMap(NV_DATA_S(ys[i]), stateref.size());
#endif

    rhsMap = data->rhs->ApplyJacobian(0, data->autonomous? 0 : 1, data->inputs, (Eigen::VectorXd)sensMap);
  }

  // the derivative of the rhs wrt the parameter
  const Eigen::MatrixXd& dfdp = data->rhs->Jacobian(0, data->autonomous? data->wrtIn : data->wrtIn+1, data->inputs);

  // now, loop through and fill in the rhs vectors stored in ySdot
  for( unsigned int i=0; i<Ns; ++i ) {
#if MUQ_HAS_PARCER==1
    Eigen::Map<Eigen::VectorXd> rhsMap(data->comm? NV_DATA_P(ySdot[i]) : NV_DATA_S(ySdot[i]), stateref.size());
    Eigen::Map<Eigen::VectorXd> sensMap(data->comm? NV_DATA_P(ys[i]) : NV_DATA_S(ys[i]), stateref.size());
#else
    Eigen::Map<Eigen::VectorXd> rhsMap(NV_DATA_S(ySdot[i]), stateref.size());
    Eigen::Map<Eigen::VectorXd> sensMap(NV_DATA_S(ys[i]), stateref.size());
#endif
    rhsMap += dfdp.col(i);
  }

  return 0;
}
