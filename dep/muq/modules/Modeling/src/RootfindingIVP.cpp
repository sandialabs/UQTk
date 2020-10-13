#include "MUQ/Modeling/RootfindingIVP.h"

// Sundials includes
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;

RootfindingIVP::RootfindingIVP(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : ODEBase(rhs, pt, algebra), root(root), maxSteps(pt.get<unsigned int>("Rootfinder.MaxSteps", (int)1e10)), maxTime(pt.get<double>("Rootfinder.MaxTime", 1.0e3)), maxErrorTests(pt.get<unsigned int>("Rootfinder.MaxErrorTests", 100)) {
  // we must know the number of inputs for both the rhs and the root and they must have at least one (the state)
  assert(rhs->numInputs>0);
  assert(root->numInputs>0);

  // we must know the number of outputs for the root
  assert(root->numOutputs>0);

  // set the input and output types
  UpdateInputOutputTypes();
}

RootfindingIVP::~RootfindingIVP() {}

Eigen::VectorXi RootfindingIVP::FindRoot(ref_vector<boost::any> const& inputs, int const wrtIn, int const wrtOut, DerivativeMode const& mode) {
  // the number of inputs must be at least the the number of inputs required by the rhs and the root
  assert(inputs.size()>=rhs->numInputs+root->numInputs-1-(autonomous? 0 : 1));

  // clear the results
  ClearResults();

  // create the state vector (have to do a hard copy --- N_Vector is a pointer to the data, the pointer has been declared const, not the data)
  N_Vector state;
  DeepCopy(state, boost::any_cast<const N_Vector&>(inputs[0]));

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, root, inputs, autonomous, wrtIn, wrtOut);
  if( !autonomous ) {
    const boost::any time = 0.0;
    data->inputs.insert(data->inputs.begin(), time);
  }

  // set solver to null
  void* cvode_mem = nullptr;

  // create the solver memory
  cvode_mem = CVodeCreate(multiStep, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));

  // initialize the solver
  CreateSolverMemory(cvode_mem, state, data);

  // sensState will hold several N_Vectors (because it's a Jacobian)
  N_Vector *sensState = nullptr;
  int paramSize = -1;

  if( wrtIn>=0 && wrtOut==0 && wrtIn<rhs->numInputs ) { // we are computing the derivative wrt one of the rhs parameters
    paramSize = algebra->Size(inputs[wrtIn]); // the dimension of the parameter

    // set up sensitivity vector
    sensState = N_VCloneVectorArray_Serial(paramSize, state);
    assert(CheckFlag((void *)sensState, "N_VCloneVectorArray_Serial", 0));
    
    // initialize the sensitivies to zero
    for( int is=0; is<paramSize; ++is ) {
      N_VConst(0.0, sensState[is]);
    }
    
    // set up solver for sensitivity
    SetUpSensitivity(cvode_mem, paramSize, sensState);

    // initalize the derivative information
    InitializeDerivative(1, NV_LENGTH_S(state), paramSize, mode);
  }

  // tell the solver how evaluate the root
  int flag = CVodeRootInit(cvode_mem, root->numOutputs, EvaluateRoot);
  assert(CheckFlag(&flag, "CVodeRootInit", 1));

  // set the maximum number of steps
  flag = CVodeSetMaxNumSteps(cvode_mem, maxSteps);
  assert(CheckFlag(&flag, "CVodeSetMaxNumSteps", 1));

  // set the maximum number of error test failures
  flag = CVodeSetMaxErrTestFails(cvode_mem, maxErrorTests);
  assert(CheckFlag(&flag, "CVodeSetMaxErrorTestFails", 1));
  
  // set the intial time
  double t = 0.0;

  // each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
  if( inputs.size()>rhs->numInputs+root->numInputs-1 ) {
    ref_vector<boost::any> outputTimes(inputs.begin()+rhs->numInputs+root->numInputs-1, inputs.end());
    std::vector<std::pair<unsigned int, unsigned int> > timeIndices = TimeIndices(outputTimes);    

    // first: the next time to integrate to, second: the output index
    std::pair<double, int> nextTime;

    while( NextTime(nextTime, timeIndices, outputTimes) ) {
      // if we have to move forward --- i.e., not at the initial time or another output does not need the current time
      if( std::fabs(nextTime.first-t)>1.0e-14 ) { // we have to move forward in time
	flag = CVode(cvode_mem, nextTime.first, state, &t, CV_NORMAL);
	assert(CheckFlag(&flag, "CVode", 1));
      }

      // if we are at the root, but the root time is not a time we asked for do not record!
      if( std::fabs(nextTime.first-t)<1.0e-14 ) { // we do not have to move forward in time
	// save the result at this timestep
	if( timeIndices[nextTime.second].second>1 ) { // the output has more than one compnent ...
	  // .. save the current state as an element in the vector
	  DeepCopy(boost::any_cast<std::vector<N_Vector>&>(outputs[nextTime.second]) [timeIndices[nextTime.second].first-1], state);
	} else { // the out has one component ...
	  // ... save the current state (not inside a vector)
	  DeepCopy(boost::any_cast<N_Vector&>(outputs[nextTime.second]), state);
	}
      }

      if( flag==CV_ROOT_RETURN ) { break; }
    }
  }

  if( flag!=CV_ROOT_RETURN ) {
    // integrate forward in time
    flag = CVode(cvode_mem, maxTime, state, &t, CV_NORMAL);
    assert(CheckFlag(&flag, "CVode", 1));
  }

  // make sure we found a root
  assert(flag==CV_ROOT_RETURN);
      
  // set the output
  outputs.insert(outputs.begin(), t);
  outputs.insert(outputs.begin(), state);

  if( sensState ) {
    // get the sensitivity
    flag = CVodeGetSens(cvode_mem, &t, sensState);
    assert(CheckFlag(&flag, "CVodeGetSens", 1));
    
    // create a new jacobian --- shallow copy
    DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian);
    for( unsigned int col=0; col<paramSize; ++col ) {
      DENSE_COL(jac, col) = NV_DATA_S(sensState[col]);
      NV_DATA_S(sensState[col]) = nullptr;
    }
    
    // delete sensState
    N_VDestroyVectorArray_Serial(sensState, paramSize);
  }
 
  // retrieve information about the root
  Eigen::VectorXi rootsfound(root->numOutputs);
  flag = CVodeGetRootInfo(cvode_mem, rootsfound.data());
  assert(CheckFlag(&flag, "CVodeGetRootInfo", 1));

  // free integrator memory (don't destory the state --- outputs is a pointer to it)
  CVodeFree(&cvode_mem);

  return rootsfound;
}

void RootfindingIVP::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // find the first root to one of the root outputs
  FindRoot(inputs);
}

void RootfindingIVP::JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) {
  if( wrtOut>0 ) {
    std::cerr << std::endl << "ERROR: Cannot compute derivative infomration of muq::Modeling::RootfindingIVP for outputs other than the state at the root (the first output)" << std::endl << std::endl;
    
    assert(wrtOut==0);
  }
  
  // find the first root to one of the root outputs
  const Eigen::VectorXi& rootsfound = FindRoot(inputs, wrtIn, wrtOut);
  unsigned int ind = 0;
  for( ; ind<rootsfound.size(); ++ind ) {
    if( rootsfound(ind)!=0 ) { break; } 
  }

  // the inputs to the rhs function
  ref_vector<boost::any> rhsIns;
  rhsIns.push_back(outputs[0]);
  rhsIns.insert(rhsIns.end(), inputs.begin()+1, inputs.begin()+rhs->numInputs);

  // the inputs to the root function
  ref_vector<boost::any> rootIns;
  rootIns.push_back(outputs[0]);
  rootIns.insert(rootIns.end(), inputs.begin()+rhs->numInputs, inputs.begin()+rhs->numInputs+root->numInputs-1);

  // the derivative of the root function at the final time --  derivative of the root wrt the state
  const boost::any& dgdy = root->Jacobian(0, ind, rootIns);
  const DlsMat& dgdyref = boost::any_cast<const DlsMat&>(dgdy); // 1 x stateSize matrix

  // evaluate the right hand side
  const std::vector<boost::any>& f = rhs->Evaluate(rhsIns); // stateSize x 1 vector
  const N_Vector& fref = boost::any_cast<N_Vector>(f[0]);
  assert(NV_LENGTH_S(fref)==dgdyref->N);

  double dum = 0.0;
  for( unsigned int i=0; i<dgdyref->N; ++i ) {
    dum += DENSE_ELEM(dgdyref, 0, i)*NV_Ith_S(fref, i);
  }

  if( wrtIn>=rhs->numInputs ) { // derivative wrt root parameter
    boost::any dgdpara = root->Jacobian(wrtIn-rhs->numInputs+1, ind, rootIns);
    const DlsMat& dgdpararef = boost::any_cast<const DlsMat&>(dgdpara); // 1 x paramSize matrix
    DenseScale(-1.0/dum, dgdpararef);

    jacobian = NewDenseMat(NV_LENGTH_S(fref), dgdpararef->N);
    DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian); // stateSize x paramSize matrix
    SetToZero(jac);

    // update the jacobian jac += f.transpose()*dgdpara (outer product)
    for( unsigned int i=0; i<jac->M; ++i ) {
      for( unsigned int j=0; j<jac->N; ++j ) {
	DENSE_ELEM(jac, i, j) = NV_Ith_S(fref, i)*DENSE_ELEM(dgdpararef, 0, j);
      }
    }

    return;
  }

  // get a reference to the jacobian
  DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian); // stateSize x paramSize matrix
  assert(jac->M==dgdyref->N);
  assert(NV_LENGTH_S(fref)==jac->M);
  
  if( wrtIn==0 ) { // derivative wrt the initial conditions
    // add the identity
    AddIdentity(jac);
  }

  DenseScale(-1.0/dum, dgdyref);

  // the derivative of the final time wrt to the parameter; dtfdpara = -jac_{state}(root)*jac/dot(jac_{state}(root), f)
  N_Vector dtfdic = N_VNew_Serial(jac->N);
  for( unsigned int i=0; i<jac->N; ++i ) {
    NV_Ith_S(dtfdic, i) = 0.0;
    for( unsigned int j=0; j<jac->M; ++j ) {
      NV_Ith_S(dtfdic, i) += DENSE_ELEM(dgdyref, 0, j)*DENSE_ELEM(jac, j, i);
    }
  }

  // update the jacobian jac += f.transpose()*dtfdpara (outer product)
  for( unsigned int i=0; i<jac->M; ++i ) {
    for( unsigned int j=0; j<jac->N; ++j ) {
      DENSE_ELEM(jac, i, j) += NV_Ith_S(fref, i)*NV_Ith_S(dtfdic, j);
    }
  }
  
  // destroy temp vectors/matrices
  N_VDestroy(dtfdic);
}

void RootfindingIVP::UpdateInputOutputTypes() {
  if( autonomous ) {
    // the type of the first input (the state) for the rhs and the root
    assert(rhs->InputType(0, false).compare(root->InputType(0, false))==0);
  } else { // non-autonomous
    // the input type of the second input of the rhs is the input of the first input of the root
    assert(rhs->InputType(1, false).compare(root->InputType(0, false))==0);
  }

  // the second set of input parameters are the parameters for the root
  for( auto intype : root->InputTypes() ) {
    if( intype.first==0 ) { // we've already set the state type
      continue;
    }

    inputTypes[rhs->numInputs+intype.first-1] = intype.second;
  }

  // the first output type is the state
  outputTypes[0] = root->InputType(0, false);

  // the second output is the time where the root is reached
  outputTypes[1] = typeid(double).name();
}

int RootfindingIVP::EvaluateRoot(realtype t, N_Vector state, realtype *root, void *user_data) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->root);

  // set the state input
  const boost::any& anyref = state;

  // the inputs the root function
  ref_vector<boost::any> rootins(data->inputs.begin()+data->rhs->numInputs, data->inputs.begin()+data->rhs->numInputs+data->root->numInputs-1);
  rootins.insert(rootins.begin(), anyref);

  // evaluate the root
  const std::vector<boost::any>& result = data->root->Evaluate(rootins);
  assert(result.size()==data->root->numOutputs);
  for( unsigned int i=0; i<result.size(); ++i ) {
    root[i] = boost::any_cast<double>(result[i]);
  }

  return 0;
}
