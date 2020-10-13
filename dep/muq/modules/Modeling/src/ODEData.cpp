#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

// construct basic ode data
#if MUQ_HAS_PARCER==1
ODEData::ODEData(std::shared_ptr<ModPiece> const& rhs, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn, std::shared_ptr<parcer::Communicator> const& comm) :
#else
ODEData::ODEData(std::shared_ptr<ModPiece> const& rhs, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn) :
#endif
rhs(rhs), autonomous(autonomous), wrtIn(wrtIn)
#if MUQ_HAS_PARCER==1
, comm(comm)
#endif
{
  inputs.resize(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs[i] = refinputs[i].get(); }
}

// construct with root function
ODEData::ODEData(std::shared_ptr<ModPiece> const& rhs, std::shared_ptr<ModPiece> const& root, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn) : rhs(rhs), root(root), autonomous(autonomous), wrtIn(wrtIn) {
  inputs.resize(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs[i] = refinputs[i].get(); }
}

void ODEData::UpdateInputs(Eigen::VectorXd const& state, double const time) {
  if( autonomous ) {
    inputs[0] = state;
    return;
  }

  inputs[0] = Eigen::VectorXd::Constant(1, time);
  inputs[1] = state;
}
