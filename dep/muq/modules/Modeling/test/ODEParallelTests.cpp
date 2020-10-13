#include "gtest/gtest.h"

#include <Eigen/Core>

#include <parcer/Communicator.h>
#include <parcer/Eigen.h>

#include "MUQ/Modeling/ODE.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;

class RHS : public ModPiece {
public:

  inline RHS(unsigned int const sizeLocal, std::shared_ptr<parcer::Communicator> comm) : ModPiece(Eigen::VectorXi::Constant(1, sizeLocal), Eigen::VectorXi::Constant(1, sizeLocal)), comm(comm) {
    A = Eigen::MatrixXd::Ones(comm->GetSize()*sizeLocal, comm->GetSize()*sizeLocal);
  }

  inline virtual ~RHS() = default;

private:

  inline Eigen::VectorXd GetGlobal(Eigen::VectorXd const& local) const {
    Eigen::VectorXd globalIn = Eigen::VectorXd::Zero(comm->GetSize()*inputSizes(0));
    for( unsigned int rank=0; rank<comm->GetSize(); ++rank ) {
      Eigen::VectorXd temp;
      if( rank==comm->GetRank() ) { temp = local; }
      comm->Bcast(temp, rank);

      globalIn.segment(rank*inputSizes(0), inputSizes(0)) = temp;
    }

    return globalIn;
  }

  inline virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    Eigen::VectorXd globalIn = GetGlobal(inputs[0]);
    globalIn = (A*globalIn).segment(comm->GetRank()*inputSizes(0), inputSizes(0));

    // set the output
    outputs.resize(1);
    outputs[0] = globalIn;
  }

  inline virtual void ApplyJacobianImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) override {
    // there is only one input and output
    assert(wrtOut==0);
    assert(wrtIn==0);

    Eigen::VectorXd globalVec = GetGlobal(vec);
    globalVec = (A*globalVec).segment(comm->GetRank()*inputSizes(0), inputSizes(0));

    jacobianAction = globalVec;
  }

  inline virtual void JacobianImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs) override {
    // there is only one input and output
    assert(wrtOut==0);
    assert(wrtIn==0);

    jacobian = A.block(comm->GetRank()*inputSizes(0), 0, inputSizes(0), comm->GetSize()*inputSizes(0));
  }

  inline virtual void GradientImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) override {
    // there is only one input and output
    assert(wrtOut==0);
    assert(wrtIn==0);

    Eigen::VectorXd globalVec = GetGlobal(vec);

    gradient = (Eigen::VectorXd)(A.transpose()*globalVec);
  }

  Eigen::MatrixXd A;

  std::shared_ptr<parcer::Communicator> comm;
};

class ParallelODETests : public::testing::Test {
  public:
  /// Default constructor
  ParallelODETests() {
    comm = std::make_shared<parcer::Communicator>();

    // create the right hand side
    rhs = std::make_shared<RHS>(sizeLocal, comm);

    // solver options
    pt.put<double>("ODE.RelativeTolerance", 1.0e-10);
    pt.put<double>("ODE.AbsoluteTolerance", 1.0e-10);
    pt.put<double>("ODE.MaxStepSize", 1.0);
    pt.put<unsigned int>("ODE.NumObservations", outTimes.size());
    pt.put<unsigned int>("ODE.MaxNumSteps", 1000);
    pt.put<unsigned int>("ODE.GlobalSize", comm->GetSize()*sizeLocal);
  }

  /// Default destructor
  virtual ~ParallelODETests() = default;

  virtual void TearDown() override {
    // create the ODE integrator
    auto ode = std::make_shared<ODE>(rhs, pt.get_child("ODE"), comm);

    // integrate the ODE
    const std::vector<Eigen::VectorXd>& result = ode->Evaluate(ic, outTimes);

    // each processor should compute the same solution
    for( unsigned int rank=0; rank<comm->GetSize(); ++rank ) {
      Eigen::VectorXd temp;
      if( rank==comm->GetRank() ) { temp = result[0]; }
      comm->Bcast(temp, rank);

      EXPECT_NEAR((temp-result[0]).norm(), 0.0, 1.0e-10);
    }

    // compute the Jacobian
    const Eigen::MatrixXd& jac = ode->Jacobian(0, 0, ic, outTimes);

    // compute the gradient
    const Eigen::VectorXd vec = Eigen::VectorXd::Ones(sizeLocal*outTimes.size());
    const Eigen::VectorXd& grad = ode->Gradient(0, 0, ic, outTimes, vec);

    // assemble the full jaocbian
    Eigen::MatrixXd fulljac(comm->GetSize()*sizeLocal*outTimes.size(), comm->GetSize()*sizeLocal);
    for( unsigned int i=0; i<comm->GetSize(); ++i ) {
      Eigen::MatrixXd temp = jac;
      comm->Bcast(temp, i);
      for( unsigned int t=0; t<outTimes.size(); ++t ) {
        fulljac.block(t*comm->GetSize()*sizeLocal+i*sizeLocal, 0, sizeLocal, comm->GetSize()*sizeLocal) = temp.block(t*sizeLocal, 0, sizeLocal, comm->GetSize()*sizeLocal);
      }
    }

    // the expected gradient is the action of the full Jacobian transpose
    const Eigen::VectorXd expectedGrad = fulljac.transpose()*Eigen::VectorXd::Ones(comm->GetSize()*sizeLocal*outTimes.size());

    // check the gradients
    EXPECT_NEAR((grad-expectedGrad).norm(), 0.0, 1.0e-4);
  }

  /// The communicator
  std::shared_ptr<parcer::Communicator> comm;

  /// The right hand side
  std::shared_ptr<RHS> rhs;

  /// Options for the ODE integrator
  pt::ptree pt;

  /// State dimension on one processor
  const unsigned int sizeLocal = 2;

  /// The initial condition
  Eigen::VectorXd ic = Eigen::VectorXd::Ones(sizeLocal);

  /// The output times
  const Eigen::VectorXd outTimes = Eigen::VectorXd::LinSpaced(10, 0.0, 1.0);
};

// note: not all parameter combinations are compatible with parallel runs.  This is because some of SUNDIAL's solvers are not compatible with the nvector_parallel library.

TEST_F(ParallelODETests, BDFIterMethod) {
  pt.put<std::string>("ODE.MultistepMethod", "BDF");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "Dense");
}

TEST_F(ParallelODETests, AdamsIterMethod) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "Dense");
}

TEST_F(ParallelODETests, SPGMR) {
  pt.put<std::string>("ODE.MultistepMethod", "BDF");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "SPGMR");
}

TEST_F(ParallelODETests, SPBCG) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "SPBCG");
}

TEST_F(ParallelODETests, SPTFQMR) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "SPTFQMR");
}
