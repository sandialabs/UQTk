#ifndef MODPIECE_H
#define MODPIECE_H

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Utilities/VariadicMacros.h"

#include <Eigen/Core>
#include <vector>

namespace muq{
  namespace Modeling{

  /**
    @ingroup Modeling
    @class ModPiece
    @brief Provides an abstract interface for defining vector-valued model components.
    @details Let \f$x_i \in \mathbb{R}^{N_{xi}}\f$ be a collection of \f$M_x\f$ vectors index by \f$i\f$.
             Each vector could have a different size.  Now let \f$ f_i(x_1,\ldots, x_{M_x}) \in \mathbb{R}^{N_{fi}}\f$
             denote a collection of \f$M_f\f$ model outputs, each of which can depend
             on all of the input vectors \f$x_i\f$.   The ModPiece class defines an
             interface for defining the mappings of the form \f$[x_1, \ldots, x_{M_x}] \rightarrow [f_1(x_1, \ldots, x_N), \ldots, f_{M_f}(x_1, \ldots, x_{M_{x}})] \f$.
             While this form is quite general, functions with a single input and a single output
             can easily be defined as a ModPiece with \f$M_x = M_f = 1\f$.

             The ModPiece class (and therefore all its children) contains two
             vectors, `inputSizes` and `outputSizes`, that are useful for investigating
             the number and size of the inputs and outputs.
             @code{.cpp}
             std::shared_ptr<ModPiece> m;

             // ... create a child of ModPiece and store in m

             int m_x = m->inputSizes.size();  // number of inputs
             int m_f = m->outputSizes.size(); // number of outputs
             int n_x1 = m->intputSizes(0); // The number of components in the first input vector
             int n_f1 = m->outputSizes(0); // The number of components in the first output vector
             @endcode

             <B>Defining a new model component</B>
             The `ModPiece` class defines a general modeling interface that allows algorithms (such as MCMC)
             to work with arbitrary models, without needing to know the model details.  Users can add
             create new models in two ways:
             <ul>
             <li> Combining existing ModPiece classes in a WorkGraph to create more complicated relationships. </li>
             <li> Creating a new ModPiece by defining a class that inherits from ModPiece. </li>
             </ul>

             As an example of creating a new ModPiece, consider two vectors \f$x_1\in\mathbb{R}^2\f$ and
             \f$x_2\in\mathbb{R}^2\f$ as well as three functions \f$f_1(x_1,x_2) = sin(x_1*x_2)\f$, \f$f_2(x_1,x_2) = cos(x_1*x_2)\f$,
             and \f$f_3(x_1,x_2) = tan(x_1*x_2)\f$.  The following code implements
             this relationship with a ModPiece.  Notice that the input sizes are \f$[2,2]\f$ and
             the output sizes (assuming componentwise products) are \f$[2,2,2]\f$.
             @code{.cpp}
            class MyTrigonometricPiece : public muq::Modeling::ModPiece
            {
            public:
              MyTrigonometricPiece() : muq::Modeling::ModPiece(2*Eigen::VectorXi::Ones(2),    // inputSizes = [2,2]
                                                               2*Eigen::VectorXi::Ones(3)){}; // outputSizes = [2,2,2]

            protected:
              virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override
              {
                Eigen::VectorXd const& x1 = inputs.at(0);
                Eigen::VectorXd const& x2 = inputs.at(1);

                Eigen::VectorXd xprod = x1.array()*x2.array();

                Eigen::VectorXd sinOut, cosOut, tanOut;
                sinOut = xprod.array().sin();
                cosOut = xprod.array().cos();
                tanOut = xprod.array().tan();


                // Resize the outputs vector (which lives in the ModPiece base class)
                outputs.resize(3);
                outputs.at(0) = sinOut;
                outputs.at(1) = cosOut;
                outputs.at(2) = tanOut;
              }
            };
             @endcode

             If MUQ was compiled with Python bindings, it is also possible to create
            new children of ModPiece in Python.  In Python, we use numpy arrays
            instead of Eigen::VectorXd.  With this, the above class in Python
            might take the form
            @code{.py}
            import pymuqModeling as mm
            import numpy as np

            class MyTrigonometricPiece(mm.PyModPiece):

                def __init__(self):
                    mm.PyModPiece.__init__(self, [2,2], [2,2,2])

                def EvaluateImpl(self, inputs):
                    x1 = inputs[0]
                    x2 = inputs[1]

                    xprod = x1*x2

                    sinOut = np.sin(xprod)
                    cosOut = np.cos(xprod)
                    tanOut = np.tan(xprod)

                    self.outputs = [sinOut, cosOut, tanOut]
            @endcode

            Once created, we can evaluate a ModPiece using the `Evaluate` function.
            Notice that the user-level `Evaluate` function and the `EvaluateImpl`
            function discussed above are different.   Direct calls to `EvaluateImpl`
            should be avoided (hence the reason for making it protected) because
            the `Evaluate` function performs additional checks to make sure
            the inputs and outputs are the correct size.  The `Evaluate` function
            also instruments the runtime and number of `EvaluateImpl` calls.

            For example, to evaluate the model defined in the MyTrigonometricPiece class,
            we might use something like this:
            @code{.cpp}
            int main()
            {

              auto piece = std::make_shared<MyTrigonometricPiece>();

              std::vector<Eigen::VectorXd> inputs(2);
              inputs.at(0) = Eigen::VectorXd::Random(2);
              inputs.at(1) = Eigen::VectorXd::Random(2);

              std::vector<Eigen::VectorXd> outputs = piece->Evaluate(inputs);
              // analogously, we can pass in individual Eigen::VectorXd's instead of the std::vector
              // outputs = piece->Evaluate(inputs.at(0), inputs.at(1));

              return 0;
            };
            @endcode
            Similarly, in python
            @code{.py}
            piece = MyTrigonometricPiece()

            inputs = [ np.random.rand(2), np.random.rand(2) ]
            outputs = piece.Evaluate(inputs)
            @endcode
            Note: in Python, a list of input vectors must be used.  We cannot pass
            each vector individually like we can in c++.

  */
  class ModPiece : public WorkPiece
  {

  public:

    ModPiece(Eigen::VectorXi const& inputSizes,
             Eigen::VectorXi const& outputSizes);

    virtual ~ModPiece() = default;


    /** @brief Get the average run time for one of the implemented methods.
     *   @details This function returns the average wall clock time (in milliseconds) for the EvaluateImpl, GradientImpl,
     * JacobianImpl, JacobianActionImpl, or HessianImp functions.
     *            If the function was never called, -1 is returned.
     *   @param[in] method The implemented function of interest.  Possible options are "Evaluate", "Gradient", "Jacobian",
     * "JacobianAction", or "Hessian"
     *   @return The average wall clock time in milli-seconds.
     */
    virtual double GetRunTime(const std::string& method="Evaluate") const override;
    virtual void ResetCallTime() override;

    /** @brief get the number of times one of the implemented methods has been called.
     *   @details This function returns the number of times the EvaluateImpl, GradientImpl, JacobianImpl,
     * JacobianActionImpl, or HessianImp functions have been called.
     *   @param[in] method The implemented function of interest.  Possible options are "Evaluate", "Gradient", "Jacobian",
     * "JacobianAction", or "Hessian"
     *   @return An integer with the number of calls.
     */
    virtual unsigned long int GetNumCalls(const std::string& method = "Evaluate") const override;

    /** @brief Evaluate the ModPiece
      @details This function takes in a vector of inputs and returns a vector of outputs.
      It calls the `EvaluateImpl` function after checking the size of the inputs.  It also
      measures the run time of `EvaluateImpl` and keeps track of the total number of `Evaluate` calls.
      @param[in] input A vector of inputs corresponding to \f$x_1,\ldots, x_{M_x}\f$
      @return A vector of outputs corresponding to \f$f_1(x_1, \ldots, x_N), \ldots, f_{M_f}(x_1, \ldots, x_{M_{x}})\f$
    */
    virtual std::vector<Eigen::VectorXd> const& Evaluate(std::vector<Eigen::VectorXd> const& input);
    virtual std::vector<Eigen::VectorXd> const& Evaluate(ref_vector<Eigen::VectorXd> const& input);
    VARIADIC_TO_REFVECTOR(Evaluate, Eigen::VectorXd, std::vector<Eigen::VectorXd> const&);

    /** @brief Compute the Gradient \f$J^Tv\f$.
      @details Consider some scalar-valued function \f$h(f_i) : \mathbb{R}^{N_{fi}} \rightarrow \mathbb{R}\f$ and let
      \f$v = \nabla_{f_i} h\f$ denote the gradient of \f$h\f$ with respect to \f$f_i\f$.
      Let \f$[x_1, \ldots, x_{M_x}]\f$ be the inputs to \f$f_i\f$, i.e.,
      \f[
      h(f_i) = h( f_i(x_1, \ldots, x_{M_x}) )
      \f]
      Given \f$v = \nabla_{f_i} h\f$, this function computes \f$ \nabla_{x_j} h\f$ for
      some user-specified index \f$j\f$.  Thus, this function computes one step
      of the chain rule.  If \f$J_{ij}\f$ is the Jacobian matrix of \f$f_i\f$ with respect
      to \f$x_j\f$.  Then, \f[ \nabla_{x_j} h = J_{ij}^T \left[ \nabla_{f_i} h \right]. \f]

      Like the `Evaluate` function, this function performs some checks on the inputs
      and then calls the `GradientImpl` function.  It also keeps track of how many
      times `Gradient` is called and the runtime of `GradientImpl`.

      If the GradientImpl function is not overriden by a child class, it will
      default to using finite differences by calling the GradientByFD function.

      @param[in] outputDimWrt The index \f$i\f$ describing the output of interest. (Think of the derivative \f$\partial f_i / \partial x_j\f$)
      @param[in] inputDimWrt The index \f$j\f$ describing which input we want to take derivatives with respect to.
      @param[in] input A vector of inputs corresponding to \f$x_1,\ldots, x_{M_x}\f$
      @param[in] sensitivity The vector \f$\nabla_{f_i} h\f$ containing the gradient of some scalar function \f$h\f$ with respect to output \f$i\f$ of this ModPiece.
      @return An Eigen::VectorXd containing the new gradient \f$\nabla_{x_j} h\f$.
    */
    virtual Eigen::VectorXd const& Gradient(unsigned int                 const  outputDimWrt,
                                            unsigned int                 const  inputDimWrt,
                                            std::vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd              const& sensitivity);

    virtual Eigen::VectorXd const& Gradient(unsigned int                const  outputDimWrt,
                                            unsigned int                const  inputDimWrt,
                                            ref_vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd             const& sensitivity);

    inline Eigen::VectorXd const& Gradient(unsigned int outWrt,
                                           unsigned int inWrt,
                                           Eigen::VectorXd const& last,
                                           Eigen::VectorXd const& sens) {
      ref_vector<Eigen::VectorXd> vec;
      vec.push_back(std::cref(last));
      return Gradient(outWrt, inWrt, vec, sens);
    }

    template<typename... Args>
    inline Eigen::VectorXd  const& Gradient(unsigned int wrtOut,
                                            unsigned int wrtIn,
                                            Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return GradientRecurse(wrtOut, wrtIn, vec, args...);
    }


    /** @brief Compute the Jacobian of this ModPiece.
      @details This function computes the Jacobian matrix of some output \f$f_i(x_1, \ldots, x_{M_x})\f$
               with respect to one of the inputs \f$x_j\f$.

      Like the `Evaluate` function, this function performs some checks on the inputs
      and then calls the `JacobianImpl` function.  It also keeps track of how many
      times `Jacobian` is called and the runtime of `JacobianImpl`.

      If a child class does not override the `JacobianImpl` function, finite
      differences are used by call the `JacobianByFD` function.

      @param[in] outputDimWrt The index \f$i\f$ describing the output of interest. (Think of the derivative \f$\partial f_i / \partial x_j\f$)
      @param[in] inputDimWrt The index \f$j\f$ describing which input we want to take derivatives with respect to.
      @param[in] input A vector of inputs corresponding to \f$x_1,\ldots, x_{M_x}\f$
      @return An Eigen::MatrixXd containing the new Jacobian matrix \f$ J_{ij}\f$.
    */
    virtual Eigen::MatrixXd const& Jacobian(unsigned int                 const  outputDimWrt,
                                            unsigned int                 const  inputDimWrt,
                                            std::vector<Eigen::VectorXd> const& input);

    virtual Eigen::MatrixXd const& Jacobian(unsigned int                const  outputDimWrt,
                                            unsigned int                const  inputDimWrt,
                                            ref_vector<Eigen::VectorXd> const& input);


    template<typename... Args>
    inline Eigen::MatrixXd const& Jacobian(unsigned int outWrt, unsigned int inWrt, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return Jacobian(outWrt, inWrt, vec, args...);
    }

    template<typename... Args>
    inline Eigen::MatrixXd JacobianByFD(unsigned int outWrt, unsigned int inWrt, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return JacobianByFD(outWrt, inWrt, vec, args...);
    }

    template<typename... Args>
    inline Eigen::MatrixXd ApplyJacobianByFD(unsigned int outWrt, unsigned int inWrt, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return ApplyJacobianByFD(outWrt, inWrt, vec, args...);
    }


    /** @brief Apply the Jacobian of this ModPiece to a vector.
      @details This function computes the matrix-vector product \f$J_{ij} v \f$,
       where \f$J_{ij}\f$ is the Jacobian matrix of output \f$f_i\f$ with respect
       to input \f$x_j\f$.

      Like the `Evaluate` function, this function performs some checks on the inputs
      and then calls the `ApplyJacobianImpl` function.  It also keeps track of how many
      times `ApplyJacobian` is called and the runtime of `ApplyJacobianImpl`.

      If a child class does not override the `ApplyJacobianImpl` function, finite
      differences are used by call the `ApplyJacobianByFD` function.  Computing
      \f$J_{ij} v \f$ with finite differences requires two evaluations regardless
      of the dimension of \f$v\f$.

      @param[in] outputDimWrt The index \f$i\f$ describing the output of interest. (Think of the derivative \f$\partial f_i / \partial x_j\f$)
      @param[in] inputDimWrt The index \f$j\f$ describing which input we want to take derivatives with respect to.
      @param[in] input A vector of inputs corresponding to \f$x_1,\ldots, x_{M_x}\f$
      @param[in] vec The vector \f$v\f$ that we want to multiply with the Jacobian matrix.
      @return An Eigen::VectorXd containing the matrix-vector product \f$ J_{ij} v\f$.
    */
    virtual Eigen::VectorXd const& ApplyJacobian(unsigned int                 const  outputDimWrt,
                                                 unsigned int                 const  inputDimWrt,
                                                 std::vector<Eigen::VectorXd> const& input,
                                                 Eigen::VectorXd              const& vec);

    virtual Eigen::VectorXd const& ApplyJacobian(unsigned int                const  outputDimWrt,
                                                 unsigned int                const  inputDimWrt,
                                                 ref_vector<Eigen::VectorXd> const& input,
                                                 Eigen::VectorXd             const& vec);

    inline Eigen::VectorXd const& ApplyJacobian(unsigned int outWrt, unsigned int inWrt, Eigen::VectorXd const& last, Eigen::VectorXd const& sens) {     \
      ref_vector<Eigen::VectorXd> vec;
      vec.push_back(std::cref(last));                                                                                               \
      return ApplyJacobian(outWrt, inWrt, vec, sens);                                                                              \
    }
    template<typename... Args>
    inline Eigen::VectorXd  const& ApplyJacobian(unsigned int wrtOut, unsigned int wrtIn, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return ApplyJacobianRecurse(wrtOut, wrtIn, vec, args...);
    }

    // virtual Eigen::VectorXd ApplyHessian(int                         const  outputDimWrt,
    //                                      int                         const  inputDimWrt1,
    //                                      int                         const  inputDimWrt2,
    //                                      ref_vector<Eigen::VectorXd> const& input,
    //                                      Eigen::VectorXd             const& sensitivity
    //                                      Eigen::VectorXd             const& vec);



    virtual Eigen::VectorXd GradientByFD(unsigned int                 const  outputDimWrt,
                                         unsigned int                 const  inputDimWrt,
                                         std::vector<Eigen::VectorXd> const& input,
                                         Eigen::VectorXd              const& sensitivity);

    virtual Eigen::VectorXd GradientByFD(unsigned int                const  outputDimWrt,
                                         unsigned int                const  inputDimWrt,
                                         ref_vector<Eigen::VectorXd> const& input,
                                         Eigen::VectorXd             const& sensitivity);

    virtual Eigen::MatrixXd JacobianByFD(unsigned int                 const  outputDimWrt,
                                         unsigned int                 const  inputDimWrt,
                                         std::vector<Eigen::VectorXd> const& input);

    virtual Eigen::MatrixXd JacobianByFD(unsigned int                const  outputDimWrt,
                                         unsigned int                const  inputDimWrt,
                                         ref_vector<Eigen::VectorXd> const& input);

    virtual Eigen::VectorXd ApplyJacobianByFD(unsigned int                 const  outputDimWrt,
                                              unsigned int                 const  inputDimWrt,
                                              std::vector<Eigen::VectorXd> const& input,
                                              Eigen::VectorXd              const& vec);

    virtual Eigen::VectorXd ApplyJacobianByFD(unsigned int                const  outputDimWrt,
                                              unsigned int                const  inputDimWrt,
                                              ref_vector<Eigen::VectorXd> const& input,
                                              Eigen::VectorXd             const& vec);

    const Eigen::VectorXi inputSizes;
    const Eigen::VectorXi outputSizes;

  protected:

    std::vector<Eigen::VectorXd> ToStdVec(ref_vector<Eigen::VectorXd> const& input);

    // The following variables keep track of how many times the Implemented functions, i.e. EvaluateImpl, GradientImpl,
    // etc... are called.
    unsigned long int numGradCalls   = 0;
    unsigned long int numJacCalls    = 0;
    unsigned long int numJacActCalls = 0;
    unsigned long int numHessCalls   = 0;

    // these variables keep track of the total wall-clock time spent in each of the Implemented functions.  They are in
    // units of milliseconds
    double gradTime   = 0;
    double jacTime    = 0;
    double jacActTime = 0;
    double hessTime   = 0;

    std::vector<Eigen::VectorXd> outputs;
    Eigen::VectorXd gradient;
    Eigen::VectorXd jacobianAction;
    Eigen::MatrixXd jacobian;

    void CheckInputs(ref_vector<Eigen::VectorXd> const& input);

    virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

    virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) = 0;

    virtual void GradientImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input,
                              Eigen::VectorXd             const& sensitivity);

    virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input);

    virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& vec);

    // virtual void ApplyHessianImpl(int                         const  outputDimWrt,
    //                       int                         const  inputDimWrt1,
    //                       int                         const  inputDimWrt2,
    //                       ref_vector<Eigen::VectorXd> const& input,
    //                       Eigen::VectorXd             const& sensitivity
    //                       Eigen::VectorXd             const& vec);


    // virtual void ApplyHessianByFD(int                         const  outputDimWrt,
    //                       int                         const  inputDimWrt1,
    //                       int                         const  inputDimWrt2,
    //                       ref_vector<Eigen::VectorXd> const& input,
    //                       Eigen::VectorXd             const& sensitivity
    //                       Eigen::VectorXd             const& vec);

  private:
    template<typename NextType, typename... Args>
    inline Eigen::VectorXd const& GradientRecurse(unsigned int outWrt,
                                                  unsigned int inWrt,
                                                  ref_vector<Eigen::VectorXd>& vec,
                                                  NextType const& ith,
                                                  Args const&... args) {

      static_assert(std::is_same<Eigen::VectorXd,
                    NextType>::value,
                    "In ModPiece::Gradient, cannot cast input to Eigen::VectorXd.");

        vec.push_back(std::cref((NextType&)ith));
        return GradientRecurse(outWrt, inWrt, vec, args...);

    }

    template<typename NextType>
    inline Eigen::VectorXd const& GradientRecurse(unsigned int outWrt,
                                                  unsigned int inWrt,
                                                  ref_vector<Eigen::VectorXd>& vec,
                                                  NextType const& last,
                                                  Eigen::VectorXd const& sens) {
        static_assert(std::is_same<Eigen::VectorXd, NextType>::value, "In ModPiece::Gradient, cannot cast input to Eigen::VectorXd.");
        vec.push_back(std::cref((NextType&)last));
        return Gradient(outWrt, inWrt, vec, sens);
    }

    template<typename... Args>
    inline Eigen::MatrixXd const& Jacobian(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& ith, Args const&... args) {     \
        vec.push_back(std::cref(ith));                                                                                               \
        return Jacobian(outWrt, inWrt, vec, args...);                                                                              \
    }
    inline Eigen::MatrixXd const& Jacobian(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& last) {
      vec.push_back(std::cref(last));
      return Jacobian(outWrt, inWrt, vec);
    }

    template<typename NextType, typename... Args>                                                                                               \
    inline Eigen::VectorXd const& ApplyJacobianRecurse(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& ith, Args const&... args) {     \
        static_assert(std::is_same<Eigen::VectorXd, NextType>::value, "In ModPiece::Gradient, cannot cast input to Eigen::VectorXd."); \
        vec.push_back(std::cref((NextType&)ith));                                                                                               \
        return ApplyJacobianRecurse(outWrt, inWrt, vec, args...);                                                                              \
    }

    template<typename NextType>                                                                                               \
    inline Eigen::VectorXd const& ApplyJacobianRecurse(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& last, Eigen::VectorXd const& sens) {     \
        static_assert(std::is_same<Eigen::VectorXd, NextType>::value, "In ModPiece::Gradient, cannot cast input to Eigen::VectorXd."); \
        vec.push_back(std::cref((NextType&)last));                                                                                               \
        return ApplyJacobian(outWrt, inWrt, vec, sens);                                                                              \
    }

    template<typename... Args>                                                                                               \
    inline Eigen::MatrixXd JacobianByFD(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& ith, Args const&... args) {     \
        vec.push_back(std::cref(ith));                                                                                               \
        return JacobianByFD(outWrt, inWrt, vec, args...);                                                                              \
    }
    inline Eigen::MatrixXd JacobianByFD(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& last) {
      vec.push_back(std::cref(last));
      return JacobianByFD(outWrt, inWrt, vec);
    }

    template<typename NextType, typename... Args>
    inline Eigen::MatrixXd ApplyJacobianByFD(unsigned int outWrt, unsigned int  inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& ith, Args const&... args) {
      vec.push_back(std::cref(ith));
      return ApplyJacobianByFD(outWrt, inWrt, vec, args...);
    }
    template<typename NextType>
    inline Eigen::MatrixXd ApplyJacobianByFD(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& last, Eigen::VectorXd const& sens) {
      vec.push_back(std::cref(last));
      return ApplyJacobianByFD(outWrt, inWrt, vec, sens);
    }

  };

  }
}



#endif // #ifndef MODPIECE_H
