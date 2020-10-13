#ifndef COVARIANCEKERNELS_H_
#define COVARIANCEKERNELS_H_

#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ConstantKernel.h"
#include "MUQ/Approximation/GaussianProcesses/LinearTransformKernel.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ProductKernel.h"
#include "MUQ/Approximation/GaussianProcesses/SquaredExpKernel.h"
#include "MUQ/Approximation/GaussianProcesses/SumKernel.h"
#include "MUQ/Approximation/GaussianProcesses/WhiteNoiseKernel.h"
#include "MUQ/Approximation/GaussianProcesses/MaternKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"

namespace muq
{
namespace Approximation
{

/**
\defgroup CovarianceKernels Covariance Kernels
\ingroup GaussianProcesses
\brief Covariance kernels for defining Gaussian Processes

Gaussian processes are defined by their mean function and covariance kernel.  In this group, we provide tools for constructing and using covariance kernels.  Several simple kernels are provided as well as tools for combining simple kernels into more complicated kernels.

Templates are used extensively in this module to maximize performance (by reducing virtual function calls).  Thus, using the c++11 "auto" keyword can result in much cleaner and readable code.

<h2>Combining Kernels</h2>
As an example of constructing a covariance kernel that is constructed from several simple kernels, consider a kernel of the form
\f[
k(x_1,x_2) = k_1(x_1,x_2)\, k_2(x_1,x_2) + k_3(x_1,x_2),
\f]
where \f$k_1\f$ is a squared exponential kernel, \f$k_2\f$ is a periodic kernel and \f$k_3\f$ is a white noise kernel.  Assume the dimension of \f$x_1\f$ and \f$x_2\f$ is stored in a variable called <code>dim</code>.  Then, the combined kernel \f$k(x,y)\f$ can be constructed with the following snippet:
\code{.cpp}
#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"

// Other setup ...

auto k1 = SquaredExpKernel(dim, var1, length1);
auto k2 = PeriodicKernel(dim, var2, length2, period);
auto k3 = WhiteNoiseKernel(dim, var3);

auto k = kernel1*kernel2 + kernel3;
\endcode
or, more succinctly, as
\code{.cpp}
auto k = SquaredExpKernel(dim, var1, length1) * PeriodicKernel(dim, var2, length2, period) + WhiteNoiseKernel(dim, var3);
\endcode

In either case, the <code>k</code> variable is an instance of "SumKernel<ProductKernel<SquaredExpKernel,PeriodicKernel>, WhiteNoiseKernel>".  The "auto" keyword allows us to avoid typing this long type.

<h2>Anisotropic Kernels</h2>
In many cases, the correlation of a GP in one dimension will be different that the correlation is some other direction.  For example, let's say we have two spatial variables \f$x\f$ and \f$y\f$ as well as a kernel of the form
\f[
k([x_1,y_1], [x_2, y_2]) = k_x(x_1, x_2)\, k_y(y_1, y_2).
\f]
Such kernels commonly arise when modeling anisotropic media (e.g., hydraulic conductivity fields).  In MUQ, it is possible to specify the dimensions that are used by a kernel.  For example, if \f$k_1\f$ and \f$k_2\f$ were both squared exponential kernels, than \f$k\f$ could be defined as
\code{.cpp}

// Only keep the 0 index
std::vector<unsigned> indsx = {0};
auto kx = SquaredExpKernel(2, indsx, varx, Lx);

// Only keep the 1 index
std::vector<unsigned> inds2 = {1};
auto ky = SquaredExpKernel(2, indsy, vary, Ly);

auto k = k1 * k2;

\endcode

It is also possible to define more complicated relationships.  For example, consider a third component \f$z\f$, and let the kernel be defined as
\f[
k([x_1,y_1,z_1], [x_2, y_2, z_2]) = k_{xy}([x_1,y_1], [x_2,y_2])\, k_z(z_1, z_2).
\f]
This kernel might be constructed with
\code{.cpp}

// Keep both the 0 and 1 indices
std::vector<unsigned> indsxy = {0,1};
auto kx = SquaredExpKernel(3, indsxy, varxy, Lxy);

// Only keep the 2 index
std::vector<unsigned> inds2 = {2};
auto kz = SquaredExpKernel(3, indsz, varz, Lz);

auto k = kxy * kz;

\endcode

*/


} // namespace Approximation
} // namespace muq



#endif // #ifndef COVARIANCEKERNELS_H_
