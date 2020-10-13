# make sure that the HDF5 library is available
set(CMAKE_REQUIRED_LIBRARIES ${NLOPT_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${NLOPT_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")

CHECK_CXX_SOURCE_COMPILES(
  "
#include <math.h>
#include <nlopt.h>

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }


int main()
{

double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
nlopt_opt opt;

opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
nlopt_set_lower_bounds(opt, lb);
nlopt_set_min_objective(opt, myfunc, NULL);

my_constraint_data data[2] = { {2,0}, {-1,1} };

nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);

nlopt_set_xtol_rel(opt, 1e-4);

double x[2] = { 1.234, 5.678 };  /* some initial guess */
double minf; /* the minimum objective value, upon return */

nlopt_optimize(opt, x, &minf);

return 0;
}  
  "
  NLOPT_COMPILES)
  