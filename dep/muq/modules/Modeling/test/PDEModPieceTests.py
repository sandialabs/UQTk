import sys
sys.path.insert(0, 'modules/Modeling/python')

# import unit testing framework
import unittest

# import numpy
import numpy as np

# import dolfin/fenics
import dolfin as dfn

import matplotlib.pyplot as plt

# import muq Modeling
import pymuqModeling as mm

class PDEModPieceTests(unittest.TestCase):
    # solve -k u_{xx} = 1, u(x=0) = 1 and -k u_{x}(x=1) = 1, with k = 0.5
    # the exact solution is u(x) = -x^2/(2k) + 1
    def test_Poisson(self):
        # create a 1D mesh
        Nelem = 150 # number of elements in the mesh
        mesh = dfn.UnitIntervalMesh(Nelem)

        # create function spaces
        Vh2 = dfn.FunctionSpace(mesh, 'Lagrange', 2) # quadratic
        Vh1 = dfn.FunctionSpace(mesh, 'Lagrange', 1) # linear
        #Vh = [Vh2, Vh1, Vh2] # function space for [state, parameter, adjoint]

        # define the weak form
        def weak_form(u, k, p):
            return k*dfn.inner(dfn.grad(u), dfn.grad(p))*dfn.dx - p*dfn.dx + p*dfn.ds

        # define the Dirichlet boundary
        def boundary(x, on_boundary):
            return on_boundary and x[0]<dfn.DOLFIN_EPS

        # define the boundary condition for the forward problem
        u0 = dfn.Constant(1.0)
        bc_fwd = dfn.DirichletBC(Vh2, u0, boundary)

        # define the boundary condition for the adjoint problem
        p0 = dfn.Constant(0.0)
        bc_adj = dfn.DirichletBC(Vh2, p0, boundary)

        # create the PDE model
        pde = mm.PDEModPiece(Vh2, Vh1, weak_form, bc_fwd, bc_adj, is_fwd_linear=True)

        # create a constant field for the parameter
        k = dfn.Function(Vh1)
        k.vector()[:] = 0.5

        # evaluate the pde model
        u = pde.Evaluate([k.vector()]) [0]

        # import the solution into a function
        soln = dfn.Function(Vh2)
        soln.vector().set_local(u)

        # check the numerical solution against the true solution
        for x in np.linspace(0.0, 1.0, num=50):
            self.assertAlmostEqual(soln([x]), 1.0-x*x)

        # create a sensitivity function
        sens = dfn.Function(Vh2)
        sens.vector()[:] = 1.0

        # comute the gradient with finite differences
        gradFD = pde.GradientByFD(0, 0, [k.vector()], sens.vector())

        # compute the gradient with the adjoint method
        grad = pde.Gradient(0, 0, [k.vector()], sens.vector())

        # check to make sure that the adjoint and FD gradients are similar
        self.assertEqual(len(grad), len(gradFD))
        for i in range(len(grad)):
            self.assertAlmostEqual(grad[i], gradFD[i], places=1)
