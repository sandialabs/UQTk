#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
try:
    import uqtkarray
    import quad as uqtkquad
    import pce as uqtkpce
    import tools as uqtktools
    import pce_tools
except ImportError:
    import PyUQTk.uqtkarray as uqtkarray
    import PyUQTk.quad as uqtkquad
    import PyUQTk.pce as uqtkpce
    import PyUQTk.tools as uqtktools
    from PyUQTk.PyPCE import pce_tools
except ImportError:
    print("PyUQTk array, quad, pce, tools or pce_tools modules not found")

try:
    import numpy as np
except ImportError:
    print('numpy needs to be installed')
#######################################################################

def gauss_adaptation(c_k, ndim, method = 0):
    '''
    Input:
           c_k: N0 dimensional numpy array of first order coefficients,
                where N0 is the # of pc terms of 1 dimensional expansion
          ndim: dimension of the problem
	    method: The method to compute the isometry (Rotation matrix), from set {0,1,2,3}. Returns isometry.
            0 : (default) By Gram-Schmidt procedure on matrix A with Gaussian coeffs (normalized) at its first row,
                and ones along diagonal zeros elsewhere for other rows.
            1 : By orthogonal decomposition of a * a.T, where a contains Gaussian coeffs
            2 : By orthogonal decomposition of the Householder matrix H = I - 2 a a.T / ||a|| ** 2
            3 : Via a Gram-Schmidt procedure on the matrix A with Gaussian coeffs (normalized) at its first row, and,
                starting from 2th row, put the ith largest Gaussian coeff on the column corresponding to its
                variable xi at (i+1) row zeros elsewhere.
    Output:
           ndim by ndim numpy array, the rotation matrix
    '''
    assert c_k.shape[0] == ndim
    a = c_k.reshape((c_k.shape[0], 1))
    if method == 0:   # Gram-Schmidt
        C = np.eye(a.shape[0])
        C[0,:] = a[:,0]
        [q,r] = np.linalg.qr(C.T)
        return q.T
    elif method == 1: # orthogonal decomposition of a * a.T
        C = np.dot(a, a.T)
        [vals, vecs] = np.linalg.eigh(C)
        return vecs[:,::-1].T
    elif method == 2: # Orthogonal decomposition of Householder matrix
        H = np.eye(ndim) - 2 * np.dot(a, a.T) / np.linalg.norm(a) ** 2 # Householder matrix
        [vals, vecs] = np.linalg.eigh(H)
        return np.real(vecs).T
    elif method == 3: # Sort by importance, recommended method
        c3 = np.argsort(np.abs(c_k))[::-1]
        C = np.zeros((ndim, ndim))
        C[:,0] = a[:,0]
        loc = 0
        for i in range(0, np.size(c_k)-1):
            C[c3[loc], i+1] = c_k[c3[loc]]
            loc += 1
        [q,r] = np.linalg.qr(C)
        q=np.array(q)
        return q.T
    else:
        raise ValueError('Method parameter must be in {0,1,2,3}')


def eta_to_xi_mapping(eta, A, zeta = None):
    '''
    Maps points from lower dimensional eta space to the xi space.
    A is isometry which serves as the rotation matrix (xi = A' [eta, zeta])
    Input:
        eta : N by d0 numpy array, eta space points, N can be # of quarature points or MC points
        A   : d by d numpy array, Rotation matrix or isometry
        zeta: (optional) N by d-d0 numpy array, matrix to expand eta to make eta N by d
    Output:
        N by d numpy array, xi's mapped from eta's
    '''

    assert A.shape[0] == A.shape[1]

    d0 = eta.shape[1]
    d = A.shape[0]
    N = eta.shape[0]
    if zeta == None:
        zeta = np.zeros((N, d-d0))
    else:
        assert eta.shape[0] == zeta.shape[0]
        assert eta.shape[1] + zeta.shape[1] == A.shape[0]
    eta_full = np.hstack([eta, zeta]) # Augment eta with zeros
    return np.dot(A.T, eta_full.T).T


def mi_terms_loc(d1, d2, nord, pc_type, param, sf, pc_alpha=0.0, pc_beta=1.0):
    '''
    Locate basis terms in Multi-index matrix
    Input:
        d1     : interger, lower dimension to be located
        d2     : interger, higher dimension to locate
        nord   : PC order
        pc_type: polynomial type
        param  : quadrature level
        sf     : "sparse" or "full" indicator
        pc_alpha, pc_beta: parameters
    Output:
        locs   : N1 numpy array (N1 is # of PC terms in d1 dimesnional expansion),
                 the coefficients locations of d1 dimensional expansion in d2 dimensional expansion
    '''
    assert d1 < d2

    # obtain pc_model of d2 dimensional expansion
    pc_model2 = uqtkpce.PCSet("NISP", nord, d2, pc_type, pc_alpha, pc_beta)
    if d2==1:
        pc_model2.SetQuadRule(pc_type, 'full', param)
    else:
        pc_model2.SetQuadRule(pc_type, sf, param)

    # obtain multi-indices of d2 dimensional expansion
    MI_uqtk = uqtkarray.intArray2D()
    pc_model2.GetMultiIndex(MI_uqtk)
    MI = np.zeros((pc_model2.GetNumberPCTerms(), d2),dtype='int64')
    MI_uqtk.getnpintArray(MI)

    # find locations where the first d1 multi-indices of d2 space agree with
    # multi-indices in d1 space, while the remaining indices equal zeros
    if d2 == d1 + 1:  # in this case, just find locations of zeros in the last entries
        return np.where(MI[:,-1] == 0)[0]
    else:
        TF = MI[:,d1:] == [0]         # find zeros at individual d2-d1 ~ d2 location and mark with bool value
        TF = np.all(TF,axis=1)        # find where locations where zeros appear at all d2-d1 ~ d2 locations
        locs = np.where(TF)[0]
        return locs


def l2_error_eta(c_1, c_2, d1, d2, nord, pc_type, param, sf, pc_alpha=0.0, pc_beta=1.0):
    '''
    l2-norm relative error function
    Return relative l2-error norm(c_1 - c_2) / norm(c_2)
    Inputs:
        c_1: N1 dimensional numpy array (N1 is # of PC terms in d1 dimesnional expansion),
             PC coefficients of lower dimension
        c_2: N2 dimensional numpy array (N2 is # of PC terms in d2 dimesnional expansion),
             PC coefficients of higher dimension
        other parameters are explained as above
    Output:
        relative l2-error norm(c_1 - c_2) / norm(c_2)
        C1 : N2 dimensional numpy array, projection of c_1 into d2 expansion space
    '''
    import math
    assert np.shape(c_2)[0] == math.factorial(d2+nord) / (math.factorial(nord) * math.factorial(d2))
    assert np.shape(c_1)[0] <= np.shape(c_2)[0]
    C1 = np.zeros(c_2.shape[0])
    # call mi_terms_loc to make projections
    C1[mi_terms_loc(d1, d2, nord, pc_type, param, sf, pc_alpha, pc_beta)] = c_1
    return (np.linalg.norm((C1 - c_2),2) / np.linalg.norm(c_2,2)), C1

def transf_coeffs_xi(coeffs, nord, ndim, pc_type, param, R, sf="sparse", pc_alpha=0.0, pc_beta=1.0):
    '''
    Transfer coefficient from eta-space to xi-space to check convergence
    eta-space is the current dimension-reduced space, while
    xi-space is the full dimension expansion space.
    Dimension of Eta is equal to dimension of Xi
    Inputs:
        coeffs: N dimensional numpy array (N is the # of PC terms of ndim dimensional expansion),
                coefficients in eta space
    Output:
        N dimensional numpy array, coefficients projected to xi space
    this function is only used to check the accuracy of adaptation method
    '''

    ##### obtain pc_model of eta space and xi space #####
    pc_model_eta = uqtkpce.PCSet("NISP", nord, ndim, pc_type, pc_alpha, pc_beta)
    pc_model_xi = uqtkpce.PCSet("NISP", nord, ndim, pc_type, pc_alpha, pc_beta)
    if ndim == 1:
        pc_model_eta.SetQuadRule(pc_type, 'full', param)
        pc_model_xi.SetQuadRule(pc_type, 'full', param)
    else:
        pc_model_eta.SetQuadRule(pc_type, sf, param)
        pc_model_xi.SetQuadRule(pc_type, sf, param)
    qdpts_eta, totquat_eta= pce_tools.UQTkGetQuadPoints(pc_model_eta)

    ##### Obtain Psi_xi at quadrature points of eta #####
    qdpts_xi = eta_to_xi_mapping(qdpts_eta, R)
    qdpts_xi_uqtk = uqtkarray.dblArray2D(totquat_eta, ndim)
    qdpts_xi_uqtk.setnpdblArray(np.asfortranarray(qdpts_xi))

    psi_xi_uqtk = uqtkarray.dblArray2D()
    pc_model_xi.EvalBasisAtCustPts(qdpts_xi_uqtk, psi_xi_uqtk)
    psi_xi = np.zeros((totquat_eta, pc_model_xi.GetNumberPCTerms()))
    psi_xi_uqtk.getnpdblArray(psi_xi)

    ##### Obtian Psi_eta at quadrature points of eta #####
    weight_eta_uqtk = uqtkarray.dblArray1D()
    pc_model_eta.GetQuadWeights(weight_eta_uqtk)
    weight_eta=np.zeros((totquat_eta,))
    weight_eta_uqtk.getnpdblArray(weight_eta)

    psi_eta_uqtk = uqtkarray.dblArray2D()
    pc_model_eta.GetPsi(psi_eta_uqtk)
    psi_eta = np.zeros( (totquat_eta, pc_model_eta.GetNumberPCTerms()) )
    psi_eta_uqtk.getnpdblArray(psi_eta)

    return np.dot(coeffs, np.dot(psi_eta.T * weight_eta, psi_xi))
