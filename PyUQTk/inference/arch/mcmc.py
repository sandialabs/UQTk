#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
import numpy as npy
import scipy.stats
import scipy.linalg
import math
import matplotlib.pyplot as plt

global Rmat,invRmat


#---------------------------------------------------------------------------------------
#  Simple Hamiltonian MCMC routine
#  Uses Leapfrog for the time stepping
#---------------------------------------------------------------------------------------
def HMCMC(U,grad_U,dt,nT,q):
    '''
    Hamiltonian MCMC routine
    
    Input:
    -----
    
    U      - potential energy function, -log(posterior)
    grad_U - gradient of potential energy function
    dt     - time step, dt, for leapfrog method
    nT     - number of time steps in leapfrog method
    q      - initial state of chain (position vector)
    
    Output:
    ------
    Next vector in the chain state

    Example:
    -------
    q_next = HMCMC(U,grad_U,1e-2,25,q_current)
    
    '''
    current_q = npy.copy(q) # save current

    # generate current p
    # propcov = 4*array([[ 0.01175383,  0.02065261],[ 0.02065261,  0.04296117]])
    p = npy.random.randn(len(current_q))
    # p = random.multivariate_normal([0,0],propcov)
    current_p = npy.copy(p) # save current p

    # make half step for momentum used for leap frog step
    p = p - dt * grad_U(q)/2.0

    for i in range(nT):
        # p = p - dt * grad_U(q)/2.0
        q = q + dt*p
        # p = p - dt * grad_U(q)/2.0
        if (i != nT-1): p = p - dt*grad_U(q)

    # make a half step for momentum at the end
    p = p - dt * grad_U(q)/2.0

    # negate the momentum to make a symmetric proposal
    p = -p

    # Evaluate potential and kinetic energy
    current_U = U(current_q)[0]
    current_K = npy.sum(current_p**2)/2.0
    proposed_U = U(q)[0]
    proposed_K = npy.sum(p**2)/2.0

    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position

    if (npy.log(npy.random.rand()) < current_U-proposed_U+current_K-proposed_K):
        return q
    else:
        return current_q


#---------------------------------------------------------------------------------------
#  Example:
#  1. Banana-shaped posterior density
#---------------------------------------------------------------------------------------
def norm_pdf_multivariate(x, mu, sigma):
    """
    Multi-variate normal pdf
    x    : list or numpy array
    mu   : 1D numpy array 
    sigma: 2D numpy array"""
    size = len(x)
    if size == len(mu) and (size, size) == sigma.shape:
        det = npy.linalg.det(sigma)
        if det == 0:
            raise NameError("The covariance matrix can't be singular")
        norm_const = 1.0/ ( math.pow((2*npy.pi),float(size)/2) * math.pow(det,1.0/2) )
        x_mu = npy.matrix(x - mu)
        inv = npy.linalg.inv(sigma)
        result = math.pow(math.e, -0.5 * (x_mu * inv * x_mu.T))
        return norm_const * result
    else:
        raise NameError("The dimensions of the input don't match")

def tranB(x1,x2,a):
    """
    Coordinate transform for banana-shaped pdf
    x1,x2: 2D numpy arrays
    a: list containing the transform factors
    """
    a1 = a[0]; a2 = a[1]; 
    y1 = a1*x1;
    y2 = x2/a1 - a2*(y1**2 + a1**2);
    return y1,y2

def invTranB(x1,x2,a):
    """ Inverse coordinate transform for banana-shaped pdf
    x1,x2: 2D numpy arrays
    a: list containing the transform factors
    """
    a1 = a[0]; a2 = a[1]; 
    y1 = x1/a1;
    y2 = x2*a1 + a1*a2*(x1**2 + a1**2);
    return y1,y2

def plotBanana():
    """
    Plot banana-shaped function; parameters are hard-wired
    """
    xb,yb = npy.mgrid[-3:3:.05, -11:1:.05]
    x, y  = invTranB(xb,yb,[1,1])
    pos = npy.empty(x.shape + (2,))
    pos[:, :, 0] = x; pos[:, :, 1] = y
    mu  = npy.array([0.0,0.0])
    cov = npy.array([[1.0, 0.9], [0.9, 1.0]])
    z = x.copy()
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            z[i,j] = norm_pdf_multivariate([x[i,j],y[i,j]], mu, cov)
    plt.contour(xb,yb,z,50)
    plt.show()
    return 

def postBanana(spl,postinfo):
    """
    Posterior density for banana-shaped pdf
    postinfo : setup for the posterior density
    
    """
    afac = postinfo['afac']
    mu   = postinfo['mu'  ]
    cov  = postinfo['cov' ]
    xb,yb = spl ;
    x, y  = invTranB(xb,yb,afac) ;
    return npy.log(norm_pdf_multivariate([x,y], mu, cov))

#---------------------------------------------------------------------------------------
#  DRAM
#---------------------------------------------------------------------------------------
def logPropRatio(iq,spls):
    """
    Gaussian n:th stage log proposal ratio
    log of q_i(y_n,..,y_n-j) / q_i(x,y_1,...,y_j)
    """
    global invRmat
    stage = len(spls)-1;
    if stage == iq:
        return (0.0); # symmetric
    else:
        iRmat = invRmat[iq-1]; # proposal^(-1/2)
        y1 = spls[0]        ; # y1
        y2 = spls[iq]       ; # y_i
        y3 = spls[stage   ] ; # y_n
        y4 = spls[stage-iq] ;  # y_(n-i)
        return (-0.5*(npy.linalg.norm(npy.dot(y4-y3,iRmat))**2-npy.linalg.norm(npy.dot(y2-y1,iRmat))**2));

def logPostRatio(p1,p2):
    return (p2-p1);

def getAlpha(spls,post):
    stage = len(spls) - 1; 
    a1 = 1.0; a2 = 1.0;
    for k in range(1,stage):
        a1 = a1*(1-getAlpha(spls[:k+1],post[:k+1]));
        a2 = a2*(1-getAlpha(spls[-1:-(k+1):-1],post[-1:-(k+1):-1]));
        if  a2 == 0.0:
            return (0.0);
    y = logPostRatio(post[0],post[-1]);
    for k in range(1,stage+1):
        y = y + logPropRatio(k,spls);
    return min(1.0, npy.exp(y)*a2/a1);

def ucov(spl,splmean,cov,lastup):
    #
    #  update covariance
    #
    if len(spl.shape) == 1:
        nspl = 1;
        ndim = spl.shape[0];
    else:
        (nspl,ndim)=spl.shape;
    if nspl>0:
        for i in range(nspl):
            iglb   = lastup+i;
            splmean = (iglb*splmean+spl[i])/(iglb+1);
            rt = (iglb-1.0)/iglb;
            st = (iglb+1.0)/iglb**2;
            cov = rt*cov+st*npy.dot(npy.reshape(spl[i]-splmean,(ndim,1)),npy.reshape(spl[i]-splmean,(1,ndim)))
    return lastup+nspl,splmean,cov

def dram_ex(method,nsteps):
    # define MCMC parameters
    cini  = npy.array([-1.0,-4.0])
    spllo = npy.array([-4.0,-12.0])
    splhi = npy.array([ 4.0,  2.0])
    cvini = npy.array([[0.1,0.0],[0.0,0.1]])
    opts={'method':method,'nsteps':nsteps,'nburn':1000,'nadapt':100,'nfinal':10000000,
          'inicov':cvini,'coveps':1.e-10,'burnsc':5,'ndr':2,'drscale':[5,4,3],
          'spllo':spllo,'splhi':splhi}
    lpinfo={'afac':[1.0,1.0],'cov': npy.array([[1,0.9],[0.9,1]]),'mu':npy.array([0.0,0.0])}
    sol=dram(opts,cini,postBanana,lpinfo)
    return sol

def dram(opts,cini,likTpr,lpinfo):
    """
    #
    # DRAM
    #
    Delayed Rejection Adaptive MCMC
    opts - dictionary of parameters for DRAM
           method : either 'am' (adaptive metropolis) or 'dram' (am+delayed rejection)
           nsteps : no. of mcmc steps
           nburn  : no. of mcmc steps for burn-in (proposal fixed to initial covariance)
           nadapt : adapt every nadapt steps after nburn
           nfinal : stop adapting after nfinal steps
           inicov : initial covariance
           coveps : small additive factor to ensure covariance matrix is positive definite
           burnsc : factor to scale up/down proposal is acceptance rate is too high/low
           ndr    : no. of delayed rejection steps (if dram is requested)
           drscale: scale factors for delayed rejection
    cini - starting mcmc state
    likTpr - log-posterior function
    lpinfo - dictionary with settings that will be passed to the log-posterior function
    
    """    
    # -------------------------------------------------------------------------------
    # Parse options
    # -------------------------------------------------------------------------------
    if 'method' in opts:
        method = opts['method']
    else:
        print 'Error in dram: method unspecified !'; quit()
    nsteps = opts['nsteps']
    nburn  = opts['nburn' ]
    nadapt = opts['nadapt']
    nfinal = opts['nfinal']
    inicov = opts['inicov']
    coveps = opts['coveps']
    burnsc = opts['burnsc']
    spllo  = opts['spllo' ]
    splhi  = opts['splhi' ]
    if method=='dram':
        ndr     = opts['ndr']
        drscale = opts['drscale']
    rej    = 0;
    rejlim = 0;
    rejsc  = 0;
    # -------------------------------------------------------------------------------
    # Pre-processing
    # -------------------------------------------------------------------------------
    cdim   = cini.shape[0];            # chain dimensionality
    cov    = npy.zeros((cdim,cdim));   # covariance matrix
    spls   = npy.zeros((nsteps,cdim)); # MCMC samples
    na     = 0;                        # counter for accepected jumps
    sigcv  = 2.4/npy.sqrt(cdim);       # covariance factor
    spls[0] = cini;                    # initial sample set
    p1 = likTpr(spls[0],lpinfo);       # and 
    pmode = p1;                        # store chain MAP
    cmode = spls[0];
    nref  = 0;
    for k in range(nsteps-1):
        #
        # Deal with covariance matrix
        #
        covMatUpd = False
        if k == 0:
            splmean   = spls[0];
            propcov   = inicov ;
            Rchol     = scipy.linalg.cholesky(propcov) ;
            lastup    = 1;      # last covariance update
            covMatUpd = True ;
        else:
            if (nadapt>0) and ((k+1)%nadapt)==0:
                if k<nburn:
                    if float(rejsc)/nref>0.95:
                        Rchol = Rchol/burnsc # scale down proposal
                        covMatUpd = True ;
                        print "Scaling down the proposal at step",k
                    elif float(rejsc)/nref<0.05:
                        Rchol = Rchol*burnsc # scale up proposal
                        covMatUpd = True ;
                        print "Scaling up the proposal at step",k
                    nref  = 0 ;
                    rejsc = 0 ;
                else:
                    lastup,splmean,cov=ucov(spls[lastup:lastup+nadapt,:],splmean,cov,lastup)
                    try:
                        Rchol = scipy.linalg.cholesky(cov)
                    except scipy.linalg.LinAlgError:
                        try:
                            # add to diagonal to make the matrix positive definite
                            Rchol = scipy.linalg.cholesky(cov+coveps*npy.identity(cdim))
                        except scipy.linalg.LinAlgError:
                            print "Covariance matrix is singular even after the correction"
                    Rchol = Rchol*sigcv
                    covMatUpd = True ;
        if (method == 'dram') and covMatUpd:
            Rmat = [Rchol]; invRmat = [scipy.linalg.inv(Rchol)]
            for i in range(1,ndr):
                Rmat.append(Rmat[i-1]/drscale[i-1])
                invRmat.append(invRmat[i-1]*drscale[i-1])
        #-Done with covariance matrix
        nref = nref + 1 ;
        #
        # generate proposal and check bounds
        #
        u  = spls[k]+npy.dot(npy.random.randn(1,cdim),Rchol)[0];
        if npy.any(npy.less(u,spllo)) or npy.any(npy.greater(u,splhi)):
            outofbound = True
            accept     = False
            p2 = -1.e6
        else:
            outofbound = False
        if not outofbound:
            p2 = likTpr(u,lpinfo);
            pr = npy.exp(p2-p1);
            if (pr>=1.0) or (npy.random.random_sample()<=pr):
                spls[k+1] = u.copy();
                p1 = p2;
                if p1 > pmode:
                    pmode = p1 ;
                    cmode = spls[k+1] ;
                accept = True
            else:
                accept = False
        #
        # See if we can do anything about a rejected proposal
        #
        if not accept:
            if (method == 'am'):
                # if 'am' then reject
                spls[k+1]=spls[k];
                rej   = rej   + 1;
                rejsc = rejsc + 1;
                if outofbound:
                    rejlim  = rejlim + 1;
            elif (method == 'dram'):
                # try delayed rejection
                tryspls = [spls[k].copy(),u.copy()]
                trypost = [p1,p2]
                jdr = 1
                while (not accept) and (jdr<ndr):
                    jdr = jdr+1;
                    u  = spls[k]+npy.dot(npy.random.randn(1,cdim),Rmat[jdr-1])[0];
                    if npy.any(npy.less(u,spllo)) or npy.any(npy.greater(u,splhi)):
                        outofbound = True
                        tryspls.append(u.copy())
                        trypost.append(-1.0e6)
                        continue
                    outofbound = False
                    p2 = likTpr(u,lpinfo);
                    tryspls.append(u.copy())
                    trypost.append(p2)
                    alpha = getAlpha(tryspls,trypost);
                    if (alpha >= 1.0) or (npy.random.random_sample() < alpha): 
                        accept = True;
                        spls[k+1] = u.copy();
                        p1 = p2;
                        if p1 > pmode:
                            pmode = p1 ;
                            cmode = spls[k+1] ;
                if not accept:
                    spls[k+1]=spls[k] ;
                    rej   = rej   + 1;
                    rejsc = rejsc + 1;
                    if outofbound:
                        rejlim  = rejlim + 1;
            else:
                print "Unknown MCMC method ",method," -> Quit\n"; quit()
            # Done with if over methods
        # Done with if over original accept
    # Done loop over all steps
    return (spls,[cmode,pmode],[1.0-float(rej)/nsteps,1.0-float(rejlim)/nsteps],[rej,rejlim])
