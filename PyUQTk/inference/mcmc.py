#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
import numpy as npy
import scipy.stats
import scipy.linalg
import math
import uuid
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
    Computes the Log of the posterior density for banana-shaped pdf

    Input:
        spl:        Current parameter set sample
        postinfo :  Contains parameters for the posterior density
    Output:
        The log of the posterior density
    """

    afac = postinfo['afac']
    mu   = postinfo['mu'  ]
    cov  = postinfo['cov' ]
    xb,yb = spl ;
    x, y  = invTranB(xb,yb,afac) ;
    return npy.log(norm_pdf_multivariate([x,y], mu, cov))

def dram_ex(method,nsteps):
    """
    Example using the DRAM sampler to explore the posterior of the banana-shaped
    posterior density.

    Input:
        method: either 'am' or 'dram' (see below under the dram function)
        nsteps: number of steps to take (samples to take)
    Output:
        A tuple with samples and other information. See the dram function for more info.
    """
    # define MCMC parameters (see below under the dram function for more info about these options)
    cini  = npy.array([-1.0,-4.0])      # Initial guesses
    spllo = npy.array([-4.0,-12.0])     # Lower bounds on samples
    splhi = npy.array([ 4.0,  2.0])     # Upper bounda on samples
    cvini = npy.array([[0.1,0.0],[0.0,0.1]])    # Initial covariance matrix of proposal distribution
    opts={'method':method,'nsteps':nsteps,'nburn':1000,'nadapt':100,'nfinal':10000000,
          'inicov':cvini,'coveps':1.e-10,'burnsc':5,'gamma':1.0,'ndr':2,'drscale':[5,4,3],
          'spllo':spllo,'splhi':splhi}
    lpinfo={'afac':[1.0,1.0],'cov': npy.array([[1,0.9],[0.9,1]]),'mu':npy.array([0.0,0.0])}
    sol=dram(opts,cini,postBanana,lpinfo)
    return sol
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
                    (only added to diagonal if covariance matrix is singular without it)
           burnsc : factor to scale up/down proposal if acceptance rate is too high/low
           gamma  : factor to multiply proposed jump size with in the chain past the burn-in phase
                    (Reduce this factor to get a higher acceptance rate.)
                    (Defaults to 1.0)
           ndr    : no. of delayed rejection steps (if dram is requested)
           drscale: scale factors for delayed rejection
           spllo  : lower bounds for chain samples
           splhi  : upper bounds for chain samples
           rnseed : Optional seed for random number generator (needs to be integer >= 0)
                    If not specified, then random number seed is not fixed and every chain will
                    be different.
           tmpchn : Optional;
                           if present, will save chain state every 'ofreq' to ascii file.
                           filaname is randomly generated if tmpchn is set to 'tmpchn', or set to
                           the string passed through this option
                           if not present, chain states are not saved during the MCMC progress
    cini    - starting mcmc state
    likTpr  - log-posterior function; it takes two input parameters as follows
                - first parameter is a 1D array containing the chain state at which the posterior
                will to be evaluated
                - the second parameter contains settings the user can pass to this function;
                see below info for 'lpinfo'
            - this function is expected to return log-Likelihood and log-Prior values (in this order)

    lpinfo - object containing settings that will be passed to the log-posterior function;
                this object can be of any type (e.g. None, scalar, list, array, dictionary, etc)
                as long as it is consistent with settings expected inside the 'likTpr' function

    Output:
      spls: chain samples (dimension nsteps x chain dimension)
      [cmode,pmode]: MAP estimate (cmode) and posterior at MAP estimate (pmode)
      [1.0-float(rej)/nsteps,
       1.0-float(rejlim)/nsteps]: acceptance ratio and fraction of samples inside the bounds
      [rej,rejlim]: total number of rejected samples and total number
                    of samples outside the bounds
      meta_info: acceptance probability and posterior probability for each sample (dimension nsteps x 2)

    To Do:
      Provide option to dump MCMC chain as the computations proceed, to avoid having such large
      files to hold all states, and so that partial output is available during the MCMC run for
      preliminary analysis.
    """
    # -------------------------------------------------------------------------------
    # Parse options
    # -------------------------------------------------------------------------------
    if 'method' in opts:
        method = opts['method']
    else:
        print('Error in dram: method unspecified !')
        return {}

    nsteps = opts['nsteps']
    nburn  = opts['nburn' ]
    nadapt = opts['nadapt']
    nfinal = opts['nfinal']
    inicov = opts['inicov']
    coveps = opts['coveps']
    burnsc = opts['burnsc']
    spllo  = opts['spllo' ]
    splhi  = opts['splhi' ]

    if 'gamma' not in opts:
        gamma = 1.0 # Default for backwards compatibility
    else:
        gamma  = opts['gamma' ]

    if method=='dram':
        ndr     = opts['ndr']
        drscale = opts['drscale']

    if 'ofreq' not in opts:
        ofreq = 10000 # Default for backwards compatibility
    else:
        ofreq  = opts['ofreq' ]

    if 'tmpchn' not in opts:
        tmp_file = 'None'
    else:
        if opts['tmpchn'] == 'tmpchn':
            tmp_file = str(uuid.uuid4())+'.dat'
        else:
            tmp_file = opts['tmpchn']
        print('Saving intermediate chains to', tmp_file)

    # If desired, fix random number seed to make chain reproducible
    if 'rnseed' in opts:
        iseed = opts['rnseed']
        if isinstance(iseed, (int,long)) and iseed >= 0:
            npy.random.seed(iseed)
            print('\nmcmc::dram Fixing the random number seed to ', iseed)
        else:
            print('\nWARNING: mcmc::dram invalid random number seed specified: ', iseed)
            print('Will proceed without fixing random number seed.\n')

    rej    = 0;                        # Counts number of samples rejected
    rejlim = 0;                        # Counts number of samples rejected as out of prior bounds
    rejsc  = 0;                        # Counts number of rejected samples since last rescaling
    # -------------------------------------------------------------------------------
    # Pre-processing
    # -------------------------------------------------------------------------------
    cdim   = cini.shape[0]             # chain dimensionality
    cov    = npy.zeros((cdim,cdim))    # covariance matrix
    spls   = npy.zeros((nsteps,cdim))  # MCMC samples
    meta_info = npy.zeros((nsteps,3))  # Column for acceptance probability and posterior prob. of current sample
    na     = 0                         # counter for accepted jumps
    sigcv  = 2.4*gamma/npy.sqrt(cdim)  # covariance factor
    spls[0] = cini                     # initial sample set
    p1LikPri = likTpr(spls[0],lpinfo)  # and posterior probability of initial sample set
    if not isinstance(p1LikPri, list):
        print('\nERROR: This version requires the model return both log-likelihood and log-prior')
        return {}
    p1 = p1LikPri[0]+p1LikPri[1]
    meta_info[0] = [0.e0,p1LikPri[0],p1LikPri[1]]     # Arbitrary initial acceptance and posterior probability of initial guess
    pmode = p1                         # store current chain MAP probability value
    cmode = spls[0]                    # current MAP parameter Set
    nref  = 0                          # Samples since last proposal rescaling
    # -------------------------------------------------------------------------------
    # Main loop
    # -------------------------------------------------------------------------------
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
                        print('Scaling down the proposal at step %d'%(k))
                    elif float(rejsc)/nref<0.05:
                        Rchol = Rchol*burnsc # scale up proposal
                        covMatUpd = True ;
                        print('Scaling up the proposal at step %d'%(k))
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
                            print('WARNING: Covariance matrix is singular even after the correction')
                    Rchol = Rchol*sigcv
                    covMatUpd = True ;
        if (method == 'dram') and covMatUpd:
            Rmat = [Rchol]; invRmat = [scipy.linalg.inv(Rchol)]
            for i in range(1,ndr):
                Rmat.append(Rmat[i-1]/drscale[i-1])
                invRmat.append(invRmat[i-1]*drscale[i-1])
        #-Done with covariance matrix
        nref = nref + 1
        #
        # generate proposal and check bounds
        #
        u  = spls[k]+npy.dot(npy.random.randn(1,cdim),Rchol)[0];
        if npy.any(npy.less(u,spllo)) or npy.any(npy.greater(u,splhi)):
            outofbound = True
            accept     = False
            p2 = -1.e100        # Arbitrarily low posterior likelihood
            pr = -1.e100        # Arbitrarily low acceptance probability
        else:
            outofbound = False
        if not outofbound:
            p2Lik,p2Pri = likTpr(u,lpinfo)
            p2 = p2Lik+p2Pri
            pr = npy.exp(p2-p1);
            if (pr>=1.0) or (npy.random.random_sample()<=pr):
                spls[k+1] = u.copy();                # Store accepted sample
                meta_info[k+1] = [pr,p2Lik,p2Pri]    # and its meta information
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
                meta_info[k+1,0] = pr               # acceptance probability of failed sample
                meta_info[k+1,1:] = meta_info[k,1:] # Posterior probability of sample k that has been retained
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
                    p2Lik,p2Pri = likTpr(u,lpinfo)
                    p2 = p2Lik+p2Pri
                    tryspls.append(u.copy())
                    trypost.append(p2)
                    alpha = getAlpha(tryspls,trypost);
                    if (alpha >= 1.0) or (npy.random.random_sample() < alpha):
                        accept = True;
                        spls[k+1] = u.copy();                # Store accepted sample
                        meta_info[k+1] = [alpha,p2Lik,p2Pri] # and its meta information
                        p1 = p2;
                        if p1 > pmode:
                            pmode = p1 ;
                            cmode = spls[k+1] ;
                if not accept:
                    spls[k+1]=spls[k] ;
                    meta_info[k+1,0]  = alpha            # acceptance probability of failed sample
                    meta_info[k+1,1:] = meta_info[k,1:]   # Posterior probability of sample k that has been retained
                    rej   = rej   + 1;
                    rejsc = rejsc + 1;
                    if outofbound:
                        rejlim  = rejlim + 1;
            else:
                print('Unknown MCMC method %s -> Quit'%(method))
                quit()
            # Done with if over methods
        # Done with if over original accept
        if ((k+1)%ofreq==0 and tmp_file != 'None'):
            print('No. steps: %d, No. of rej:%d'%(k+1,rej))
            fout = open(tmp_file, 'ab')
            npy.savetxt(fout, spls[k-ofreq+1:k+1,:], fmt='%.8e',delimiter=' ', newline='\n')
            fout.close()
    # Done loop over all steps

    # return output dictionary: samples, MAP sample and its posterior probability, overall acceptance probability
    # and probability of having sample inside prior bounds, overall number of samples rejected, and rejected
    # due to being out of bounds.
    mcmcRes={}
    mcmcRes['chain' ] = spls                     # chain
    mcmcRes['cmap'  ] = cmode                    # MAP state
    mcmcRes['pmap'  ] = pmode                    # MAP log posterior
    mcmcRes['accr'  ] = 1.0-float(rej)/nsteps    # acceptance rate (overall)
    mcmcRes['accb'  ] = 1.0-float(rejlim)/nsteps # samples inside bounds
    mcmcRes['rejAll'] = rej                      # overall no. of samples rejected
    mcmcRes['rejAll'] = rejlim                   # no. of samples rejected due to being outside bounds
    mcmcRes['minfo' ] = meta_info                # acceptance probability and log posterior for each state
    return mcmcRes
