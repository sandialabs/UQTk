import numpy as npy
import scipy.stats 

def bayesevid_ChibJeliazkov(spls_star,spls,llp_star,llp,likTpr,lpinfo,prop_cov,coveps,gamma):
    '''
    Implements the one block sampling procedure from Chib and Jeliazkov 2001 JASA paper
    Inputs: 
        spls_star - a sample close to the MAP estimate
        spls      - MCMC samples from the posterior
        llp_star  - log-likelihood + log-prior at spls_star
        llp.      - log-likelihood + log-prior at spls
        likTpr    - function retrieving log-likelihood + log-prior for a custom sample
        lpinfo    - auxiliary info needed by likTpr
        prop_cov  - proposal covariance (from MCMC)
        coveps    - factor to be added to diagonal to ensure matrix is stricly positive definite
        gamma     - proposal covariance scaling factor
    '''
    cdim   = spls_star.shape[0]
    npSpls = spls.shape[0]
    # prepare covariance
    sigcv   = 2.4*gamma/npy.sqrt(cdim)
    covPeps = sigcv**2*(prop_cov + coveps * npy.identity(cdim))
    Rchol   = scipy.linalg.cholesky(covPeps)
    # compute numerator
    num = 0.0
    for i in range(npSpls):
        if llp_star-llp[i]<0.0:
            alpha = npy.exp(llp_star-llp[i])
        else:
            alpha = 1.0
        #print(spls[i], llp_star, covPeps)
        val = scipy.stats.multivariate_normal.pdf(spls[i], mean=spls_star, cov=covPeps)
        num = num+val*alpha
    # compute denominator
    den = 0.0
    for i in range(npSpls):
        u = spls_star + npy.dot(npy.random.randn(1,cdim),Rchol)[0];
        p2Lik,p2Pri = likTpr(u,lpinfo)
        if p2Lik + p2Pri - llp_star < 0.0:
            alpha = npy.exp(p2Lik + p2Pri - llp_star)
        else:
            alpha = 1.0
        den = den + alpha
    return llp_star-npy.log(num)+npy.log(den)

def log_lxp(x,info):
    '''
    log-likelihood & log-prior for test cases: normalizing factors for MVNs
    '''
    mu  = info['mu']
    sgm = info['cov']
    return npy.log(scipy.stats.multivariate_normal.pdf(x, mean=mu, cov=sgm)*npy.sqrt(2*npy.pi*npy.linalg.det(sgm))),0.0

def checkBE(nspls,mu,sgm):
    '''
    setting Bayes evidence calculations for canonical test cases
    '''
    # settings
    pcov = sgm.copy()
    lpinfo = {'mu':mu,'cov':sgm}    
    # generate samples from a 2D MVN
    spls = npy.random.multivariate_normal(mu, sgm, nspls)
    llp  = npy.zeros(nspls)
    for i in range(nspls):
        ll,lp = log_lxp(spls[i],lpinfo)
        llp[i] = ll+lp
    ll,lp = log_lxp(mu,lpinfo)
    llp_star = ll+lp
    return bayesevid_ChibJeliazkov(mu,spls,llp_star,llp,log_lxp,lpinfo,pcov,1.e-12,0.8)

def checkErr(dims):
    '''
    convergence of Bayes evidence calculations for canonical test cases
    '''
    # default settings
    m=npy.array([1.1*i  for i in range(dims)])
    b=npy.array([1.1**i for i in range(dims)])
    cvm = npy.identity(dims)
    for i in range(dims):
        cvm[i] = cvm[i]*b[i]
    truth=npy.log(npy.sqrt(2*npy.pi*npy.prod(b)))
    errAll=[]
    for nspl in [1000,10000,100000]:
        err = []
        for k in range(10):
            err.append(checkBE(nspl,m,cvm))
        err = (npy.array(err)-truth)/truth
        print(nspl,err)
        errAll.append(err)
    return errAll

