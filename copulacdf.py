#!/usr/bin/env python

import math
import numpy as np

from scipy.stats import mvn                     # contains inverse CDF of Multivariate Gaussian
from scipy.stats import norm                    # contains PDF of Gaussian

# TODO: make this more like copularnd in matlab/octave, one function w/ options
# which selects

"""
copulacdf.py contains routines which provide Copula CDF values 
"""

# TODO: consolidate all other functions into this, instead of having
# a million functions all do similar things, switch w/ family argument
# similar to octave and matlab
def copulacdf(family, u, **kwargs):
    pass

def gaussian_copula_cdf(u, rho):
    """ Generates values of the Gaussian copula
    
    Inputs:
    u -- u is an N-by-P matrix of values in [0,1], representing N
         points in the P-dimensional unit hypercube.  
    rho -- a P-by-P correlation matrix.
    
    Outputs:
    y -- the value of the Gaussian Copula
    """

    # TODO: consider adding some error checking for input arguments
    
    n  = u.shape[0]
    p  = u.shape[1]
    lo = np.full((1,p), -10)
    hi = norm.ppf(u)
    
    mu = np.zeros(p)
    
    # need to use ppf(q, loc=0, scale=1) as replacement for norminv
    # need to use mvn.mvnun as replacement for mvncdf
    # the upper bound needs to be the output of the ppf call, right now it is set to random above
    y = np.zeros(n)
    # I don't know if mvnun is vectorized, I couldn't get that to work
    for ii in np.arange(n):
        # do some error checking.  if we have any -inf or inf values,
        # 
        p,i = mvn.mvnun(lo, hi[ii,:], mu, rho)
        y[ii] = p
    
    return y

def t_copula_cdf(u, rho, n):
    pass

def clayton_copula_cdf(u, alpha):
    if(alpha<0):
        # TODO: some sort of error here?
        return None
    elif(alpha==0):
        y = np.prod(u,1)
    else:
        tmp1 = np.power(u, -alpha) 
        tmp2 = np.sum(tmp1,1) - 1
        y = np.power(tmp2, -1.0/alpha)
        
    return y

def frank_copula_cdf(u, alpha):
    if(alpha==0):
        y = np.prod(u,1)
    else:
        tmp1 = np.power(u, -alpha) 
        tmp2 = np.sum(tmp1,1) - 1
        y = np.power(tmp2, -1.0/alpha)
        
        #p = -log((exp(-alpha) + (exp(-alpha.*sum(u,2)) - sum(exp(-alpha.*u),2))) ./ expm1(-alpha)) ./ alpha;
        
    return y

def gumbel_copula_cdf(u, alpha):
    if(alpha<0):
        # TODO: some sort of error here?
        return None
    elif(alpha==0):
        y = np.prod(u,1)
    else:
        pass

if __name__=='__main__':
    import scipy.io
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    n = 10
    p = 2
    
    eps = np.finfo(float).eps
    u = np.linspace(0+eps,1-eps,n)
    UU = np.meshgrid(u,u)
    U1 = np.reshape(UU[0], (UU[0].shape[0]*UU[0].shape[1], 1))
    U2 = np.reshape(UU[1], (UU[1].shape[0]*UU[1].shape[1], 1))
    U = np.concatenate((U1,U2),axis=1)
    
    rho = 0.8
    Rho = np.array([[1,rho],[rho,1]])
    
    # test the Gaussian Copula against Matlab
    # TODO: make python execute the matlab script which generates these samples
    matlab_data = scipy.io.loadmat('matlab/copula_cdf_test.mat')
    
    gaussian_copula_cdf_python = gaussian_copula_cdf(U,Rho)
    gaussian_copula_cdf_matlab = matlab_data['gaussian_copula_cdf']
    gaussian_copula_cdf_matlab = gaussian_copula_cdf_matlab[:,0]
    
    # compare the two
    gaussian_copula_test_result = np.allclose(gaussian_copula_cdf_python,gaussian_copula_cdf_matlab)
    if(gaussian_copula_test_result):
        print 'Gaussian Copula Python calculation matches Matlab!'
    else:
        print 'Gaussian Copula Python calculation does NOT match Matlab!'
    # plot the Guassian Copula for fun
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = UU[0]
    Y = UU[1]
    Z = np.reshape(gaussian_copula_cdf_python,UU[0].shape)
    
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.xlabel('U1')
    plt.ylabel('U2')
    plt.title('Gaussian Copula CDF')
    plt.show()