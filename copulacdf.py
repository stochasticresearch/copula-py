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
def copulacdf(family, u, *args):
    """ Generates values of the Gaussian copula
    
    Inputs:
    u -- u is an N-by-P matrix of values in [0,1], representing N
         points in the P-dimensional unit hypercube.  
    
    rho -- a P-by-P correlation matrix, the first argument required for the gaussian copula
    
    Outputs:
    y -- the value of the Gaussian Copula
    """

    num_var_args = len(args)
    family_lc = family.lower()
    if(family_lc=='gaussian'):
        if(num_var_args!=1):
            raise ValueError("Gaussian family requires one additional argument -- rho (correlation matrix) [P x P]")
        rho = args[1]
        # TODO: some error checking to see if rho is a P x P matrix?
        y = gaussian_copula_cdf(u, rho)
        
    elif(family_lc=='t'):
        pass
    elif(family_lc=='clayton'):
        if(num_var_args!=1):
            raise ValueError("Clayton family requires one additional argument -- alpha [scalar]")
        alpha = args[1]
        # TODO: make sure alpha is a scalar
        y = clayton_copula_cdf(u, alpha)
    elif(family_lc=='frank'):
        if(num_var_args!=1):
            raise ValueError("Frank family requires one additional argument -- alpha [scalar]")
        alpha = args[1]
        # TODO: make sure alpha is a scalar
        y = frank_copula_cdf(u, alpha)
    elif(family_lc=='gumbel'):
        if(num_var_args!=1):
            raise ValueError("Gumbel family requires one additional argument -- alpha [scalar]")
        alpha = args[1]
        # TODO: make sure alpha is a scalar
        y = gumbel_copula_cdf(u, alpha)
    else:
        raise ValueError("Unrecognized family of copula")
    
    return y

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
    # C(u1,u2) = (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)
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
    # C(u1,u2) = -(1/alpha)*log(1 + (exp(-alpha*u1)-1)*(exp(-alpha*u2)-1)/(exp(-alpha)-1))
    if(alpha==0):
        y = np.prod(u,1)
    else:
        tmp1 = np.exp(-alpha*np.sum(u,1)) - np.sum(np.exp(-alpha*u),1)
        y = -np.log( (np.exp(-alpha) + tmp1) / np.expm1(-alpha)) / alpha;
        
    return y

def gumbel_copula_cdf(u, alpha):
    # C(u1,u2) = exp(-((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))
    if(alpha<1):
        # TODO: some sort of error here?
        return None
    elif(alpha==1):
        y = np.prod(u,1)
    else:
        pass

def test_python_vs_matlab(family):
    # generate U1, U2
    n = 10
    p = 2
    
    # generate all u1,u2 combinations
    eps = np.finfo(float).eps
    u = np.linspace(0+eps,1-eps,n)
    UU = np.meshgrid(u,u)
    U1 = np.reshape(UU[0], (UU[0].shape[0]*UU[0].shape[1], 1))
    U2 = np.reshape(UU[1], (UU[1].shape[0]*UU[1].shape[1], 1))
    U = np.concatenate((U1,U2),axis=1)
    
    rho = 0.8
    Rho = np.array([[1,rho],[rho,1]])
    
    alpha = 0.3
    
    # test the Gaussian Copula against Matlab
    # TODO: make python execute the matlab script which generates these samples
    matlab_data = scipy.io.loadmat('matlab/copula_cdf_test.mat')
    
    if(family.lower()=='gaussian'):
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
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(gaussian_copula_cdf_python,UU[0].shape)
        
        plot_3d(X,Y,Z, 'Gaussian Copula CDF')
        
    elif(family.lower()=='clayton'):
        clayton_copula_cdf_python = clayton_copula_cdf(U,alpha)
        clayton_copula_cdf_matlab = matlab_data['clayton_copula_cdf']
        clayton_copula_cdf_matlab = clayton_copula_cdf_matlab[:,0]
        
        # compare the two
        clayton_copula_test_result = np.allclose(clayton_copula_cdf_python,clayton_copula_cdf_matlab)
        if(clayton_copula_test_result):
            print 'Clayton Copula Python calculation matches Matlab!'
        else:
            print 'Clayton Copula Python calculation does NOT match Matlab!'
            
        # plot the Clayton Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(clayton_copula_cdf_python,UU[0].shape)
        
        plot_3d(X,Y,Z, 'Clayton Copula CDF')
        
    elif(family.lower()=='frank'):
        frank_copula_cdf_python = frank_copula_cdf(U,alpha)
        frank_copula_cdf_matlab = matlab_data['frank_copula_cdf']
        frank_copula_cdf_matlab = frank_copula_cdf_matlab[:,0]
        
        # compare the two
        frank_copula_test_result = np.allclose(frank_copula_cdf_python,frank_copula_cdf_matlab)
        if(frank_copula_test_result):
            print 'Frank Copula Python calculation matches Matlab!'
        else:
            print 'Frank Copula Python calculation does NOT match Matlab!'
            
        # plot the Clayton Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(frank_copula_cdf_python,UU[0].shape)
        
        plot_3d(X,Y,Z, 'Frank Copula CDF')

def plot_3d(X,Y,Z, titleStr):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.xlabel('U1')
    plt.ylabel('U2')
    plt.title(titleStr)
    plt.show()

if __name__=='__main__':
    import scipy.io
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    test_python_vs_matlab('Gaussian')
    test_python_vs_matlab('Clayton')
    test_python_vs_matlab('Frank')
    