#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#******************************************************************************
#* 
#* Copyright (C) 2015  Kiran Karra <kiran.karra@gmail.com>
#*
#* This program is free software: you can redistribute it and/or modify
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with this program.  If not, see <http://www.gnu.org/licenses/>.
#******************************************************************************

import math
import numpy as np

import sys

from scipy.stats import norm                    # contains PDF of Gaussian
from scipy.stats import multivariate_normal
from scipy.stats import uniform
from scipy.stats import gamma
from scipy.stats import logser
from scipy.stats import t
from stable import stable
from statsmodels.sandbox.distributions import multivariate as mvt

"""
copularnd.py contains routines which provide samples of a copula density
"""

def copularnd(family, M, *args):
    """ Generates values of the Gaussian copula
    
    Inputs:
    family -- Should be either 'gaussian', 't', 'clayton', 'frank', or 'gumbel'
    M      -- the number of samples to generate
    args   -- variable number of arguments depending on which type of copula you are trying to simulate
                Gaussian -- should be provided a NxN rho matrix as a numpy array datatype
                t        -- should be provided a NxN rho matrix and a nu value
                Clayton/Frank/Gumbel - should be provided a N for the dimensionality, and a scalar alpha value
    
    Outputs:
    U -- a M x N matrix of samples from the copula density chosen
    """

    num_var_args = len(args)
    family_lc = family.lower()
    if(family_lc=='gaussian'):
        if(num_var_args!=1):
            raise ValueError("Gaussian family requires one additional argument -- rho (correlation matrix) [P x P]")
        rho = args[0]
        shape0 = rho.shape[0]
        shape1 = rho.shape[1]
        if(shape0!=shape1):
            raise ValueError("Gaussian family requires rho to be of type numpy.ndarray with shape=[P x P]")
        U = _gaussian(M, rho)
    elif(family_lc=='t'):
        if(num_var_args!=2):
            raise ValueError("T family requires two additional argument -- rho (correlation matrix) [P x P] and nu [scalar]")
        rho = args[0]
        shape0 = rho.shape[0]
        shape1 = rho.shape[1]
        if(shape0!=shape1):
            raise ValueError("T family requires rho to be of type numpy.ndarray with shape=[P x P]")
        nu = args[1]
        U = _t(M, rho, nu)
    elif(family_lc=='clayton'):
        if(num_var_args!=2):
            raise ValueError("Clayton family requires two additional arguments -- N, alpha [scalar]")
        N     = args[0]
        alpha = args[1]
        U = _clayton(M, N, alpha)
    elif(family_lc=='frank'):
        if(num_var_args!=2):
            raise ValueError("Frank family requires two additional arguments -- N, alpha [scalar]")
        N     = args[0]
        alpha = args[1]
        U = _frank(M, N, alpha)
    elif(family_lc=='gumbel'):
        if(num_var_args!=2):
            raise ValueError("Gumbel family requires two additional arguments -- N, alpha [scalar]")
        N     = args[0]
        alpha = args[1]
        U = _gumbel(M, N, alpha)
    else:
        raise ValueError("Unrecognized family of copula")
    
    return U

def _gaussian(M, Rho):
    """
    Generates samples from the Gaussian Copula, w/ dependency
    matrix described by Rho.  Rho should be a numpy square matrix.
    It is assumed that we have a 0 mean.
    """
    N = Rho.shape[0]
    mu = np.zeros(N)
    y = multivariate_normal(mu,Rho)
    mvnData = y.rvs(size=M)
    U = norm.cdf(mvnData)
    
    return U
    
def _t(M, Rho, nu):
    N = Rho.shape[0]
    mu = np.zeros(N)        # zero mean
    x = mvt.multivariate_t_rvs(mu,Rho,nu,M) # generate T RV's
    U = t.cdf(x, nu)
    
    return U

# We generate the Archimedean Copula's as follows:
# Random pairs from these copulae can be generated sequentially: first
# generate u1 as a uniform r.v.  Then generate u2 from the conditional
# distribution F(u2 | u1; alpha) by generating uniform random values, then
# inverting the conditional CDF.
# This method is outlined in Nelsen's Introduction to Copula's

def _clayton(M, N, alpha):
    if(N<2):
        raise ValueError('Dimensionality Argument [N] must be an integer >= 2')
    elif(N==2):
        u1 = uniform.rvs(size=M)
        p = uniform.rvs(size=M)
        if(alpha<np.spacing(1)):
            u2 = p
        else:
            u2 = u1*np.power((np.power(p,(-alpha/(1.0+alpha))) - 1 + np.power(u1,alpha)),(-1.0/alpha))
        
        U = np.column_stack((u1,u2))
        
        return U
    else:
        # Algorithm 1 described in both the SAS Copula Procedure, as well as the
        # paper: "High Dimensional Archimedean Copula Generation Algorithm"
        U = np.empty(M,N)
        for ii in range(0,M):
            shape = 1.0/alpha
            loc = 0
            scale = 1
            v = gamma(shape, loc, scale)
            
            # sample N independent uniform random variables
            x_i = uniform.rvs(size=N)

            # TODO: we can probably vectorize this operation
            for jj in range(0,N):
                t = -1*np.log(x_ii[jj]/v)
                U[ii][jj] = np.power((1+t), -1.0/alpha)
            
        return U

def _frank(M, N, alpha):
    if(N<2):
        raise ValueError('Dimensionality Argument [N] must be an integer >= 2')
    elif(N==2):        
        u1 = uniform.rvs(size=M)
        p = uniform.rvs(size=M)
        if abs(alpha) > math.log(sys.float_info.max):
            u2 = (u1 < 0).astype(int) + np.sign(alpha)*u1  # u1 or 1-u1
        elif abs(alpha) > math.sqrt(np.spacing(1)):
            u2 = -1*np.log((np.exp(-alpha*u1)*(1-p)/p + np.exp(-alpha))/(1 + np.exp(-alpha*u1)*(1-p)/p))/alpha
        else:
            u2 = p
        
        U = np.column_stack((u1,u2))
        return U
    else:
        # Algorithm 1 described in both the SAS Copula Procedure, as well as the
        # paper: "High Dimensional Archimedean Copula Generation Algorithm"
        U = np.empty(M,N)
        for ii in range(0,M):
            p = 1 - np.exp(-1*alpha)
            v = logser.rvs(p, size=1)
            
            # sample N independent uniform random variables
            x_i = uniform.rvs(size=N)

            # TODO: we can probably vectorize this operation
            for jj in range(0,N):
                t = -1*np.log(x_ii[jj]/v)
                U[ii][jj] = -1*np.log( 1-np.exp(t)*(1-np.exp(-1*alpha)) )/alpha
            
        return U
    

def _gumbel(M, N, alpha):
    if alpha < 1:
        raise ValueError('Bad Gumbel Dependency Parameter!')
    if(N<2):
        raise ValueError('Dimensionality Argument [N] must be an integer >= 2')
    elif(N==2):
        if alpha < (1 + math.sqrt(np.spacing(1))):
            u1 = uniform.rvs(size=M);
            u2 = uniform.rvs(size=M);
        else:
            # use the Marshal-Olkin method
            # Generate gamma as Stable(1/alpha,1), c.f. Devroye, Thm. IV.6.7
            u = (uniform.rvs(size=M) - .5) * math.pi # Generate M uniformly distributed RV's between -pi/2 and pi/2
            u2 = u + math.pi/2
            e  = -1*np.log(uniform.rvs(size=M))
            t = np.cos(u - u2/alpha)/e
            gamma = np.power(np.sin(u2/alpha)/t,(1.0/alpha)) * t/np.cos(u);
            
            # Frees&Valdez, eqn 3.5
            u1 = np.exp(-1* (np.power(-1*np.log(uniform.rvs(size=M)), 1.0/alpha) / gamma) )
            u2 = np.exp(-1* (np.power(-1*np.log(uniform.rvs(size=M)), 1.0/alpha) / gamma) )
            
        U = np.column_stack((u1,u2))
        return U
    else:
        # Algorithm 1 described in both the SAS Copula Procedure, as well as the
        # paper: "High Dimensional Archimedean Copula Generation Algorithm"
        U = np.empty(M,N)
        for ii in range(0,M):
            a  = 1.0/alpha
            b  = 1
            g  = np.pow(np.cos(math.pi/(2.0*alpha)), alpha)
            d  = 0
            pm = 1
            s = stable(a, b, g, d, pm)
            v = s.rvs()
            
            # sample N independent uniform random variables
            x_i = uniform.rvs(size=N)

            # TODO: we can probably vectorize this operation
            for jj in range(0,N):
                t = -1*np.log(x_ii[jj]/v)
                U[ii][jj] = np.exp(np.power(-1*t, 1.0/alpha))
            
        return U


if __name__=='__main__':
    import matplotlib.pyplot as plt
    M = 1000
    rh = 0.6
    Rho = np.array([[1,rh],[rh,1]])
    nu = 2
    N = 2
    alpha = 5
    
    Ug = copularnd('gaussian', M, Rho)
    Ut = copularnd('t', M, Rho, nu)
    Uc  = copularnd('clayton', M, N, alpha)
    Uf  = copularnd('frank', M, N, alpha)
    Ugu = copularnd('gumbel', M, N, alpha)
    
    plt.scatter(Ug[:,0],Ug[:,1])
    plt.title('Gaussian Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Ut[:,0],Ut[:,1])
    plt.title('T Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Uc[:,0],Uc[:,1])
    plt.title('Clayton Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Uf[:,0],Uf[:,1])
    plt.title('Frank Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Ugu[:,0],Ugu[:,1])
    plt.title('Gumbel Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    