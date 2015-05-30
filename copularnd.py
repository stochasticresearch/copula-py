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

"""
copularnd.py contains routines which provide samples of a copula density
"""

def copularnd(family, n, *args):
    """ Generates values of the Gaussian copula
    
    Inputs:
    family -- Should be either 'gaussian', 't', 'clayton', 'frank', or 'gumbel'
    n      -- the number of samples to generate
    args   -- variable number of arguments depending on which type of copula you are trying to simulate
                Gaussian -- should be provided a 2x2 rho matrix as a numpy array datatype
                t        -- should be provided a rho matrix and a nu value
                Clayton/Frank/Gumbel - should be provided a scalar alpha value
    
    Outputs:
    U -- a 2xn matrix of samples from the copula density chosen
    """

    num_var_args = len(args)
    family_lc = family.lower()
    if(family_lc=='gaussian'):
        if(num_var_args!=1):
            raise ValueError("Gaussian family requires one additional argument -- rho (correlation matrix) [P x P]")
        rho = args[0]
        rho_expected_shape = (2,2)
        if(type(rho)!=np.ndarray or rho.shape!=rho_expected_shape):
            raise ValueError("Gaussian family requires rho to be of type numpy.ndarray with shape=[P x P]")
        U = _gaussian(n, rho)
        
    elif(family_lc=='t'):
        return None
    elif(family_lc=='clayton'):
        if(num_var_args!=1):
            raise ValueError("Clayton family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        U = _clayton(n, alpha)
    elif(family_lc=='frank'):
        if(num_var_args!=1):
            raise ValueError("Frank family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        U = _frank(n, alpha)
    elif(family_lc=='gumbel'):
        if(num_var_args!=1):
            raise ValueError("Gumbel family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        U = _gumbel(n, alpha)
    else:
        raise ValueError("Unrecognized family of copula")
    
    return U

def _gaussian(n, Rho):
    """
    Generates samples from the Gaussian Copula, w/ dependency
    matrix described by Rho.  Rho should be a numpy 2x2 matrix.
    It is assumed that we have a 0 mean.
    """
    mu = [0,0]
    # TODO: some error checking on Rho to make sure it is 2x2?
    y = multivariate_normal(mu,Rho)
    
    # generate samples of the multivariate normal distribution
    # and apply normal cdf to generate U
    U = np.empty((2,n))
    
    for ii in range(0,n):
        mvnData = y.rvs()
        U[0][ii] = norm.cdf(mvnData[0])
        U[1][ii] = norm.cdf(mvnData[1])
    
    return U
    
# TODO :)
def _t():
    pass

# We generate the Archimedean Copula's as follows:
# Random pairs from these copulae can be generated sequentially: first
# generate u1 as a uniform r.v.  Then generate u2 from the conditional
# distribution F(u2 | u1; alpha) by generating uniform random values, then
# inverting the conditional CDF.
# This method is outlined in Nelsen's Introduction to Copula's

def _clayton(n, alpha):
    u1 = uniform.rvs(size=n)
    p = uniform.rvs(size=n)
    if(alpha<np.spacing(1)):
        u2 = p
    else:
        u2 = u1*np.power((np.power(p,(-alpha/(1.0+alpha))) - 1 + np.power(u1,alpha)),(-1.0/alpha))
        
    U = np.vstack((u1,u2))
    return U

def _frank(n, alpha):
    u1 = uniform.rvs(size=n)
    p = uniform.rvs(size=n)
    if abs(alpha) > math.log(sys.float_info.max):
        u2 = (u1 < 0).astype(int) + np.sign(alpha)*u1  # u1 or 1-u1
    elif abs(alpha) > math.sqrt(np.spacing(1)):
        u2 = -1*np.log((np.exp(-alpha*u1)*(1-p)/p + np.exp(-alpha))/(1 + np.exp(-alpha*u1)*(1-p)/p))/alpha
    else:
        u2 = p
    
    U = np.vstack((u1,u2))
    return U

def _gumbel(n, alpha):
    if alpha < 1:
        raise ValueError('Bad Gumbel Dependency Parameter!')
    if alpha < (1 + math.sqrt(np.spacing(1))):
        u1 = uniform.rvs(size=n);
        u2 = uniform.rvs(size=n);
    else:
        # use the Marshal-Olkin method
        # Generate gamma as Stable(1/alpha,1), c.f. Devroye, Thm. IV.6.7
        u = (uniform.rvs(size=n) - .5) * math.pi # Generate n uniformly distributed RV's between -pi/2 and pi/2
        u2 = u + math.pi/2
        e  = -1*np.log(uniform.rvs(size=n))
        t = np.cos(u - u2/alpha)/e
        gamma = np.power(np.sin(u2/alpha)/t,(1.0/alpha)) * t/np.cos(u);
        
        # Frees&Valdez, eqn 3.5
        u1 = np.exp(-1* (np.power(-1*np.log(uniform.rvs(size=n)), 1.0/alpha) / gamma) )
        u2 = np.exp(-1* (np.power(-1*np.log(uniform.rvs(size=n)), 1.0/alpha) / gamma) )
        
    U = np.vstack((u1,u2))
    return U

if __name__=='__main__':
    import matplotlib.pyplot as plt
    n = 1000
    rh = 0.6
    Rho = np.array([[1,rh],[rh,1]])
    alpha = 5
    
    Ug = copularnd('gaussian', n, Rho)
    Uc  = copularnd('clayton', n,alpha)
    Uf  = copularnd('frank', n,alpha)
    Ugu = copularnd('gumbel', n,alpha)
    
    plt.scatter(Ug[0,:],Ug[1,:])
    plt.title('Gaussian Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Uc[0,:],Uc[1,:])
    plt.title('Clayton Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Uf[0,:],Uf[1,:])
    plt.title('Frank Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    plt.scatter(Ugu[0,:],Ugu[1,:])
    plt.title('Gumbel Copula Samples')
    plt.grid()
    plt.axis((0,1,0,1))
    plt.show()
    
    