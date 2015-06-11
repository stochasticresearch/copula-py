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

from scipy.stats import mvn                     # contains inverse CDF of Multivariate Gaussian
from scipy.stats import norm                    # contains PDF of Gaussian
from scipy.stats import t
from statsmodels.sandbox.distributions import multivariate as mvt

"""
copulacdf.py contains routines which provide Copula CDF values 
"""

def copulacdf(family, u, *args):
    """ Generates values of a requested copula family
    
    Inputs:
    u -- u is an N-by-P matrix of values in [0,1], representing N
         points in the P-dimensional unit hypercube.  
    
    rho -- a P-by-P correlation matrix, the first argument required for the Gaussian copula
    alpha -- a scalar argument describing the dependency for Frank, Gumbel, and Clayton copula's
    
    Outputs:
    y -- the value of the Gaussian Copula
    """
    n  = u.shape[0]
    p  = u.shape[1]

    num_var_args = len(args)
    family_lc = family.lower()
    if(family_lc=='gaussian'):
        if(num_var_args!=1):
            raise ValueError("Gaussian family requires one additional argument -- rho (correlation matrix) [P x P]")
        rho = args[0]
        rho_expected_shape = (p,p)
        if(type(rho)!=np.ndarray or rho.shape!=rho_expected_shape):
            raise ValueError("Gaussian family requires rho to be of type numpy.ndarray with shape=[P x P]")
        y = _gaussian(u, rho)
        
    elif(family_lc=='t'):
        if(num_var_args!=2):
            raise ValueError("T family requires two additional arguments -- rho (correlation matrix) [P x P] and nu [scalar]")
        rho = args[0]
        nu  = args[1]
        rho_expected_shape = (p,p)
        if(type(rho)!=np.ndarray or rho.shape!=rho_expected_shape):
            raise ValueError("T family requires rho to be of type numpy.ndarray with shape=[P x P]")
        y = _t(u, rho, nu)        
    elif(family_lc=='clayton'):
        if(num_var_args!=1):
            raise ValueError("Clayton family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        y = _clayton(u, alpha)
    elif(family_lc=='frank'):
        if(num_var_args!=1):
            raise ValueError("Frank family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        y = _frank(u, alpha)
    elif(family_lc=='gumbel'):
        if(num_var_args!=1):
            raise ValueError("Gumbel family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        y = _gumbel(u, alpha)
    else:
        raise ValueError("Unrecognized family of copula")
    
    return y

def _gaussian(u, rho):
    """ Generates values of the Gaussian copula
    
    Inputs:
    u -- u is an N-by-P matrix of values in [0,1], representing N
         points in the P-dimensional unit hypercube.  
    rho -- a P-by-P correlation matrix.
    
    Outputs:
    y -- the value of the Gaussian Copula
    """
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

def _t(u, rho, nu):
    """ Generates values of the T copula
    
    Inputs:
    u -- u is an N-by-P matrix of values in [0,1], representing N
         points in the P-dimensional unit hypercube.  
    rho -- a P-by-P correlation matrix.
    nu  -- degrees of freedom for T Copula
    
    Outputs:
    y -- the value of the T Copula
    """
    n  = u.shape[0]
    p  = u.shape[1]
    loIntegrationVal = -40
    lo = np.full((1,p), loIntegrationVal)        # more accuracy, but slower :/
    hi = t.ppf(u, nu)
    
    mu = np.zeros(p)
    
    y = np.zeros(n)
    for ii in np.arange(n):
        x = hi[ii,:]
        x[x<-40] = -40
        p = mvt.mvstdtprob(lo[0], x, rho, nu)
        y[ii] = p
    
    return y
    
def _clayton(u, alpha):
    # C(u1,u2) = (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)
    if(alpha<0):
        raise ValueError("Clayton family -- invalid alpha argument. alpha must be >=0")
    elif(alpha==0):
        y = np.prod(u,1)
    else:
        tmp1 = np.power(u, -alpha) 
        tmp2 = np.sum(tmp1,1) - 1
        y = np.power(tmp2, -1.0/alpha)
        
    return y

def _frank(u, alpha):
    # C(u1,u2) = -(1/alpha)*log(1 + (exp(-alpha*u1)-1)*(exp(-alpha*u2)-1)/(exp(-alpha)-1))
    if(alpha==0):
        y = np.prod(u,1)
    else:
        tmp1 = np.exp(-alpha*np.sum(u,1)) - np.sum(np.exp(-alpha*u),1)
        y = -np.log( (np.exp(-alpha) + tmp1) / np.expm1(-alpha)) / alpha;
        
    return y

def _gumbel(u, alpha):
    # C(u1,u2) = exp(-( (-log(u1))^alpha + (-log(u2))^alpha )^(1/alpha))
    n = u.shape[0]
    p = u.shape[1]
    
    if(alpha<1):
        raise ValueError("Gumbel family -- invalid alpha argument. alpha must be >=1")
    elif(alpha==1):
        y = np.prod(u,1)
    else:
        # TODO: NaN checking like Matlab here would be nice :)
        exparg = np.zeros(n)
        for ii in np.arange(p):
            tmp1 = np.power(-1*np.log(u[:,ii]), alpha)
            exparg = exparg + tmp1
        exparg = np.power(exparg, 1.0/alpha)
        y = np.exp(-1*exparg)
        
    return y

def test_python_vs_matlab(family):
    # generate U1, U2
    n = 10
    p = 2
    
    # generate all u1,u2 combinations
    eps = np.finfo(float).eps
    u = np.linspace(0+eps,1-eps,n)
    UU = np.meshgrid(u,u)
    U2 = np.reshape(UU[0], (UU[0].shape[0]*UU[0].shape[1], 1))
    U1 = np.reshape(UU[1], (UU[1].shape[0]*UU[1].shape[1], 1))
    U = np.concatenate((U1,U2),axis=1)
        
    rho = 0.8
    nu = 2
    Rho = np.array([[1,rho],[rho,1]])
    
    alpha = 0.3
    
    # test the python data against Matlab
    # TODO: make python execute the matlab script which generates these samples
    matlab_data = scipy.io.loadmat('matlab/copulacdf_test.mat')
    
    if(family.lower()=='gaussian'):
        gaussian_copula_cdf_python = copulacdf(family,U,Rho)
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
        
        plot_utils.plot_3d(X,Y,Z, 'Gaussian Copula CDF')
    
    elif(family.lower()=='t'):
        t_copula_cdf_python = copulacdf(family,U,Rho,nu)
        t_copula_cdf_matlab = matlab_data['t_copula_cdf']
        t_copula_cdf_matlab = t_copula_cdf_matlab[:,0]
        
        # compare the two
        t_copula_test_result = np.allclose(t_copula_cdf_python,t_copula_cdf_matlab,atol=0.01)   # a high tolerance required
                                                                                                # b/c of way that mvt
                                                                                                # is implemented in python
        if(t_copula_test_result):
            print 'T Copula Python calculation matches Matlab!'
        else:
            print 'T Copula Python calculation does NOT match Matlab!'
            
        # plot the Guassian Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(t_copula_cdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'T Copula CDF')
    
    elif(family.lower()=='clayton'):
        clayton_copula_cdf_python = copulacdf(family,U,alpha)
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
        
        plot_utils.plot_3d(X,Y,Z, 'Clayton Copula CDF')
        
    elif(family.lower()=='frank'):
        frank_copula_cdf_python = copulacdf(family,U,alpha)
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
        
        plot_utils.plot_3d(X,Y,Z, 'Frank Copula CDF')
        
    elif(family.lower()=='gumbel'):
        alpha = 1.5
        gumbel_copula_cdf_python = copulacdf(family,U,alpha)
        gumbel_copula_cdf_matlab = matlab_data['gumbel_copula_cdf']
        gumbel_copula_cdf_matlab = gumbel_copula_cdf_matlab[:,0]
        
        # compare the two
        gumbel_copula_test_result = np.allclose(gumbel_copula_cdf_python,gumbel_copula_cdf_matlab)
        if(gumbel_copula_test_result):
            print 'Gumbel Copula Python calculation matches Matlab!'
        else:
            print 'Gumbel Copula Python calculation does NOT match Matlab!'
            
        # plot the Clayton Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(gumbel_copula_cdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'Gumbel Copula CDF')

if __name__=='__main__':
    import scipy.io
    import plot_utils
    
    #test_python_vs_matlab('Gaussian')
    test_python_vs_matlab('T')
    #test_python_vs_matlab('Clayton')
    #test_python_vs_matlab('Frank')
    #test_python_vs_matlab('Gumbel')
    