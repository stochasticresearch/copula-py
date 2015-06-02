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

from numpy.linalg import solve
from numpy.linalg import cholesky
from numpy.linalg import LinAlgError

"""
copulapdf.py contains routines which provide Copula PDF values 
"""

def copulapdf(family, u, *args):
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
        # TODO: fix me :)
        return None
    elif(family_lc=='clayton'):
        if(num_var_args!=1):
            raise ValueError("Clayton family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        if(type(alpha)!=float):
            raise ValueError('Clayton family requires a scalar alpha value')
        y = _clayton(u, alpha)
    elif(family_lc=='frank'):
        if(num_var_args!=1):
            raise ValueError("Frank family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        if(type(alpha)!=float):
            raise ValueError('Clayton family requires a scalar alpha value')
        y = _frank(u, alpha)
    elif(family_lc=='gumbel'):
        if(num_var_args!=1):
            raise ValueError("Gumbel family requires one additional argument -- alpha [scalar]")
        alpha = args[0]
        if(type(alpha)!=float):
            raise ValueError('Clayton family requires a scalar alpha value')
        y = _gumbel(u, alpha)
    else:
        raise ValueError("Unrecognized family of copula")
    
    return y

def _gaussian(u, rho):
    matlab_data = scipy.io.loadmat('matlab/copulapdf_test.mat')
    
    try:
        R = cholesky(rho)
    except LinAlgError:
        raise ValueError('Provided Rho matrix is not Positive Definite!')
    
    x = norm.ppf(u)
    logSqrtDetRho = np.sum(np.log(np.diag(R)))
    z = solve(R,x.T)
    z = z.T    
    y = np.exp(-0.5 * np.sum(  np.power(z,2) - np.power(x,2) ,  axis=1  ) - logSqrtDetRho)
    
    return y

def _t(u, rho, nu):
    return None

def _clayton(u, alpha):
    return None

def _frank(u, alpha):
    return None

def _gumbel(u, alpha):
    return None

def test_python_vs_matlab(family):
    # generate U1, U2
    n = 10
    p = 2
    
    # generate all u1,u2 combinations
    eps = np.finfo(float).eps
    u = np.linspace(0.1,0.9,n)
    UU = np.meshgrid(u,u)
    U2 = np.reshape(UU[0], (UU[0].shape[0]*UU[0].shape[1], 1))
    U1 = np.reshape(UU[1], (UU[1].shape[0]*UU[1].shape[1], 1))
    U = np.concatenate((U1,U2),axis=1)
    
    rho = 0.8
    Rho = np.array([[1,rho],[rho,1]])
    
    alpha = 0.3
    
    # test the python data against Matlab
    # TODO: make python execute the matlab script which generates these samples
    matlab_data = scipy.io.loadmat('matlab/copulapdf_test.mat')
    
    if(family.lower()=='gaussian'):
        gaussian_copula_pdf_python = copulapdf(family,U,Rho)
        gaussian_copula_pdf_matlab = matlab_data['gaussian_copula_pdf'][:,0]
        
        # compare the two
        gaussian_copula_test_result = np.allclose(gaussian_copula_pdf_python,gaussian_copula_pdf_matlab)
        if(gaussian_copula_test_result):
            print 'Gaussian Copula Python calculation matches Matlab!'
        else:
            print 'Gaussian Copula Python calculation does NOT match Matlab!'
            
        # plot the Guassian Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(gaussian_copula_pdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'Gaussian Copula PDF')
        
if __name__=='__main__':
    import scipy.io
    import plot_utils
    
    test_python_vs_matlab('Gaussian')