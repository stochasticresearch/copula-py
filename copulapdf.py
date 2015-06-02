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
from scipy.special import gammaln

from numpy.linalg import solve
from numpy.linalg import cholesky
from numpy.linalg import LinAlgError

from copulacdf import copulacdf

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
        rho = args[0]
        rho_expected_shape = (p,p)
        if(type(rho)!=np.ndarray or rho.shape!=rho_expected_shape):
            raise ValueError("T family requires rho to be of type numpy.ndarray with shape=[P x P]")
        nu = int(args[1])       # force to be an integer
        if(nu<1):
            raise ValueError("T family Degrees of Freedom argument must be an integer >= 1")
        return _t(u, rho, nu)
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
    try:
        R = cholesky(rho)
    except LinAlgError:
        raise ValueError('Provided Rho matrix is not Positive Definite!')
    
    x = norm.ppf(u)
    z = solve(R,x.T)
    z = z.T
    logSqrtDetRho = np.sum(np.log(np.diag(R)))
    y = np.exp(-0.5 * np.sum(  np.power(z,2) - np.power(x,2) ,  axis=1  ) - logSqrtDetRho)
    
    return y

def _t(u, rho, nu):
    d = u.shape[1]
    nu = float(nu)
    
    try:
        R = cholesky(rho)
    except LinAlgError:
        raise ValueError('Provided Rho matrix is not Positive Definite!')
    
    ticdf = t.ppf(u, nu)
    
    z = solve(R,ticdf.T)
    z = z.T
    logSqrtDetRho = np.sum(np.log(np.diag(R)))
    const = gammaln((nu+d)/2.0) + (d-1)*gammaln(nu/2.0) - d*gammaln((nu+1)/2.0) - logSqrtDetRho
    sq = np.power(z,2)
    summer = np.sum(np.power(z,2),axis=1)
    numer = -((nu+d)/2.0) * np.log(1.0 + np.sum(np.power(z,2),axis=1)/nu)
    denom = np.sum(-((nu+1)/2) * np.log(1 + (np.power(ticdf,2))/nu), axis=1)
    y = np.exp(const + numer - denom)
    
    return y

def _clayton(u, alpha):
    n = u.shape[0]
    d = u.shape[1]
    if(d>2):
        raise ValueError('Maximum dimensionality supported is 2 for the Clayton Copula Family')
    if alpha<0:
        raise ValueError('Dependency parameter for Clayton copula must be >= 0')
    elif alpha==0:
        y = np.ones((n,1))
    else:
        # below is the closed form of d2C/dudv of the Clayton copula
        y = (alpha + 1) * np.power( u[:,0]*u[:,1], -1*(alpha+1) ) * np.power( np.power(u[:,0], -alpha) + np.power(u[:,1], -alpha) - 1, -1*(2*alpha+1)/alpha )
    
    return y

def _frank(u, alpha):
    if alpha == 0:
        y = ones(n,1);
    else:
        summer = np.sum(u,1)
        differ = np.diff(u,1,1)
        differ = differ[:,0]
        denom = np.power(np.cosh(alpha*differ/2)*2 - np.exp(alpha*(summer-2)/2) - np.exp(-alpha*summer/2), 2)
        y = alpha*(1-np.exp(-alpha)) / denom
        
    return y

def _gumbel(U, alpha):
    n = U.shape[0]
    d = U.shape[1]
    if(d>2):
        raise ValueError('Maximum dimensionality supported is 2 for the Gumbel Copula Family')
    
    if(alpha < 1):
        raise ValueError('Bad dependency parameter for Gumbel copula')
    elif alpha==1:
        y = np.ones((n,1))
    else:
        # below is the closed form of d2C/dudv of the Gumbel copula
        C = copulacdf('Gumbel', U, alpha)
        u = U[:,0]
        v = U[:,1]
        p1 = C*1.0/(u*v)*np.power(np.power(-1*np.log(u),alpha) + np.power(-1*np.log(v),alpha), -2.0 + 2.0/alpha)*np.power(np.log(u)*np.log(v),alpha-1.0)
        p2 = 1.0 + (alpha - 1.0)*np.power(np.power(-1*np.log(u),alpha) + np.power(-1*np.log(v),alpha), -1.0/alpha )
        y = p1*p2
    return y

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
    nu = 2
    
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
            
        # plot the Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(gaussian_copula_pdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'Gaussian Copula PDF')
    
    elif(family.lower()=='t'):
        t_copula_pdf_python = copulapdf(family,U,Rho,nu)
        t_copula_pdf_matlab = matlab_data['t_copula_pdf'][:,0]
        
        # compare the two
        t_copula_test_result = np.allclose(t_copula_pdf_python,t_copula_pdf_matlab)
        if(t_copula_test_result):
            print 'T Copula Python calculation matches Matlab!'
        else:
            print 'T Copula Python calculation does NOT match Matlab!'
            
        # plot the Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(t_copula_pdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'T Copula PDF')
        
    elif(family.lower()=='clayton'):
        clayton_copula_pdf_python = copulapdf(family,U,alpha)
        clayton_copula_pdf_matlab = matlab_data['clayton_copula_pdf'][:,0]
        
        # compare the two
        clayton_copula_test_result = np.allclose(clayton_copula_pdf_python,clayton_copula_pdf_matlab)
        if(clayton_copula_test_result):
            print 'Clayton Copula Python calculation matches Matlab!'
        else:
            print 'Clayton Copula Python calculation does NOT match Matlab!'
            
        # plot the Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(clayton_copula_pdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'Clayton Copula PDF')
        
    elif(family.lower()=='gumbel'):
        alpha = 1.5
        gumbel_copula_pdf_python = copulapdf(family,U,alpha)
        gumbel_copula_pdf_matlab = matlab_data['gumbel_copula_pdf'][:,0]
        
        # compare the two
        gumbel_copula_test_result = np.allclose(gumbel_copula_pdf_python,gumbel_copula_pdf_matlab)
        if(gumbel_copula_test_result):
            print 'Gumbel Copula Python calculation matches Matlab!'
        else:
            print 'Gumbel Copula Python calculation does NOT match Matlab!'
            
        # plot the Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(gumbel_copula_pdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'Gumbel Copula PDF')
        
    elif(family.lower()=='frank'):
        frank_copula_pdf_python = copulapdf(family,U,alpha)
        frank_copula_pdf_matlab = matlab_data['frank_copula_pdf'][:,0]
        
        # compare the two
        frank_copula_test_result = np.allclose(frank_copula_pdf_python,frank_copula_pdf_matlab)
        if(frank_copula_test_result):
            print 'Frank Copula Python calculation matches Matlab!'
        else:
            print 'Frank Copula Python calculation does NOT match Matlab!'
            
        # plot the Copula for fun
        X = UU[0]
        Y = UU[1]
        Z = np.reshape(frank_copula_pdf_python,UU[0].shape)
        
        plot_utils.plot_3d(X,Y,Z, 'Frank Copula PDF')
    
if __name__=='__main__':
    import scipy.io
    import plot_utils
    
    test_python_vs_matlab('Gaussian')
    test_python_vs_matlab('T')
    test_python_vs_matlab('Clayton')
    test_python_vs_matlab('Gumbel')
    test_python_vs_matlab('Frank')