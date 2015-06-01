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

import multivariate_stats
from invcopulastat import invcopulastat
from scipy.stats import kendalltau

"""
copulafit.py contains routines which provide use various techniques, as
specified by the user to fit data to a family of copula (i.e. find the
dependency parameter).
"""

def copulafit(family, X, algorithm):
    """
    Attempts to determine the dependency parameter of the copula family
    type specified, using the algorithm that is specified for the data
    given by the matrix X
    
    Inputs:
      family -- the copula family to fit to, must be:
        'Gaussian'
        't'
        'Clayton'
        'Gumbel'
        'Frank'
      X -- the data to determine the copula dependency parameter for. Must be
           a numpy array of shape = M x N, where M is the number of samples 
           and N is the dimensionality of the data
      algorithm -- must be one of the following strings:
        'MLE'  - Maximum Likelihood method
        'AMLE' - Approximate Maximum Likelihood method
        'PKTE' - Use's Pairwise Kendall's Tau estimator's relationship to the 
                       copula family's dependency parameter (only applicalble
                       to Clayton, Gumbel, or Frank copula's currently)
                       
    Outputs:
      the dependency parameter for the copula
                       
    """
    algorithm_lc = algorithm.lower()
    family_lc    = family.lower()
    dep_param    = None
    if(algorithm_lc=='MLE'):
        raise ValueError('MLE method not yet supported!')
    elif(algorithm_lc=='AMLE'):
        raise ValueError('Approximate MLE method not yet supported!')
    elif(algorithm_lc=='PKTE'):
        if(family_lc=='gaussian'):
            dep_param = _gaussian_PKTE(X)
        elif(family_lc=='t'):
            dep_param = _t_PKTE(X)
        elif(family_lc=='clayton'):
            dep_param = _clayton_PKTE(X)
        elif(family_lc=='gumbel'):
            dep_param = _gumbel_PKTE(X)
        elif(family_lc=='frank'):
            dep_param = _frank_PKTE(X)
    else:
        raise ValueError('Unsupported Algorithm or options!')
    
    return dep_param
    
def _guassian_PKTE(X):
    # the algorithm for this comes from the paper:
    # "Gaussian Copula Precision Estimation with Missing Values" 
    # by Huahua Wang, Faridel Fazayeli, Soumyadeep Chatterjee, Arindam Banerjee
    N = X.shape[1]
    sigma_hat = np.ones((N,N))
    for dim1 in range(0,N-1):
        for dim2 in range(dim1+1,N):
            rho = np.sin(math.pi/2 * kendalltau(X[:,dim1],X[:,dim2]))
            # correlation matrix is symmetric
            sigma_hat[dim1][dim2] = rho
            sigma_hat[dim2][dim1] = rho

# TODO: T copula stuff
def _t_PKTE(X):
    return None
    
def _clayton_PKTE(X):
    # calculate empirical kendall's tau
    ktau = multivariate_stats.kendalls_tau(X)
    # inverse to find dependency parameter
    alpha = invcopulastat('Clayton', 'kendall', ktau)
    
    return alpha

def _gumbel_PKTE(X):
    # calculate empirical kendall's tau
    ktau = multivariate_stats.kendalls_tau(X)
    # inverse to find dependency parameter
    alpha = invcopulastat('Gumbel', 'kendall', ktau)
    
    return alpha


def _frank_PKTE(X):
    # calculate empirical kendall's tau
    ktau = multivariate_stats.kendalls_tau(X)
    # inverse to find dependency parameter
    alpha = invcopulastat('Frank', 'kendall', ktau)
    
    return alpha
