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

from scipy.stats import spearmanr
from scipy.stats import kendalltau
from scipy.misc  import comb

"""
Encompasses calculation of Spearman's Rho and Kendall's Tau (and other statistical
measures to be added in the future) for data with dimensionality >= 2
"""

def spearmans_rho(X):
    """
    Calculates a generalized Spearman's rho for a data set given by X, as 
    described by "Multivariate Extensions of Spearman's Rho and Related Statistics"
    Inputs:
      X - the input data, should be a numpy array of shape = M x N, where
          M is the number of samples, and N is the dimensionality of the data
    """
    M = X.shape[0]
    N = X.shape[1]
    if N<2:
        raise ValueError('To calculate Spearman\'s Rho, need data of dimensionality >= 2')
    
    srho = 0.0
    for dim1 in range(0,N-1):
        for dim2 in range(dim1+1,N):
            (r,p) = spearmanr(X[:,dim1],X[:,dim2])
            srho = srho + r
    # normalize
    srho = srho / comb(N,2)
    return srho
    
def kendalls_tau(X):
    """
    Calculates a generalized Kendall's tau for a data set given by X, as 
    described by "Multivariate Extensions of Spearman's Rho and Related Statistics"
    
    Inputs:
      X - the input data, should be a numpy array of shape = M x N, where
          M is the number of samples, and N is the dimensionality of the data
    """
    M = X.shape[0]
    N = X.shape[1]
    if N<2:
        raise ValueError('To calculate Kendall\'s Tau, need data of dimensionality >= 2')
    
    ktau = 0.0
    for dim1 in range(0,N-1):
        for dim2 in range(dim1+1,N):
            (t,p) = kendalltau(X[:,dim1],X[:,dim2])
            ktau = ktau + t
    # normalize
    ktau = ktau / comb(N,2)
    return ktau

if __name__=='__main__':
    X = np.array([[12,1,-3],
                  [2,4,-4],
                  [1,7,-6],
                  [12,1,2],
                  [2,0,1]])
    srho = spearmans_rho(X)
    ktau = kendalls_tau(X)
    
    print srho, ktau