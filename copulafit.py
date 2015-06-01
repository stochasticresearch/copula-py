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
        'IFM' - Inference Functions for Margin's method
        'CML' - Canonical Maximum Likelihood
        'Dependency' - Use's Kendall's Tau or Spearman's Rho relationships to the 
                       copula family's dependency parameter 
                       
    Outputs:
      the dependency parameter for the copula
                       
    """
    algorithm_lc = algorithm.lower()
    family_lc    = family.lower()
    if(algorithm_lc=='ifm'):
        raise Exception('IFM method not yet supported!')
    elif(algorithm_lc=='cml'):
        raise Exception('CML method not yet supported!')
    elif(algorithm_lc=='dependency'):
        if(family_lc=='gaussian'):
            _gauss_dependency(X)
        elif(family_lc=='t'):
            _t_dependency(X)
        elif(family_lc=='clayton'):
            _clayton_dependency(X)
        elif(family_lc=='gumbel'):
            _gumbel_dependency(X)
        elif(family_lc=='frank'):
            _frank_dependency(X)
    else:
        raise Exception('Unsupported Algorithm!')
    
def _gauss_dependency(X):
    pass
    
def _t_dependency(X):
    pass

def _clayton_dependency(X):
    pass

def _gumbel_dependency(X):
    pass

def _frank_dependency(X):
    pass