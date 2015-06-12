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

from invcopulastat import invcopulastat
from copulacdf import copulacdf

def cvolume(family, u1v1, u1v2, u2v1, u2v2, *args):
    """
    Computes the C-Volume of a specified copula family with dependency parameter
    defined in the args.
      family - the copula type, must be:
        'Gaussian'
        'T'
        'Clayton'
        'Frank'
        'Gumbel'
      u1v1 - a N x 2 matrix of values between [0,1] that represents the bottom
             left coordinate of the grid for which the C-Volume is desired
      u1v2 - a N x 2 matrix of values between [0,1] that represent the top
             left coordinate of the grid for which the C-Volume is desired
      u2v1 - a N x 2 matrix of values between [0,1] that represent the bottom
             right coordinate of the grid for which the C-volume is desired
      u2v2 - a N x 2 matrix of values between [0,1] that represents the top
             right coordinate of the grid for which the C-Volume is desired
      args - must be atleast of length 2, for which the first element in args
             is expected to be a string which describes the dependency value
             being provided, must be one of the following:
        'kendall' - means kendall's Tau is being provided
        'spearman' - means spearman's rho is being provided
        'native' - means that the dependency parameter of the copula family
                   itself is being provided directly
             the second argmuent  must be the value of the dependency type 
            provided. For kendall and spearman, a scalar value is expected.  
            For native, if the family type is Frank, Gumbel, or Clayton, then 
            a scalar value is expected, which represents the dependency
            parameter.  If the family type is Gaussian, then a 2 x 2 numpy array
            is expected, which represents the correlation matrix defining the
            Gaussian copula.  If the family is T, then the 2nd argument is the
            2x2 numpy array representing the correlation matrix, and the 3rd
            argument is the degrees of freedom
    """    
    family_lc = family.lower()
    if(family_lc=='gaussian'):
        if(len(args)<2):
            raise ValueError("Gaussian Family expects 2 variable arguments, the dependency type and value")
        if(args[0]=='kendall' or args[0]=='spearman'):
            # get the correlation parameter
            r = invcopulastat(family, args[0], args[1])
        else:
            r = args[1]
        
        cvol = _gaussian(u1v1, u1v2, u2v1, u2v2, r)
    elif(family_lc=='t'):
        if(len(args)<2):
            raise ValueError("T Family expects atleast 2 variable arguments, the dependency type and value")
        
        if(args[0]=='kendall' or args[0]=='spearman'):
            raise ValueError('T Family does not accept Kendalls Tau or Spearmans Rho, only native parameters')
        else:
            r = args[1]
            nu = args[2]
            
            cvol = _gaussian(u1v1, u1v2, u2v1, u2v2, r, nu)
            
    elif(family_lc=='clayton'):
        if(len(args)<2):
            raise ValueError("Clayton Family expects 2 variable arguments, the dependency type and value")
        
        if(args[0]=='kendall' or args[0]=='spearman'):
            # get the correlation parameter and degrees of freedom
            alpha = invcopulastat(family, args[0], args[1])
        else:
            alpha = args[1]
        
        cvol = _clayton(u1v1, u1v2, u2v1, u2v2, alpha)
        
    elif(family_lc=='frank'):
        if(len(args)<2):
            raise ValueError("Frank Family expects 2 variable arguments, the dependency type and value")
        if(args[0]=='kendall' or args[0]=='spearman'):
            # get the correlation parameter and degrees of freedom
            alpha = invcopulastat(family, args[0], args[1])
        else:
            alpha = args[1]
        
        cvol = _frank(u1v1, u1v2, u2v1, u2v2, alpha)

    elif(family_lc=='gumbel'):
        if(len(args)<2):
            raise ValueError("Gumbel Family expects 2 variable arguments, the dependency type and value")
        if(args[0]=='kendall' or args[0]=='spearman'):
            # get the correlation parameter and degrees of freedom
            alpha = invcopulastat(family, args[0], args[1])
        else:
            alpha = args[1]
        
        cvol = _gumbel(u1v1, u1v2, u2v1, u2v2, alpha)

    return cvol

    
def _gaussian(u1v1, u1v2, u2v1, u2v2, r):
    # generate the Rho matrix from r
    Rho = np.ones((2,2))
    Rho[0][1] = r
    Rho[1][0] = r
    
    # this is the equation for C Volume as defined by Nelsen
    cvol = copulacdf('Gaussian', u2v2, Rho) - \
           copulacdf('Gaussian', u2v1, Rho) - \
           copulacdf('Gaussian', u1v2, Rho) + \
           copulacdf('Gaussian', u1v1, Rho) 
    
    return cvol

def _t(u1v1, u1v2, u2v1, u2v2, r, nu):
    # generate the Rho matrix from r
    Rho = np.ones((2,2))
    Rho[0][1] = r
    Rho[1][0] = r
    
    # this is the equation for C Volume as defined by Nelsen
    cvol = copulacdf('T', u2v2, Rho, nu) - \
           copulacdf('T', u2v1, Rho, nu) - \
           copulacdf('T', u1v2, Rho, nu) + \
           copulacdf('T', u1v1, Rho, nu) 
    
    return cvol
    
    return None

def _clayton(u1v1, u1v2, u2v1, u2v2, alpha):
    
    # this is the equation for C Volume as defined by Nelsen
    cvol = copulacdf('Clayton', u2v2, alpha) - \
           copulacdf('Clayton', u2v1, alpha) - \
           copulacdf('Clayton', u1v2, alpha) + \
           copulacdf('Clayton', u1v1, alpha) 
    
    return cvol

def _frank(u1v1, u1v2, u2v1, u2v2, alpha):
    
    # this is the equation for C Volume as defined by Nelsen
    cvol = copulacdf('Frank', u2v2, alpha) - \
           copulacdf('Frank', u2v1, alpha) - \
           copulacdf('Frank', u1v2, alpha) + \
           copulacdf('Frank', u1v1, alpha) 
    
    return cvol

def _gumbel(u1v1, u1v2, u2v1, u2v2, alpha):
    
    # this is the equation for C Volume as defined by Nelsen
    cvol = copulacdf('Gumbel', u2v2, alpha) - \
           copulacdf('Gumbel', u2v1, alpha) - \
           copulacdf('Gumbel', u1v2, alpha) + \
           copulacdf('Gumbel', u1v1, alpha) 
    
    return cvol
