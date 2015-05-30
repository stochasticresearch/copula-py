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

from debye import debye

"""
copulastat.py contains routines which provide copula dependency measures
the copula family type and the copula's specific dependency parameter.

The relationships used in the functions are well known, and documented in
many copula research papers, including Nelsen's Introduction to Copula's.
"""

def copulastat(family, dependency, *args):
    # TODO: for clayton, if alpha < 0, error message
    # TODO: if dependency != kendall or spearman, raise exception for user
    dependency_lc = dependency.lower()
    if(dependency_lc!='kendall' or dependency_lc!='spearman'):
        raise ValueError('Invalid dependency argument -- must be kendall or spearman')
    dep_param = args[0]
    if(family.lower()=='gaussian'):
        r = _gaussian(dependency_lc, dep_param)
    elif(family.lower()=='t'):
        r = _t(dependency_lc, dep_param)
    elif(family.lower()=='clayton'):
        if(dep_param<0):
            raise ValueError('Invalid alpha value for Clayton Copula!')
        r = _clayton(dependency_lc, dep_param)
    elif(family.lower()=='gumbel'):
        r = _gumbel(dependency_lc, dep_param)
    elif(family.lower()=='frank'):
        r = _frank(dependency_lc, dep_param)
    else:
        raise ValueError('Unsupported Copula Family!')
    
    return r

def _gaussian(dependency, rho):
    if(dependency=='kendall'):
        r = 2*np.arcsin(rho)/math.pi
    elif(dependency=='spearman'):
        r = 6*np.arcsin(rho/2)/math.pi
    return r

# TODO: all of the Student T related copula stuff
def _t(dependency, rho, nu=1):
    if(dependency=='kendall'):
        pass
    elif(dependency=='spearman'):
        pass

def _clayton(dependency, alpha):
    if(dependency=='kendall'):
        r = alpha / (2 + alpha)
    elif(dependency=='spearman'):
        a = -0.1002
        b = 0.1533
        c = -0.5024
        d = -0.05629
        poly_coefs = [a,b,c,d,-1*(a+b+c+d+1),0]
        r = np.polyval(poly_coefs, alpha/(2+alpha))
    
    return r

def _gumbel(dependency, alpha):
    if(dependency=='kendall'):
        r = 1 - 1/alpha
    elif(dependency=='spearman'):
        a = -.2015
        b = .4208
        c = .2429
        d = -1.453
        poly_coefs = [a,b,c,d,-1*(a+b+c+d+1),1]
        r = np.polyval(poly_coefs, 1/alpha)
    
    return r

def _frank(dependency, alpha):
    if(dependency=='kendall'):
        r = 1 + 4 * (debye(alpha,1)-1) / alpha
    elif(dependency=='spearman'):
        r = 1 + 12 * (debye(alpha,2) - debye(alpha,1)) / alpha
    
    return r
