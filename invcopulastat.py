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
from scipy.optimize import fsolve

from copulastat import copulastat

"""
invcopulastat.py contains routines which provide the inverse copula dependency 
measures the copula family type and the copula's specific dependency parameter.

The relationships used in the functions are well known, and documented in
many copula research papers, including Nelsen's Introduction to Copula's.
"""

def invcopulastat(family, dependency, val):
    dependency_lc = dependency.lower()
    if(dependency_lc!='kendall' and dependency_lc!='spearman'):
        raise ValueError('Invalid dependency argument -- must be kendall or spearman')
    if(family.lower()=='gaussian'):
        r = _gaussian(dependency_lc, val)
    elif(family.lower()=='t'):
        r = None
    elif(family.lower()=='clayton'):
        r = _clayton(dependency_lc, val)
    elif(family.lower()=='gumbel'):
        r = _gumbel(dependency_lc, val)
    elif(family.lower()=='frank'):
        r = _frank(dependency_lc, val)
    else:
        raise ValueError('Unsupported Copula Family!')
    
    return r

def _gaussian(dependency, val):
    if(dependency=='kendall'):
        r = np.sin(math.pi/2.0*val)
    elif(dependency=='spearman'):
        r = 2*np.sin(math.pi/6.0*val)
    return r

# TODO: all studnet-t related stuff
def _t(dependency, val):
    return None

def _clayton(dependency, val):
    if(dependency=='kendall'):
        d = 2.0*val/(1.0-val)
    elif(dependency=='spearman'):
        raise ValueError('Spearmans Rho currently unsupported for Clayton Copula family!')
    
    return d

def _gumbel(dependency, val):
    if(dependency=='kendall'):
        d = 1.0/(1.0-val)
    elif(dependency=='spearman'):
        raise ValueError('Spearmans Rho currently unsupported for Gumbel Copula family!')
    
    return d

def _frank_kendall_fopt(alpha, tau):
    return 4*( debye(alpha,1) - 1 )/alpha + 1 - tau

def _frank(dependency, val):
    if(dependency=='kendall'):
        return fsolve(_frank_kendall_fopt, 1, args=(val))
    elif(dependency=='spearman'):
        # TODO --  use function solvers in scipy to invert debye function for the closed form solution
        raise ValueError('Spearmans Rho currently unsupported for Frank Copula family!')
    
    return r

def test_python_vs_matlab(family):
    # DISCLAIMER: this code assumes copulastat is working properly and tested
    
    if(family.lower()=='gaussian'):
        dependency = 'kendall'
        rho = 0.3
        rho_calc = invcopulastat(family, dependency, copulastat(family, dependency, rho))
        p1 = np.isclose(rho, rho_calc)
        
        rho = 0.7
        rho_calc = invcopulastat(family, dependency, copulastat(family, dependency, rho))
        p2 = np.isclose(rho, rho_calc)
        
        dependency = 'spearman'
        rho = 0.3
        rho_calc = invcopulastat(family, dependency, copulastat(family, dependency, rho))
        p3 = np.isclose(rho, rho_calc)
        
        rho = 0.7
        rho_calc = invcopulastat(family, dependency, copulastat(family, dependency, rho))
        p4 = np.isclose(rho, rho_calc)
        
        if(p1 and p2 and p3 and p4):
            print 'Gaussian CopulaStat tests PASSED!'
        else:
            print 'Gaussian CopulaStat tests FAILED!'
        
    elif(family.lower()=='t'):
        pass
    
    elif(family.lower()=='clayton' or family.lower()=='gumbel' or family.lower()=='frank'):
        dependency = 'kendall'
        alpha = 0.3
        tau = copulastat(family, dependency, alpha)
        alpha_calc = invcopulastat(family, dependency, tau)
        p1 = np.isclose(alpha, alpha_calc)
        
        alpha = 0.7
        tau = copulastat(family, dependency, alpha)
        alpha_calc = invcopulastat(family, dependency, tau)
        p2 = np.isclose(alpha, alpha_calc)
        
        if(p1 and p2):
            print family + ' CopulaStat tests PASSED!'
        else:
            print family + ' CopulaStat tests FAILED!'
    
if __name__=='__main__':    
    test_python_vs_matlab('Gaussian')
    test_python_vs_matlab('Clayton')
    test_python_vs_matlab('Gumbel')
    test_python_vs_matlab('Frank')