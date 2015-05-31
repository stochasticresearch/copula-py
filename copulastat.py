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
    dependency_lc = dependency.lower()
    if(dependency_lc!='kendall' and dependency_lc!='spearman'):
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
        poly_coefs = [a,b,c,d,-1*(a+b+c+d-1),0]
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

def test_python_vs_matlab(family):
    # test the python data against Matlab
    # TODO: make python execute the matlab script which generates these samples
    matlab_data = scipy.io.loadmat('matlab/copulastat_test.mat')
    
    if(family.lower()=='gaussian'):
        rho = 0.3
        gauss_ktau_rho_0_3_python = copulastat(family,'kendall',rho)
        gauss_srho_rho_0_3_python = copulastat(family,'spearman',rho)
        rho = 0.7
        gauss_ktau_rho_0_7_python = copulastat(family,'kendall',rho)
        gauss_srho_rho_0_7_python = copulastat(family,'spearman',rho)
        rho = 1.0
        gauss_ktau_rho_1_0_python = copulastat(family,'kendall',rho)
        gauss_srho_rho_1_0_python = copulastat(family,'spearman',rho)
        
        p1 = np.isclose(gauss_ktau_rho_0_3_python, matlab_data['gauss_ktau_rho_0_3'])
        p2 = np.isclose(gauss_srho_rho_0_3_python, matlab_data['gauss_srho_rho_0_3'])
        p3 = np.isclose(gauss_ktau_rho_0_7_python, matlab_data['gauss_ktau_rho_0_7'])
        p4 = np.isclose(gauss_srho_rho_0_7_python, matlab_data['gauss_srho_rho_0_7'])
        p5 = np.isclose(gauss_ktau_rho_1_0_python, matlab_data['gauss_ktau_rho_1_0'])
        p6 = np.isclose(gauss_srho_rho_1_0_python, matlab_data['gauss_srho_rho_1_0'])
        
        if(p1 and p2 and p3 and p4 and p5 and p6):
            print 'Gaussian CopulaStat tests PASSED!'
        else:
            print 'Gaussian CopulaStat tests FAILED!'
    elif(family.lower()=='t'):
        pass
    elif(family.lower()=='clayton'):
        alpha = 0.3
        clayton_ktau_alpha_0_3_python = copulastat(family,'kendall',alpha)
        clayton_srho_alpha_0_3_python = copulastat(family,'spearman',alpha)
        alpha = 0.7
        clayton_ktau_alpha_0_7_python = copulastat(family,'kendall',alpha)
        clayton_srho_alpha_0_7_python = copulastat(family,'spearman',alpha)
        alpha = 1.0
        clayton_ktau_alpha_1_0_python = copulastat(family,'kendall',alpha)
        clayton_srho_alpha_1_0_python = copulastat(family,'spearman',alpha)
        
        p1 = np.isclose(clayton_ktau_alpha_0_3_python, matlab_data['clayton_ktau_alpha_0_3'])
        p2 = np.isclose(clayton_srho_alpha_0_3_python, matlab_data['clayton_srho_alpha_0_3'])
        p3 = np.isclose(clayton_ktau_alpha_0_7_python, matlab_data['clayton_ktau_alpha_0_7'])
        p4 = np.isclose(clayton_srho_alpha_0_7_python, matlab_data['clayton_srho_alpha_0_7'])
        p5 = np.isclose(clayton_ktau_alpha_1_0_python, matlab_data['clayton_ktau_alpha_1_0'])
        p6 = np.isclose(clayton_srho_alpha_1_0_python, matlab_data['clayton_srho_alpha_1_0'])
        
        if(p1 and p2 and p3 and p4 and p5 and p6):
            print 'Clayton CopulaStat tests PASSED!'
        else:
            print 'Clayton CopulaStat tests FAILED!'
    elif(family.lower()=='gumbel'):
        alpha = 1.0
        gumbel_ktau_alpha_1_0_python = copulastat(family,'kendall',alpha)
        gumbel_srho_alpha_1_0_python = copulastat(family,'spearman',alpha)
        alpha = 3.0
        gumbel_ktau_alpha_3_0_python = copulastat(family,'kendall',alpha)
        gumbel_srho_alpha_3_0_python = copulastat(family,'spearman',alpha)
        
        p1 = np.isclose(gumbel_ktau_alpha_1_0_python, matlab_data['gumbel_ktau_alpha_1_0'])
        p2 = np.isclose(gumbel_srho_alpha_1_0_python, matlab_data['gumbel_srho_alpha_1_0'])
        p3 = np.isclose(gumbel_ktau_alpha_3_0_python, matlab_data['gumbel_ktau_alpha_3_0'])
        p4 = np.isclose(gumbel_srho_alpha_3_0_python, matlab_data['gumbel_srho_alpha_3_0'])
        
        if(p1 and p2 and p3 and p4):
            print 'Gumbel CopulaStat tests PASSED!'
        else:
            print 'Gumbel CopulaStat tests FAILED!'
    elif(family.lower()=='frank'):
        alpha = 0.3
        frank_ktau_alpha_0_3_python = copulastat(family,'kendall',alpha)
        frank_srho_alpha_0_3_python = copulastat(family,'spearman',alpha)
        alpha = 0.7
        frank_ktau_alpha_0_7_python = copulastat(family,'kendall',alpha)
        frank_srho_alpha_0_7_python = copulastat(family,'spearman',alpha)
        alpha = 1.0
        frank_ktau_alpha_1_0_python = copulastat(family,'kendall',alpha)
        frank_srho_alpha_1_0_python = copulastat(family,'spearman',alpha)
        
        p1 = np.isclose(frank_ktau_alpha_0_3_python, matlab_data['frank_ktau_alpha_0_3'])
        p2 = np.isclose(frank_srho_alpha_0_3_python, matlab_data['frank_srho_alpha_0_3'])
        p3 = np.isclose(frank_ktau_alpha_0_7_python, matlab_data['frank_ktau_alpha_0_7'])
        p4 = np.isclose(frank_srho_alpha_0_7_python, matlab_data['frank_srho_alpha_0_7'])
        p5 = np.isclose(frank_ktau_alpha_1_0_python, matlab_data['frank_ktau_alpha_1_0'])
        p6 = np.isclose(frank_srho_alpha_1_0_python, matlab_data['frank_srho_alpha_1_0'])
        
        if(p1 and p2 and p3 and p4 and p5 and p6):
            print 'Frank CopulaStat tests PASSED!'
        else:
            print 'Frank CopulaStat tests FAILED!'
    
if __name__=='__main__':
    import scipy.io
    
    test_python_vs_matlab('Gaussian')
    test_python_vs_matlab('Clayton')
    test_python_vs_matlab('Gumbel')
    test_python_vs_matlab('Frank')
