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

from scipy.stats import uniform
from scipy.stats import expon

"""
Algorithms copied directly from R source code of the copula package
    - rstable1.R
    - retstable.c
"""

# delta is assumed to be 0
def rstable1(n, alpha, beta, gamma=1, delta=0, pm=1):
    return _rstable_c(n, alpha) * gamma + delta

def _rstable_c(n, alpha):
    return np.power(np.cos(math.pi/2.0*alpha), -1.0/alpha) * _rstable0(alpha)
    
def _rstable0(alpha):
    U = uniform.rvs(size=1)
    while True:
        # generate non-zero exponential random variable
        W = expon.rvs(size=1)
        if(W!=0):
            break
    return np.power(_A(math.pi*U,alpha)/np.power(W,1.0-alpha),1.0/alpha)

def _A(x, alpha):
    Ialpha = 1.0-alpha
    return _A_3(x, alpha, Ialpha)

def _A_3(x, alpha, Ialpha):
    return np.power(Ialpha* np.sinc(Ialpha*x/math.pi), Ialpha) * \
            np.power(alpha * np.sinc(alpha *x/math.pi), alpha) / np.sinc(x/math.pi)
