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

import numpy as np
import scipy.integrate as integrate

def debye(x, n):
    """
    Evaluate the Debye function.
    See http://en.wikipedia.org/wiki/Debye_function for details
    """
    
    # ensure n is a float
    n = float(n)
    
    sol = integrate.quad( lambda xi: pow(xi,n)/(np.exp(xi)-1.0) , 0.0, x)
    return n*sol[0]/pow(x,n)