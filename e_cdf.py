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

"""
e_cdf.py contains routines which help perform empirical CDF Estimation.
"""

def e_cdf(x_i, npoints):
    """ Generates an Empirical CDF using the indicator function.
    
    Inputs:
    x_i -- the input data set, should be a numpy array
    npoints -- the number of desired points in the empirical CDF estimate
     
    Outputs:
    y -- the empirical CDF
    """
    # define the points over which we will generate the kernel density estimate
    x = np.linspace(min(x_i), max(x_i), npoints)
    n = float(x_i.size)
    y = np.zeros(npoints)
    
    for ii in np.arange(x.size):
        idxs = np.where(x_i<=x[ii])
        y[ii] = np.sum(idxs[0].size)/n
    
    return (x,y)

def kde_integral(kde):
    """ Generates a "smoother" Empirical CDF by integrating the KDE.  For this,
        the user should first generate the KDE using kde.py, and then pass the
        density estimate to this function
        
        Inputs:
        kde -- the kernel density estimate
        
        Outputs:
        y -- the smoothed CDF estimate
    """
    y = np.cumsum(kde)/sum(kde)
    
    return y

if __name__=='__main__':
    import matplotlib.pyplot as plt
    import sys
    import kde
    
    # test the E_CDF estimation
    N1 = 100 # number of data in data set 1
    m1 = -1  # mean value
    s1 = 0.1 # % variance 

    N2 = 500 # number of data in data set 2
    m2 = 2   # mean value
    s2 = 0.5 # variance 
    
    h = 0.1       # bandwidth
    npoints = 100 # number of abscis points in kde

    x1 = math.sqrt(s1)*np.random.randn(N1,1) + m1
    x2 = math.sqrt(s2)*np.random.randn(N2,1) + m2
    x = np.concatenate((x1,x2),axis=0)
    
    # Kernel Density Estimate
    (xx,kde_estimate) = kde.kde(x,'Gaussian',h, npoints)
    plt.plot(xx,kde_estimate, 'r', label='Kernel Density Estimate')
    
    # the histogram of the data
    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75, label='Histogram')
    
    # empirical CDF
    (xx,pp) = e_cdf(x, npoints)
    plt.plot(xx,pp, 'k', label='Empirical CDF')
    
    # Smooth Empirical CDF (KDE Integral)
    kde_integral = kde_integral(kde_estimate)
    plt.plot(xx,kde_integral, 'm', label='Smooth Empirical CDF')
    plt.legend(loc='upper left')
    plt.show()