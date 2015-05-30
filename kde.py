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
kde.py contains routines which help perform Kernel Density Estimation (KDE).
"""

def kde(x_i, kernel, h, n_points):
    """ Perform Kernel Density Estimation on a given set of points.

    Inputs:
    x_i -- the input data set, should be a numpy array
    kernel -- the kernel to use, must be a string of one of the choices:
                Uniform
                Triangular
                Epanechnikov
                Quartic
                Triweight
                Tricube
                Gaussian
                Cosine
                Logistic
                Silverman
    h -- the kernel bandwidth setting
    npoints -- the number of desired points in the kernel density estimate
     
    Outputs:
    y -- the kernel density estimate
    """
    # define the points over which we will generate the kernel density estimate
    x = np.linspace(min(x_i), max(x_i), n_points)
    n = x_i.size
    y = np.zeros(n_points)
    
    for ii in np.arange(n_points):
        # apply the kernel to the point of interest
        if(kernel.lower()=='uniform'):
            y[ii] = 1.0/(n*h) * np.sum(uniform_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='triangular'):
            y[ii] = 1.0/(n*h) * np.sum(triangle_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='epanechnikov'):
            y[ii] = 1.0/(n*h) * np.sum(epanechnikov_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='quartic'):
            y[ii] = 1.0/(n*h) * np.sum(quartic_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='triweight'):
            y[ii] = 1.0/(n*h) * np.sum(triweight_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='tricube'):
            y[ii] = 1.0/(n*h) * np.sum(tricube_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='gaussian'):
            y[ii] = 1.0/(n*h) * np.sum(gaussian_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='cosine'):
            y[ii] = 1.0/(n*h) * np.sum(cosine_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='logistic'):
            y[ii] = 1.0/(n*h) * np.sum(logistic_kernel( (x[ii]-x_i)/h ) )
        elif(kernel.lower()=='silverman'):
            y[ii] = 1.0/(n*h) * np.sum(silverman_kernel( (x[ii]-x_i)/h ) )
        else:
            print 'In here:)'

    return (x,y)
    
def uniform_kernel(u):
    """
    %UNIFORM_KDE - the uniform kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = 1.0/2.0
    
    return y

def triangle_kernel(u):
    """
    %TRIANGLE_KDE - the triangular kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = 1.0-abs(u[idxs[0]])
    
    return y

def epanechnikov_kernel(u):
    """
    %EPANECHNIKOV_KDGE - the epanechnikov kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = 3.0/4.0*(1-np.power(u[idxs[0]],2))
    
    return y

def quartic_kernel(u):
    """
    %QUARTIC_KDE - the quartic kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = 15.0/16.0*np.power((1-np.power(u[idxs[0]],2)),2)
    
    return y

def triweight_kernel(u):
    """
    %QUARTIC_KDE - the triweight kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = 35.0/32.0*np.power((1-np.power(u[idxs[0]],2)),3)
    
    return y

def tricube_kernel(u):
    """
    %QUARTIC_KDE - the quartic kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = 70.0/81.0*np.power((1-np.power(abs(u[idxs[0]]),3)),3)
    
    return y

def gaussian_kernel(u):
    """
    %GAUSSIAN_KDE - the gaussian kernel
    """
    y = 1.0/math.sqrt(2*math.pi) * np.exp(-np.power(u,2)/2.0)
    
    return y

def cosine_kernel(u):
    """
    %COSINE_KDE - the cosine kernel
    """
    idxs = np.where(abs(u)<=1)
    y = np.zeros(u.size)
    y[idxs[0]] = math.pi/4.0*np.cos(math.pi/2.0*u[idxs[0]])
    
    return y

def logistic_kernel(u):
    """
    %LOGISTIC_KDE - the logistic kernel
    """
    y = 1.0/(np.exp(u) + 2.0 + np.exp(-u))
    
    return y

def silverman_kernel(u):
    """
    %SILVERMAN_KDE - the silverman kernel
    """
    y = 1.0/2.0 * np.exp(-abs(u)/math.sqrt(2)) * np.sin(abs(u)/2 + math.pi/4)
    
    return y


if __name__=='__main__':
    import matplotlib.pyplot as plt
    import sys
    
    # TODO: put in argument to allow user to test windows and change
    # the if from if false to that if condition
    
    if(False):
        # Plot the Uniform Kernel
        x_i = np.linspace(-2,2,100)
        y = uniform_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Uniform Kernel')
        plt.show()
        
        # Plot the Triangle Kernel
        x_i = np.linspace(-2,2,100)
        y = triangle_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Triangle Kernel')
        plt.show()
        
        # Plot the Epanechnikov Kernel
        x_i = np.linspace(-2,2,100)
        y = epanechnikov_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Epanechnikov Kernel')
        plt.show()
        
        # Plot the Quartic Kernel
        x_i = np.linspace(-2,2,100)
        y = quartic_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Quartic Kernel')
        plt.show()
        
        # Plot the Triweight Kernel
        x_i = np.linspace(-2,2,100)
        y = triweight_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Triweight Kernel')
        plt.show()
        
        # Plot the Tricube Kernel
        x_i = np.linspace(-2,2,100)
        y = tricube_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Tricube Kernel')
        plt.show()
        
        # Plot the Gaussian Kernel
        x_i = np.linspace(-2,2,100)
        y = gaussian_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Gaussian Kernel')
        plt.show()
        
        # Plot the Cosine Kernel
        x_i = np.linspace(-2,2,100)
        y = cosine_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Cosine Kernel')
        plt.show()
        
        # Plot the Logistic Kernel
        x_i = np.linspace(-2,2,100)
        y = logistic_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Logistic Kernel')
        plt.show()
        
        # Plot the Silverman Kernel
        x_i = np.linspace(-2,2,100)
        y = silverman_kernel(x_i)
        plt.plot(x_i,y)
        plt.title('Silverman Kernel')
        plt.show()
    
    # test the KDE estimation
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
    
    # the histogram of the data
    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
    
    (xx,pp) = kde(x,'Gaussian',h, npoints)
    plt.plot(xx,pp, 'b')
    plt.title('Kernel Density Estimate')
    plt.show()