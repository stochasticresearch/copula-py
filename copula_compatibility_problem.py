#!/usr/bin/env python

import math
import numpy as np

from scipy.stats import mvn                     # contains inverse CDF of Multivariate Gaussian
from scipy.stats import norm                    # contains PDF of Gaussian

import copulacdf
import plot_utils
import matplotlib.pyplot as plt

"""
This file showcases what is known as the compatibility problem with copula's.
"""


"""
Here, we showcase the compatibility problem by using a 3-Copula and calculating
the two marginal's.
"""
def ex1():
    n = 10
    eps = np.finfo(float).eps
    u = np.linspace(0+eps,1-eps,n)
    UU = np.meshgrid(u,u)
    U2 = np.reshape(UU[0], (UU[0].shape[0]*UU[0].shape[1], 1))
    U3 = np.reshape(UU[1], (UU[1].shape[0]*UU[1].shape[1], 1))
    U1 = np.ones(U2.shape)*(1-eps)
    U = np.concatenate((U1,U2,U3),axis=1)
    
    R1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
    R2 = np.array([[1,0.6,-0.3],[0.6,1,0.4],[-0.3,0.4,1]])
    
    C1_twomarginal = copulacdf.copulacdf('Gaussian',U,R1)
    C2_twomarginal = copulacdf.copulacdf('Gaussian',U,R2)
    
    #X = UU[0]
    #Y = UU[1]
    #Z = np.reshape(C1_twomarginal,UU[0].shape)
    #plot_utils.plot_3d(X,Y,Z, 'C1 Two-Marginal')
    
    # compute error between C1_twomarginal and C2_twomarginal
    sq_err_vec = (C2_twomarginal-C1_twomarginal)**2
    
    X = UU[0]
    Y = UU[1]
    Z = np.reshape(sq_err_vec,UU[0].shape)
    plot_utils.plot_3d(X,Y,Z, 'Two-Marginal Error')
    
def ex2():
    n = 10
    eps = np.finfo(float).eps
    u = np.linspace(0+eps,1-eps,n)
    UU = np.meshgrid(u,u)
    U3 = np.reshape(UU[1], (UU[1].shape[0]*UU[1].shape[1], 1))
    U1 = np.ones(U3.shape)*(1-eps)
    U2 = np.ones(U3.shape)*(1-eps)
    U = np.concatenate((U1,U2,U3),axis=1)
    
    R1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
    R2 = np.array([[1,0.6,-0.3],[0.6,1,0.4],[-0.3,0.4,1]])
    
    C1_onemarginal = copulacdf.copulacdf('Gaussian',U,R1)
    C2_onemarginal = copulacdf.copulacdf('Gaussian',U,R2)
    sq_err_vec = (C2_onemarginal-C1_onemarginal)**2
    
    X = UU[0]
    Y = UU[1]
    Z = np.reshape(sq_err_vec,UU[0].shape)
    plot_utils.plot_3d(X,Y,Z, 'One-Margin Error')

if __name__=='__main__':
    ex1()
    ex2()