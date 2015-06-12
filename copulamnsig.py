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

from cvolume import cvolume
import multivariate_stats

from ecdf import probability_integral_transform
from scipy.stats import entropy

def copulamnsig(family, K, *args):
    """
    Computes the copula multinomial signature as described in the paper
    "Highly Efficient Learning of Mixed Copula Networks" for a specified 
    copula family.  Essentially, it breaks up the unit grid into a K x K boxes, 
    and computes the probability of a sample from that copula pdf falling in 
    that grid.  This is then aggregated into a multinomial probability 
    distribution.  This so called "multinomial" signature of a copula is then 
    used to efficiently determine the structure of the Bayesian network, as well
    as the copula which would describe the dependency between the nodes.
    
    The grid over the unit cube is numbered as follows, for a 4 x 4 grid
      ___________________
      | 4 | 8 | 12 | 16 | 
      |---|---|----|----| 
      | 3 | 7 | 11 | 15 |
      |-----------------|
      | 2 | 6 | 10 | 14 |
      |-----------------|
      | 1 | 5 |  9 | 13 |
      |___|___|____|____|
    
    Currently, this computes the multinomial signature for a specified copula
    family of 2 dimensions.  It would be nice to expand this to multiple
    dimensions, and we can use the general formula for C-volume
    
      family - the copula type, must be:
        'Gaussian'
        'T'
        'Clayton'
        'Frank'
        'Gumbel'
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
    coords_list = _makeCoordsList(K)
            
    # mnsig is a list of dictionaries.  The (list index+1) corresponds to the
    # grid of interest in the unit cube.  In the dictionary, the actual lower
    # left coordinates of the box and the upper right coordinates of the box
    # are stored as keys 'u1v1' and 'u2v2', and then the actual value of the 
    # multinomial signature in that grid is stored as 'val'
    #mnsig = []
    
    mnsig = []
    for coord in coords_list:
        # compute the C-volume and store
        u1v1 = coord[0]
        u1v2 = coord[1]
        u2v1 = coord[2]
        u2v2 = coord[3]
        val = cvolume(family, u1v1, u1v2, u2v1, u2v2, *args)

        mnsig.append(val[0])
    
    return mnsig

def empirical_copulamnsig(X, K):
    """
    Computes an empirical copula multinomial signature based on the dataset
    provided by U.  U must be a numpy array of dimensions [M x N], where M is 
    the number of data points in the dataset and, N is the dimensionality of the
    data
    """
    M = X.shape[0]
    N = X.shape[1]
    
    # convert X to U by using the probability integral transform:  F(X) = U
    U = probability_integral_transform(X)
    
    # generate the coordinates so we can then compare and see where each sample
    # falls into in the unit cube
    coords_list = _makeCoordsList(K)
    
    # this will be a list of dictionaries which has all the combinations of the
    # empirical binomial signature 
    esig = []
    
    # for all i < j, compute pairwise bivariate multinomial signature
    for dim1 in range(0,N-1):
        for dim2 in range(dim1+1,N):
            # to compute the pairwise bivariate multinomial signature, what
            # we do is essentially grid as before, and compute a histogram 
            # for each grid .. whcih is our empirical estimate
            # the grid is lay-ed out in the exact same way as described before,
            # so the index of mnsig from copulamnsig and the index of the value
            # generated here will be directly comparable
            #     ___________________
            #     | 4 | 8 | 12 | 16 | 
            #     |---|---|----|----| 
            #     | 3 | 7 | 11 | 15 |
            #     |-----------------|
            #     | 2 | 6 | 10 | 14 |
            #     |-----------------|
            #     | 1 | 5 |  9 | 13 |
            #     |___|___|____|____|
            tmp = {}
            # RV 1 that we are comparing
            tmp['rv1'] = dim1+1
            # RV 2 that we are comparing
            tmp['rv2'] = dim2+1
            # the value for the zone -- initialize to 0
            esig_vec = np.zeros(K*K)
            
            # there is probably a more efficient way to do this than to loop
            # over each value, but this is a first cut at implementing this
            u = U[:,dim1]
            v = U[:,dim2]
            
            for ii in range(0,M):
                # find which zone this specific (u,v) sample falls in
                for jj in range(0,K*K):
                    u1 = coords_list[jj][0][0][0]
                    v1 = coords_list[jj][0][0][1]
                    u2 = coords_list[jj][3][0][0]
                    v2 = coords_list[jj][3][0][1]
                    
                    if(u[ii] >= u1 and u[ii] < u2 and 
                       v[ii] >= v1 and v[ii] < v2):
                        # add one to the zone that it falls into
                        esig_vec[jj] = (esig_vec[jj] + 1.0/M)
                        # process the next pair by kicking out of this loop
                        break
            tmp['esig'] = esig_vec
            
            esig.append(tmp)
        
    return esig

def _makeCoordsList(K):
    eps = np.finfo(float).eps
    u = np.linspace(0+eps, 1-eps, K+1)
    v = np.linspace(0+eps, 1-eps, K+1)
    
    coords_list = []
    for ii in range(0,len(u)-1):
        for jj in range(0,len(v)-1):
            u1 = u[ii]
            u2 = u[ii+1]
            v1 = v[jj]
            v2 = v[jj+1]
            u1v1 = np.array([[u1,v1]])
            u1v2 = np.array([[u1,v2]])
            u2v1 = np.array([[u2,v1]])
            u2v2 = np.array([[u2,v2]])
            x = []
            x.append(u1v1)
            x.append(u1v2)
            x.append(u2v1)
            x.append(u2v2)
            coords_list.append(x)
    
    return coords_list

############################## TODO ########################
# the master function, which computes the correct copula family to choose from
# will compare the empirical signatures to the actual signature for refernence
# will do the following:
#  1.) compute the empirical kendall's tau
#  2.) load the precomputed multinomial signature for that kendall's tau
#      for all the copula families
#  3.) minimize the distance metric
def optimalCopulaFamily(X, K=4, family_search=['Gaussian', 'Clayton', 'Gumbel', 'Frank']):
    """
    This function, given a multivariate data set X, computes the best copula family which fits
    the data, using the procedure described in the paper "Highly Efficient Learning of Mixed
    Copula Networks," by Gal Elidan
      
      X - the multivariate dataset for which we desire the copula.  Must be a numpy array of 
          dimension [M x N], where M is the number of data points, and N is the dimensionality
          of the dataset
      K - the square root of the number of grid points (for now, we assume square gridding of the
          unit cube)
      family_search - a list of all the copula families to search.  Currently, what is supported is
          Gaussian, Clayton, Gumbel, and Frank.  As more copula's are added, the default list will
          be expanded.
    """
    # compute the empirical Kendall's Tau
    tau_hat = multivariate_stats.kendalls_tau(X)
    
    # compute empirical multinomial signature
    empirical_mnsig = empirical_copulamnsig(X, K)
    empirical_mnsig = empirical_mnsig[0]['esig']
    
    # compute the multinomial signature for each of the copula families specified
    # and simultaneously compute the kullback leibler divergence between the empirical
    # and the computed, and store that info
    distances = {}
    for family in family_search:
        mnsig = copulamnsig(family,K,'kendall',tau_hat)
            
        # compute KL divergence, see
        # http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.stats.entropy.html
        distances[family] = entropy(mnsig, empirical_mnsig)
    
    # search for the minimum distance, that is the optimal copula family to use
    minDistance = 1000000
    for family, distance in distances.iteritems():
        if distance<minDistance:
            minDistance = distance
            optimalFamily = family
    
    depParams = invcopulastat(optimalFamily, 'kendall', tau_hat)
    
    return (optimalFamily, depParams, tau_hat)

if __name__=='__main__':

    from copularnd import copularnd
    from invcopulastat import invcopulastat
    from scipy.stats import norm
    from scipy.stats import expon

    # some tests on the copula multinomial signature
    tau = 0.4
    K = 4
    mnsig = copulamnsig('Gumbel',K,'kendall',tau)
    # iterate through mnsig to make sure we add upto 1 as a simple sanity check
    val_total = 0
    for ii in range(0,len(mnsig)):
        val_total = val_total + mnsig[ii]  #['val']
        
    if(np.isclose(val_total, 1.0)):
        print 'CopulaMNSig total probability check passed!'
    else:
        print 'CopulaMNSig total probability check failed!'
    
    M = 1000
    N = 2
    
    ###################### GAUSSIAN COPULA EXPERIMENT #######################
    # generate samples of the Gaussian copula with tau same as the
    # empirical signature we calculated above
    r = invcopulastat('Gaussian', 'kendall', tau)
    Rho = np.array([[1.0,r],[r,1.0]])
    U = copularnd('Gaussian', M, Rho)
    
    X1 = norm.ppf(U[:,0])       # assume mean=0, var=1
    X2 = expon.ppf(U[:,1])      # assume mean=0, var=1
    
    # combine X and Y into the joint distribution w/ the copula
    X = np.vstack((X1,X2))
    X = X.T
        
    ret = optimalCopulaFamily(X)
    print 'Input: Gaussian Copula. Output: ' + ret[0] + ' copula with tau_hat=' + str(ret[2])

    ###################### GUMBEL COPULA EXPERIMENT #######################
    # generate samples of the gumbel copula with tau same as the
    # empirical signature we calculated above
    alpha = invcopulastat('Gumbel', 'kendall', tau)
    U = copularnd('Gumbel', M, N, alpha)
    
    X1 = norm.ppf(U[:,0])       # assume mean=0, var=1
    X2 = expon.ppf(U[:,1])      # assume mean=0, var=1
    
    # combine X and Y into the joint distribution w/ the copula
    X = np.vstack((X1,X2))
    X = X.T
        
    ret = optimalCopulaFamily(X)
    print 'Input: Gumbel Copula. Output: ' + ret[0] + ' copula with tau_hat=' + str(ret[2])
    
    ###################### FRANK COPULA EXPERIMENT #######################
    alpha = invcopulastat('Frank', 'kendall', tau)
    U = copularnd('Frank', M, N, alpha)
    
    X1 = norm.ppf(U[:,0])       # assume mean=0, var=1
    X2 = expon.ppf(U[:,1])      # assume mean=0, var=1
    
    # combine X and Y into the joint distribution w/ the copula
    X = np.vstack((X1,X2))
    X = X.T
        
    ret = optimalCopulaFamily(X)
    print 'Input: Frank Copula. Output: ' + ret[0] + ' copula with tau_hat=' + str(ret[2])
    
    ###################### CLAYTON COPULA EXPERIMENT #######################
    alpha = invcopulastat('Clayton', 'kendall', tau)
    U = copularnd('Clayton', M, N, alpha)
    
    X1 = norm.ppf(U[:,0])       # assume mean=0, var=1
    X2 = expon.ppf(U[:,1])      # assume mean=0, var=1
    
    # combine X and Y into the joint distribution w/ the copula
    X = np.vstack((X1,X2))
    X = X.T
        
    ret = optimalCopulaFamily(X)
    print 'Input: Clayton Copula. Output: ' + ret[0] + ' copula with tau_hat=' + str(ret[2])
    
    # TODO: test for multivariate dataset, where # vars > 2
