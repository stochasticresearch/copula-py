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
    
    mnsig = []
    for coord in coords_list:
        # compute the C-volume and store
        u1v1 = coord[0]
        u1v2 = coord[1]
        u2v1 = coord[2]
        u2v2 = coord[3]
        try:
            val = cvolume(family, u1v1, u1v2, u2v1, u2v2, *args)
        except ValueError:
            val = np.array([-1])        # for compatibility we put the numpy wrapper

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
    # replace any 0 values w/ smallest possible float value
    empirical_mnsig[empirical_mnsig==0] = np.spacing(1)
    
    # compute the multinomial signature for each of the copula families specified
    # and simultaneously compute the kullback leibler divergence between the empirical
    # and the computed, and store that info
    distances = {}
    for family in family_search:
        # because the Clayton and Gumbel Copula's have restrictions for the valid values of
        # Kendall's tau, we do checks here to ensure those restrictions are met, because there
        # will be a certain variance associated with the tau_hat measurement
        
        if(family.lower()=='clayton'):
            # here we add some additional optimizatons as follows.  We know that the Clayton copula
            # captures only positive concordance.  Like any estimator, tau_hat will have some variance
            # associated with it.  Thus, the optimization we make is as follows, if tau_hat is within
            # a configurable amount less than 0, then we will set tau_hat to 0 and continue processing.  
            # However, if tau_hat is greater than that, we theoretically wouldn't have to test against 
            # the Clayton copula model, so we set the KL-divergence to be infinity to exclude 
            # this family from being selected
            if(tau_hat<-0.05):
                distances[family] = np.inf
                continue
            elif(tau_hat>=-0.05 and tau_hat<0):
                tau_hat = 0
            elif(tau_hat>=1):
                tau_hat = 1 - np.spacing(1)     # as close to 1 as possible in our precision
        elif(family.lower()=='gumbel'):
            # here we add some additional optimizatons as follows.  We know that the Gumbel copula
            # captures only positive concordance.  Like any estimator, tau_hat will have some variance
            # associated with it.  Thus, the optimization we make is as follows, if tau_hat is within
            # a configurable amount less than 0, then we will set tau_hat to 0 and continue processing.  
            # However, if tau_hat is greater than that, we theoretically wouldn't have to test against 
            # the Gumbel copula model, so we set the KL-divergence to be infinity to exclude 
            # this family from being selected
            if(tau_hat<-0.05):
                distances[family] = np.inf
                continue
            elif(tau_hat>=-0.05 and tau_hat<0):
                tau_hat = 0
            elif(tau_hat>=1):
                tau_hat = 1 - np.spacing(1)     # as close to 1 as possible in our precision
        # any other copula families with restrictions can go here
        
        mnsig = copulamnsig(family,K,'kendall',tau_hat)
        # replace any 0 values w/ smallest possible float value
        mnsig[mnsig==0] = np.spacing(1)
        
        # compute KL divergence, see
        # http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.stats.entropy.html
        distances[family] = entropy(mnsig, empirical_mnsig)
        
    # search for the minimum distance, that is the optimal copula family to use
    minDistance = np.inf
    for family, distance in distances.iteritems():
        if distance<minDistance:
            minDistance = distance
            optimalFamily = family
    
    depParams = invcopulastat(optimalFamily, 'kendall', tau_hat)
    
    return (optimalFamily, depParams, tau_hat)

def testHELM(tau, M, N, familyToTest, numMCSims, copulaFamiliesToTest):
    results = {}
    for fam in copulaFamiliesToTest:
        results[fam.lower()] = 0
    
    for ii in range(0,numMCSims):
        # generate samples of the requested copula with tau same as the
        # empirical signature we calculated above
        if(familyToTest.lower()=='gaussian'):
            r = invcopulastat(familyToTest, 'kendall', tau)
            
            Rho = np.empty((N,N))
            for jj in range(0,N):
                for kk in range(0,N):
                    if(jj==kk):
                        Rho[jj][kk] = 1
                    else:
                        Rho[jj][kk] = r
            try:
                U = copularnd(familyToTest, M, Rho)
            except ValueError:
                # copularnd will throw a ValueError if Rho is not a positive semidefinite matrix
                return results      # return 0, which will then be ignored by tests
                
        else:       # assume Clayton, Frank, or Gumbel
            try:
                alpha = invcopulastat(familyToTest, 'kendall', tau)
                U = copularnd(familyToTest, M, N, alpha)
            except ValueError:
                continue
            
        lst = []
        for jj in range(0,N):
            U_conditioned = U[:,jj]
            # if there are any 1's, condition it
            U_conditioned[U_conditioned==1] = 0.99
            if(jj%2==0):
                lst.append(norm.ppf(U_conditioned))
            else:
                lst.append(expon.ppf(U_conditioned))
        
        # combine X and Y into the joint distribution w/ the copula
        X = np.vstack(lst)
        X = X.T
                    
        ret = optimalCopulaFamily(X, family_search=copulaFamiliesToTest)
        ret_family = ret[0].lower()
        # aggregate results
        results[ret_family] = results[ret_family] + 1.0
        
        # display some progress
        sys.stdout.write("\rComputing " + str(familyToTest) + " Copula (DIM=%d) (tau=%f)-- %d%%" % (N,tau,ii+1))
        sys.stdout.flush()
    
    sys.stdout.write("\r")
    
    # convert results to percentage
    for fam in copulaFamiliesToTest:
        results[fam.lower()] = results[fam.lower()]/float(numMCSims) * 100
    
    return results

def plotPieChartResults(results, family, title):
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']      # for the pie chart
    # explode the Gaussian portion fo the pychart
    expTup = [0,0,0,0]
    expTup[results.keys().index(family.lower())] = 0.1
    plt.pie(results.values(), explode=expTup, labels=results.keys(),
            colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
    plt.title(title)
    plt.show()
    

def testHELM_parametric(K,M,N,tauVec,families):
    # some tests on the copula multinomial signature
    
    # Monte-Carlo style simulations to test each copula generation
    numMCSims = 1000
    
    resultsAggregate = {}
    
    for family in families:    
        famResults = {}
        for tau in tauVec:
            results = testHELM(tau, M, N, family, numMCSims, families)
            famResults[tau] = results
        resultsAggregate[family] = famResults
    
    return resultsAggregate

def visualizeMNSig():
    # some tests on the copula multinomial signature
    
    K = 4
    M = 1000
    N = 3
    tauVec = np.arange(-0.9,0.95,0.05)
    # the families to test against and pick optimal copula
    families = ['Gaussian', 'Clayton', 'Gumbel', 'Frank']
    
    helmAccuracyResults = testHELM_parametric(K,M,N,tauVec,families)
        
    resultsAggregate = {}
    
    for family in families:
        famResults = {}
        for tau in tauVec:
            mnsig = copulamnsig(family,K,'kendall',tau)
            famResults[tau] = mnsig
        resultsAggregate[family] = famResults

    # visualize the results
    for tau in tauVec:
        # we would also like to visualize this copula on the side, to try to 
        # understand what may be a better way todo model selection
        try:
            r = invcopulastat('Gaussian', 'kendall', tau)
        except ValueError:
            r = -1
        Rho = np.empty((N,N))
        for jj in range(0,N):
            for kk in range(0,N):
                if(jj==kk):
                    Rho[jj][kk] = 1
                else:
                    Rho[jj][kk] = r
        
        try:
            alpha_clayton = invcopulastat('Clayton', 'kendall', tau)
        except ValueError:
            alpha_clayton = -1
        
        try:
            alpha_gumbel  = invcopulastat('Gumbel', 'kendall', tau)
        except ValueError:
            alpha_gumbel = -1
            
        try:
            alpha_frank   = invcopulastat('Frank', 'kendall', tau)
        except ValueError:
            alpha_frank   = -1
        
        if(r!=-1):
            try:
                U_gauss   = copularnd('Gaussian', M, Rho)
            except ValueError:
                U_gauss   = np.zeros((M,N))
        if(alpha_clayton!=-1):
            try:
                U_clayton = copularnd('Clayton', M, N, alpha_clayton)
            except ValueError:
                U_clayton   = np.zeros((M,N))
        if(alpha_frank!=-1):
            try:
                U_frank   = copularnd('Frank', M, N, alpha_frank)
            except ValueError:
                U_frank   = np.zeros((M,N))
        if(alpha_gumbel!=-1):
            try:
                U_gumbel  = copularnd('Gumbel', M, N, alpha_gumbel)
            except ValueError:
                U_gumbel  = np.zeros((M,N))
        
        # get each family's MN signature and plot it
        plt.figure(figsize=(30,20))
        
        plt.subplot(231)
        if(np.sum(resultsAggregate['Gaussian'][tau])>0):
            plt.plot(np.arange(1,K*K+1), resultsAggregate['Gaussian'][tau], 'b.-', label='Gaussian Copula')
        if(np.sum(resultsAggregate['Clayton'][tau])>0):
            plt.plot(np.arange(1,K*K+1), resultsAggregate['Clayton'][tau], 'g.-', label='Clayton Copula')
        if(np.sum(resultsAggregate['Gumbel'][tau])>0):
            plt.plot(np.arange(1,K*K+1), resultsAggregate['Gumbel'][tau], 'r.-', label='Gumbel Copula')
        if(np.sum(resultsAggregate['Frank'][tau])>0):
            plt.plot(np.arange(1,K*K+1), resultsAggregate['Frank'][tau], 'k.-', label='Frank Copula')
        
        plt.title(r'Copula Multinomial Signature $\tau$=' + "{0:.2f}".format(tau) + ' K=' + str(K))
        plt.legend()
        plt.grid()
        
        plt.subplot(232)
        if(r!=-1):
            plt.scatter(U_gauss[:,0], U_gauss[:,1])
        plt.grid()
        plt.title(r'Gaussian Copula, $\rho$=' + "{0:.2f}".format(r) + r' $\tau$=' + "{0:.2f}".format(tau))
        
        plt.subplot(233)
        if(alpha_clayton!=-1):
            plt.scatter(U_clayton[:,0], U_clayton[:,1])
        plt.grid()
        plt.title(r'Clayton Copula, $\alpha$=' + "{0:.2f}".format(alpha_clayton) + r' $\tau$=' + "{0:.2f}".format(tau))
        
        plt.subplot(235)
        if(alpha_frank!=-1):
            plt.scatter(U_frank[:,0], U_frank[:,1])
        plt.grid()
        plt.title(r'Frank Copula, $\alpha$=' + "{0:.2f}".format(alpha_frank) + r' $\tau$=' + "{0:.2f}".format(tau))
        
        plt.subplot(236)
        if(alpha_gumbel!=-1):
            plt.scatter(U_gumbel[:,0], U_gumbel[:,1])
        plt.grid()
        plt.title(r'Gumbel Copula, $\alpha$=' + "{0:.2f}".format(alpha_gumbel) + r' $\tau$=' + "{0:.2f}".format(tau))
        
        plt.subplot(234)
        # index manually to ensure accuracy
        cla = np.array([helmAccuracyResults['Clayton'][tau]['clayton'],
                        helmAccuracyResults['Gaussian'][tau]['clayton'],
                        helmAccuracyResults['Gumbel'][tau]['clayton'],
                        helmAccuracyResults['Frank'][tau]['clayton']])
        gau = np.array([helmAccuracyResults['Clayton'][tau]['gaussian'],
                        helmAccuracyResults['Gaussian'][tau]['gaussian'],
                        helmAccuracyResults['Gumbel'][tau]['gaussian'],
                        helmAccuracyResults['Frank'][tau]['gaussian']])
        gum = np.array([helmAccuracyResults['Clayton'][tau]['gumbel'],
                        helmAccuracyResults['Gaussian'][tau]['gumbel'],
                        helmAccuracyResults['Gumbel'][tau]['gumbel'],
                        helmAccuracyResults['Frank'][tau]['gumbel']])
        fra = np.array([helmAccuracyResults['Clayton'][tau]['frank'],
                        helmAccuracyResults['Gaussian'][tau]['frank'],
                        helmAccuracyResults['Gumbel'][tau]['frank'],
                        helmAccuracyResults['Frank'][tau]['frank']])
        ind = np.arange(4)
        width = 0.2
        p1 = plt.bar(ind,cla,width,color='b')
        p2 = plt.bar(ind,gau,width,color='g',bottom=cla)
        p3 = plt.bar(ind,gum,width,color='k',bottom=cla+gau)
        p4 = plt.bar(ind,fra,width,color='r',bottom=cla+gau+gum)
        plt.xticks(ind+width/2.,('Clayton', 'Gaussian', 'Gumbel', 'Frank'))
        plt.legend( (p1[0], p2[0], p3[0], p4[0]), ('Clayton', 'Gaussian', 'Gumbel', 'Frank') )

        plt.grid()
        plt.savefig(os.path.join('figures/HELM_performance/', 
                     'HELM_DIM_' + str(N) + '_tau_' + "{0:.2f}".format(tau) + ' _K_' + str(K) + '.png'))
        
        plt.close()


if __name__=='__main__':
    from copularnd import copularnd
    from invcopulastat import invcopulastat
    from scipy.stats import norm
    from scipy.stats import expon
    import sys
    import matplotlib.pyplot as plt
    import os

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
    
    # Monte-Carlo style simulations to test each copula generation
    numMCSims = 100
    # the families to test against and pick optimal copula
    families = ['Gaussian', 'Clayton', 'Gumbel', 'Frank']
    
    """
    for family in families:
        title = 'Reference Bivariate ' + str(family) + ' Copula - HELM Identification Breakdown'
        results = testHELM(tau, M, N, family, numMCSims, families)
        plotPieChartResults(results, family, title)
    
    N = 3
    for family in families:
        title = 'Reference Bivariate ' + str(family) + ' Copula - HELM Identification Breakdown'
        results = testHELM(tau, M, N, family, numMCSims, families)
        plotPieChartResults(results, family, title)
    """
    #tauVec = np.arange(-0.9,0.95,0.05)
    #resultsAggregate = testHELM_parametric(K,M,N,tauVec)
    
    visualizeMNSig()
    
    
    