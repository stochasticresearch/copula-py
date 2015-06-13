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
from scipy.stats import cauchy

class stable:
    """
    Defines the stable distribution.  See 
    http://academic2.american.edu/~jpnolan/stable/chap1.pdf
    for more details.  This code was ported from the open source R implementation
    at http://cran.r-project.org/web/packages/stabledist/index.html
    """
    def __init__(self, alpha, beta, gamma, delta, pm):
        # do some parameter checking to be within stable distribution limits
        if(alpha <= 0 or alpha > 2 or
           beta  <-1  or beta  > 1 or
           gamma < 0  or 
           pm    < 0  or pm > 2):
            raise ValueError('Invalid specification of Stable Distribution!')
        
        if(pm==2):
            # NOT YET IMPLEMENTED ... I'm lazy and haven't implemented the DSTABLE 
            # function yet.  I don't need it for my purposes right now, but maybe
            # in the future we can add it in, or someone else can :)
            raise NotImplementedError('Generating RVs for Parametrization of 2 not yet implemented!')
        
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.pm = pm
        
        self.testMode = False
        
    def enableTestMode(self):
        self.testMode = True
        
    def rvs(self, size=1):
        """
        Generates size number of random variates of the stable distribution
        """
        n = size        # convenience var
        
        delta_old = self.delta
        gamma_old = self.gamma
        
        if(self.pm==1):
            self.delta = self.delta + self.beta*self.gamma * _om(self.gamma,self.alpha)
        elif(self.pm==2):
            # NOT YET IMPLEMENTED ... I'm lazy and haven't implemented the DSTABLE 
            # function yet.  I don't need it for my purposes right now, but maybe
            # in the future I can add it in, or someone else can :)
            raise NotImplementedError('Generating RVs for Parametrization of 2 not yet implemented!')
        # else pm is 0
        
        ## Calculate uniform and exponential distributed random numbers:
        uni1 = uniform.rvs(size=n)
        uni2 = uniform.rvs(size=n)
        cau1 = cauchy.rvs(size=n)
        if(self.testMode):
            uni1 = 0.2
            uni2 = 0.3
            cau1 = 0.2
        
        theta = math.pi * (uni1-1.0/2.0)
        w = -1*np.log(uni2)
        
        if(self.alpha==1 and beta==0):
            retVal = cau1
        else:
            b_tan_pa = self.beta * np.tan(math.pi/2.0 * self.alpha)
            theta0 = np.min(np.array([np.max(np.array([-math.pi/2.0, np.arctan(b_tan_pa)/self.alpha])), math.pi/2.0]))
            c = np.power((1 + b_tan_pa * b_tan_pa ), 1.0/(2.0*self.alpha))
            a_tht = self.alpha * (theta+theta0)
            r = ( c*np.sin(a_tht)/np.power(np.cos(theta),(1.0/self.alpha)) ) * \
                np.power( (np.cos(theta-a_tht)/w), ((1.0-self.alpha)/self.alpha) )
            
            retVal = r - b_tan_pa
            
        retVal = retVal * self.gamma + self.delta
        
        # reset class var's
        self.delta = delta_old
        self.gamma = gamma_old
        
        return retVal
    
##' @title omega() according to Lambert & Lindsey (1999), p.412
##' @param gamma [dpqr]stable()'s scale parameter, > 0         -- of length 1
##' @param alpha [dpqr]stable()'s "main" parameter, in [0, 2]  -- of length 1
##' @return omega(.) = tan(pi/2 alpha) if alpha != 1 ...
def _om(alpha, gamma):
    if(alpha != round(alpha)):
        return np.tan(math.pi/2.0 * alpha)
    elif(alpha==1):
        return (2.0/math.pi) * np.log(gamma)
    else:
        return 0

if __name__=='__main__':
    # TODO: we should automatically run the R script to generate the test values
    # that way, when we add more tests, it just magically works :)

    # TODO: should test this implementation of stable with rStable1 which is
    # specifically in the copula package

    testFilename = 'R/stabledist/stabledist_test.txt'
    testIdx = 1
    with open(testFilename) as f:
        for line in f:
            vals = line.split(',')
            alpha = float(vals[0].rstrip())
            beta  = float(vals[1].rstrip())
            n     = float(vals[2].rstrip())
            gamma = float(vals[3].rstrip())
            delta = float(vals[4].rstrip())
            pm    = float(vals[5].rstrip())
            R_output = float(vals[6].rstrip())
            
            s = stable(alpha, beta, gamma, delta, pm)
            s.enableTestMode()
            
            pyOutput = s.rvs()
            # compare against R calculations
            if(np.isclose(R_output, pyOutput)):
                print 'Test ' + str(testIdx) + ' Passed!'
            else:
                print 'Test ' + str(testIdx) + ' Failed!'
            testIdx = testIdx + 1
