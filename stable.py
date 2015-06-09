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
            # in the future we can add it in, or someone else can :)
            raise NotImplementedError('Generating RVs for Parametrization of 2 not yet implemented!')
        # else pm is 0
        
        ## Calculate uniform and exponential distributed random numbers:
        theta = math.pi * (uniform.rvs(size=n)-1.0/2.0)
        w = -1*np.log(uniform.rvs(size=n))
        
        if(alpha==1 and beta==0):
            retVal = cauchy.rvs(size=n)
        else:
            b_tan_pa = self.beta * np.tan(math.pi/2.0 * self.alpha)
            theta0 = np.min(np.array([np.max(np.array([-math.pi/2.0, np.arctan(b_tan_pa)/self.alpha])), math.pi/2.0]))
            c = np.power((1 + b_tan_pa * b_tan_pa ), 1.0/(2.0*self.alpha))
            a_tht = self.alpha * (theta+theta0)
            r = ( c*np.sin(a_tht)/np.power(np.cos(theta),(1/alpha)) ) * \
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
        return np.tan(math.pi/2.0 * self.alpha)
    elif(alpha==1):
        return (2.0/math.pi) * np.log(gamma)
    else:
        return 0

if __name__=='__main__':
    alpha = 0.4
    beta = 0.8
    gamma = 1
    delta = 0
    pm = 1
    s = stable(alpha, beta, gamma, delta, pm)
    print s.rvs()
    
    # TODO add test mode where we input the same uniform random variables and cauchy RV and see
    # if we get the same value as matlab?