## Part of R package 'stabledist' (part of the Rmetrics project).

## The stabledist R package is free software; you can redistribute it and/or
## modify it under the terms of the GNU Library General Public
## License as published by the Free Software Foundation; either
## version 2 of the License, or (at your option) any later version.
##
## This R package is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Library General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/


################################################################################
# FUNCTIONS:		DESCRIPTION:
#  stableMode		 Computes the mode of the stable DF
################################################################################

##' Computes the mode of the alpha stable distribution
##' @title Mode of the stable distribution
##' @param alpha
##' @param beta
##' @param beta.max for numerical purposes, values of beta too close to 1,
##'  are set to beta.max
##' @param tol numerical tolerance used in optimize()
##' @return a number, the stable mode
##' @author Diethelm Wuertz and Martin Maechler
stableMode <- function(alpha, beta, beta.max = 1 - 1e-11,
		       tol = .Machine$double.eps^0.25)
{
    stopifnot(0 < alpha, alpha <= 2, length(alpha) == 1,
	      -1 <= beta, beta <= 1, length(beta) == 1,
              length(beta.max) == 1)
    # Notes:
    #	# Test for values close to beta = 1
    #	alpha <- seq(0, 2, by = 0.1)
    #	ans <- matrix(NA, nrow=length(alpha), ncol = 4)
    #	for (i in 1:seq_along(alpha)) {
    #	  ans[i,] <- c(
    #	    stableMode(alpha = alpha[i], beta = 0.99 ),
    #	    stableMode(alpha = alpha[i], beta = 0.99999 ),
    #	    stableMode(alpha = alpha[i], beta = 0.99999999 ),
    #	    stableMode(alpha = alpha[i], beta = 0.99999999999))
    #	}
    #	cbind(alpha, ans),
    #
    #	alpha	       0.99	  0.99999    0.99999999 0.99999999999
    #	0.0    0.000000e+00  0.000000e+00  0.000000e+00	 0.000000e+00
    #	0.2   -3.214142e-01 -3.246759e-01 -3.246787e-01 -3.246788e-01
    #	0.4   -6.105318e-01 -6.158562e-01 -6.158616e-01 -6.158616e-01
    #	0.6   -6.550106e-01 -6.594746e-01 -6.594790e-01 -6.594790e-01
    #	0.8   -5.558811e-01 -5.590032e-01 -5.590063e-01 -5.590063e-01
    #	1.0   -4.271033e-01 -4.293078e-01 -4.293099e-01 -4.293099e-01
    #	1.2   -3.074015e-01 -3.090820e-01 -3.090804e-01 -3.090804e-01
    #	1.4   -2.050956e-01 -2.063979e-01 -2.063951e-01 -2.063951e-01
    #	1.6   -1.199623e-01 -1.208875e-01 -1.208853e-01 -1.208853e-01
    #	1.8   -5.098617e-02 -5.145758e-02 -5.145639e-02 -5.145639e-02
    #	2.0   -7.487432e-05 -7.487432e-05 -7.487432e-05 -7.487432e-05

    if(alpha * beta == 0)
	return(0)
    ## else
    if(beta > beta.max) beta <- beta.max

    optimize(dstable, interval = c(-0.7, 0)*sign(beta),
	     alpha = alpha, beta = beta, pm = 0,
	     maximum = TRUE, tol = tol)$maximum
}
