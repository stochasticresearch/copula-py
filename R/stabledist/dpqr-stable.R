# Part of R package 'stabledist' (part of the Rmetrics project).

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
## FUNCTIONS:		 DESCRIPTION:
##  dstable		  Returns density for stable DF
##  pstable		  Returns probabilities for stable DF
##  qstable		  Returns quantiles for stable DF
##  rstable		  Returns random variates for stable DF
## UTILITY FUNCTION	 DESCRIPTION:
##  .integrate2  	  Integrates internal functions for *stable
################################################################################


##==============================================================================

### MM TODO:

## 0) All d, p, q, q -- have identical parts
##  a) 'parametrization' (pm) check
##  b) checking, alpha, beta,
##  c) subdivisions etc {all but rstable}
## ---	to do: "Fix" these in dstable(), then copy/paste to others

##==============================================================================

pi2 <- pi/2 # - we use it so often


##' @title omega() according to Lambert & Lindsey (1999), p.412
##' @param gamma [dpqr]stable()'s scale parameter, > 0         -- of length 1
##' @param alpha [dpqr]stable()'s "main" parameter, in [0, 2]  -- of length 1
##' @return omega(.) = tan(pi/2 alpha) if alpha != 1 ...
.om <- function(gamma,alpha) {
    if(alpha != round(alpha)) # non-integer usual case
	tan(pi2*alpha)# not tanpi2() !
    else if(alpha == 1)
	(2/pi)*log(gamma)
    else 0 # for alpha = 0 or = 2
}

##' @title C_alpha - the tail constant
##' @param alpha numeric vector of stable tail parameters, in [0,2]
##' @return
##' @author Martin Maechler
C.stable.tail <- function(alpha, log = FALSE) {
    stopifnot(0 <= alpha, alpha <= 2)
    r <- alpha
    i0 <- alpha == 0
    r[i0] <- if(log) -log(2) else 0.5
    al <- alpha[!i0]
    r[!i0] <-
        if(log) lgamma(al)-log(pi)+ log(sin(al*pi2))
        else gamma(al)/pi * sin(al*pi2)
    if(any(a2 <- alpha == 2)) r[a2] <- if(log) -Inf else 0
    r
}

##' @title tan(pi/2*x), for x in [-1,1] with correct limits
##'   i.e. tanpi2(-/+ 1) == -/+ Inf
##' @param x numeric vector
##' @return numeric vector of values tan(pi/2*x)
##' @author Martin Maechler
tanpi2 <- function(x) {
    r <- x
    if(any(i <- x & x == round(x)))# excluding 0
	r[i] <- (2 - (x[i] %% 4))*Inf
    io <- which(!i)
    r[io] <- tan(pi2* x[io])
    r
}

##' @title cos(pi/2*x), for x in [-1,1] with correct limits
##'   i.e. cospi2(+- 1) == 0
##' @param x numeric vector
##' @return numeric vector of values cos(pi/2*x)
##' @author Martin Maechler
cospi2 <- function(x) {
    r <- x
    if(any(i <- x == round(x)))
	r[i] <- as.numeric(x[i] == 0)# 1 or 0 - iff x \in [-1,1] !
    io <- which(!i)
    r[io] <- cos(pi2* x[io])
    r
}

##' According to Nolan's  "tail.pdf" paper, where he takes *derivatives*
##' of the tail approximation 1-F(x) ~ (1+b) C_a x^{-a}  to prove
##' that    f(x) ~  a(1+b) C_a x^{-(1+a)} ...
##'
##' @title tail approximation density for dstable()
##' @param x
##' @param alpha
##' @param beta
##' @param log if true, return  log(f(.))
##' @return
##' @author Martin Maechler
dPareto <- function(x, alpha, beta, log = FALSE) {
    if(any(neg <- x < 0)) { ## left tail
	x   [neg] <- -x	  [neg]
        beta <- rep(beta, length.out=length(x))
	beta[neg] <- -beta[neg]
    }
    if(log)
	log(alpha)+ log1p(beta)+ C.stable.tail(alpha, log=TRUE) -(1+alpha)*log(x)
    else
	alpha*(1+beta)* C.stable.tail(alpha)* x^(-(1+alpha))
}

pPareto <- function(x, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
    if(any(neg <- x < 0)) { ## left tail
	x   [neg] <- -x	  [neg]
        beta <- rep(beta, length.out=length(x))
	beta[neg] <- -beta[neg]
        stop("FIXME --- pPareto() is not correct for negative x")## switch  1-iF / iF
    }
    if(log.p) {
	if(lower.tail) ## log(1 - iF)
	    log1p(-(1+beta)* C.stable.tail(alpha)* x^(-alpha))
	else ## log(iF)
	    log1p(beta)+ C.stable.tail(alpha, log=TRUE) - alpha*log(x)
    } else {
	iF <- (1+beta)* C.stable.tail(alpha)* x^(-alpha)
	if(lower.tail) 1-iF else iF
    }
}

dstable <- function(x, alpha, beta,
		    gamma = 1, delta = 0, pm = 0, log = FALSE,
		    tol = 64*.Machine$double.eps, zeta.tol= NULL,
                    subdivisions = 1000)
{
    ## Original implemented by Diethelm Wuertz;
    ## Changes for efficiency and accuracy by Martin Maechler

    ## Description:
    ##	 Returns density for stable DF

    ## Details:
    ##	 The function uses the approach of J.P. Nolan for general
    ##	 stable distributions. Nolan derived expressions in form
    ##	 of integrals based on the charcteristic function for
    ##	 standardized stable random variables. These integrals
    ##	 can be numerically evaluated.

    ## Arguments:
    ##	 alpha = index of stability, in the range (0,2]
    ##	 beta  = skewness, in the range [-1, 1]
    ##	 gamma = scale, in the range (0, infinity)
    ##	 delta = location, in the range (-infinity, +infinity)
    ##	 param = type of parmeterization

    ## Note: S+ compatibility no longer considered (explicitly)

    ## Parameter Check:
    ## NB: (gamma, delta) can be *vector*s (vectorized along x)
    stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
	      -1 <= beta, beta	<= 1, length(beta) == 1,
	      0 <= gamma, length(pm) == 1, pm %in% 0:2,
	      tol > 0, subdivisions > 0)
    ## not an official argument {no doc!}:
    verbose <- getOption("dstable.debug", default=FALSE)
    ## Parameterizations:
    if (pm == 1) {
	delta <- delta + beta*gamma * .om(gamma,alpha)
    } else if (pm == 2) {
	delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
	gamma <- alpha^(-1/alpha) * gamma
    } ## else pm == 0

    ## Shift and Scale:
    x <- (x - delta) / gamma
    ans <-
	## Special Cases:
	if (alpha == 2) {
	    dnorm(x, mean = 0, sd = sqrt(2), log=log)
	} else if (alpha == 1 && beta == 0) {
	    dcauchy(x, log=log)
	} else {
	    ## General Case
	    if (alpha != 1) { ## 0 < alpha < 2	&  |beta| <= 1 from above
		tanpa2 <- tan(pi2*alpha)
                betan <- beta * tanpa2
		zeta <- -betan
		theta0 <- min(max(-pi2, atan(betan) / alpha), pi2)
		if(verbose) cat(sprintf(
		    "dstable(., alpha=%g, beta=%g,..): --> theta0=%g, zeta=%g,",
		    alpha, beta, theta0, zeta))
		if(is.null(zeta.tol)) {
		    zeta.tol <-
			if(betan == 0) .4e-15
			else if(1-abs(beta) < .01 || alpha < .01) 2e-15 else 5e-5
		    if(verbose) cat(sprintf(" --> zeta.tol= %g", zeta.tol))
		}
		else stopifnot(is.numeric(zeta.tol), zeta.tol >= 0)
		if(verbose) cat("\n")

		## Loop over all x values ( < , = , or >  zeta):
		vapply(x, .fct1, 0.,
		       zeta=zeta, alpha=alpha, beta=beta, theta0=theta0, log=log,
		       verbose=verbose,
		       tol=tol, zeta.tol=zeta.tol, subdivisions=subdivisions)
	    }
	    ## Special Case alpha == 1	and  -1 <= beta <= 1 (but not = 0) :
	    else { ## (alpha == 1)  and	 0 < |beta| <= 1  from above
		## Loop over all x values:
		vapply(x, function(z) {
		    if (z >= 0) {
			.fct2( z , beta, log=log, tol=tol, subdivisions=subdivisions)
		    } else {
			.fct2(-z, -beta, log=log, tol=tol, subdivisions=subdivisions)
		    }
		}, 0.)
	    }
	}

    i0 <- ans == (if(log)-Inf else 0) # --> we can do better using asymptotic:
    if(any(i0)) {
	d <- dPareto(x[i0], alpha, beta, log=log)
	## do recycle correctly:
	if(length(gamma) > 1)
	    gamma <- rep(gamma, length.out=length(x))[i0]
	ans[i0] <- if(log) d - log(gamma) else d/gamma
    }
    if(any(io <- !i0)) {
	d <- ans[io]
	if(length(gamma) > 1)
	    gamma <- rep(gamma, length.out=length(x))[io]
	ans[io] <- if (log) d - log(gamma) else d/gamma
    }
    ans
}## {dstable}

## ------------------------------------------------------------------------------

.large.exp.arg <- -(.Machine$double.min.exp * log(2)) ## == 708.396...
##' @title  x*exp(-x)  numerically stably, with correct limit 0 for x --> Inf
##' @param x  numeric
##' @return x*exp(x)
##' @author Martin Maechler
x.exp.m.x <- function(x) {
    r <- x*exp(-x)
    if(any(nax <- is.na(x)))
	r[nax] <- NA_real_
    if(any(lrg <- !nax & x > .large.exp.arg))# e.g. x == Inf
	r[lrg] <- 0
    r
}

.e.plus <- function(x, eps) x + eps* abs(x)
.e.minus<- function(x, eps) x - eps* abs(x)
pi2.. <- function(eps) pi2 * (1 - eps) ## == .e.minus(pi/2, eps), slight more efficiently

##' dstable() for very small alpha > 0
##' ok only for  x > zeta := - beta * tan(pi/2 *alpha)
dstable.smallA <- function(x, alpha, beta, log=FALSE) {
    r <- log(alpha) + log1p(beta) - (1 + log(2*x + pi*alpha*beta))
    if(log) r else exp(r)
}

## 1e-17: seems "good", but not "optimized" at all -- hidden for now
.alpha.small.dstable <- 1e-17

.fct1 <- function(x, zeta, alpha, beta, theta0, log,
                  tol, subdivisions, zeta.tol,
                  verbose = getOption("dstable.debug", default=FALSE))
{
    ## --- dstable(x, alpha, beta, ..)  for alpha < 2 ---

    ## For  x = zeta, have special case formula [Nolan(1997)];
    ## need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
    stopifnot(is.finite(zeta))
    x.m.zet <- abs(x - zeta)
    f.zeta <- function(log)
        if(log)
            lgamma(1+1/alpha)+ log(cos(theta0)) - (log(pi)+ log1p(zeta^2)/(2*alpha))
        else
            gamma(1+1/alpha)*cos(theta0) / (pi*(1+zeta^2)^(1/(2*alpha)))
    ## Modified: originally was	 if (z == zeta),
    ## then (D.W.)   if (x.m.zet < 2 * .Machine$double.eps)
    ## then (M.M.)   if (x.m.zet <= 1e-5 * abs(x))
    if(is.finite(x) && x.m.zet <= zeta.tol * (zeta.tol+ max(abs(x),abs(zeta)))) {
        if(verbose)
	    cat(sprintf(".fct1(%.11g, %.10g,..): x ~= zeta => using f.zeta()\n",
                        x, zeta))
	return(f.zeta(log))
    }
    ## the real check should be about the feasibility of g() below, or its integration

    smallAlpha <- (alpha < .alpha.small.dstable)
    if(x < zeta) {
	theta0 <- -theta0 # see Nolan(1997), Thm.1 (c)
	if(smallAlpha) {
	    beta <- -beta
	    x <- -x
	}
    }

    if(smallAlpha) {
        ## here, *MUST* have  __ x > zeta __
	if(verbose)
	    cat(sprintf(".fct1(%.11g, %.10g,..): small alpha=%g\n",
			x, zeta, alpha))
	return(dstable.smallA(x, alpha, beta, log=log))
    }

    ## constants ( independent of integrand g1(th) = g*exp(-g) ):
    ## zeta <- -beta * tan(pi*alpha/2)
    ## theta0 <- (1/alpha) * atan( beta * tan(pi*alpha/2))
    ## x.m.zet <- abs(x - zeta)
    ##-------->>>  identically as in .FCT1() for pstable() below: <<<-----------
    a_1 <- alpha - 1
    cat0 <- cos(at0 <- alpha*theta0)
    ##' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
    ##'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf
    g <- function(th) {
	r <- th
	## g(-pi/2) or g(pi/2) could become  NaN --> work around
	i.bnd <- abs(pi2 -sign(a_1)*th) < 64*.Machine$double.eps
	r[i.bnd] <- 0
	th <- th[io <- !i.bnd]
	att <- at0 + alpha*th ## = alpha*(theta0 + theta)
	r[io] <- (cat0 * cos(th) * (x.m.zet/sin(att))^alpha)^(1/a_1) * cos(att-th)
	r
    }
    ## Function to integrate: dstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1(u) du
    g1 <- function(th) {
	## g1 :=  g(.) exp(-g(.))
	x.exp.m.x( g(th) )
    }
    c2 <- ( alpha / (pi*abs(a_1)*x.m.zet) )

    ## Now, result = c2 * \int_{-t0}^{pi/2}  g1(u) du  ,  we "only" need the integral
    ## where however, g1(.) may look to be (almost) zero almost everywhere and just have a small peak
    ## ==> Find the peak, split the integral into two parts of for intervals  (t0, t_max) + (t_max, pi/2)

    ## However, this may still be bad, e.g., for dstable(71.61531, alpha=1.001, beta=0.6),
    ## or  dstable(1.205, 0.75, -0.5)
    ##   the 2nd integral was "completely wrong" (basically zero, instead of ..e-5)

    ## NB: g() is monotone, see above
    if((alpha >= 1 &&
	((!is.na(g. <- g( pi2	)) && g. > .large.exp.arg) || identical(g(-theta0), 0))) ||
       (alpha  < 1 &&
	((!is.na(g. <- g(-theta0)) && g. > .large.exp.arg) || identical(g(pi2), 0)))) {
	## g() is numerically too large *or* 0 even where it should be inf
	## ===>	 g() * exp(-g()) is 0 everywhere
	if(verbose)
	    cat(sprintf(".fct1(%.11g, %.10g,..): g() is 'Inf' (or 0) ==> result 0", x,zeta))
	return(if(log)-Inf else 0)
    }

    g. <- if(alpha >= 1) g(.e.plus(-theta0, 1e-6)) else g(pi2..(1e-6))
    if(is.na(g.))# g() is not usable --- FIXME rather use *asymptotic dPareto()?
	if(max(x.m.zet, x.m.zet / abs(x)) < .01)
	    return(f.zeta(log))

    if(verbose)
	cat(sprintf(".fct1(%.11g, %.10g,..): c2*sum(r[1:4])= %.11g*", x,zeta, c2))

    Int <- function(a,b)
	.integrate2(g1, lower = a, upper = b,
                    subdivisions=subdivisions, rel.tol= tol, abs.tol= tol)

    ## We know that the maximum of g1(.) is = exp(-1) = 0.3679  "at" g(.) == 1
    ## find that by uniroot :
    ## g(.) == 1  <==>  log(g(.)) == 0   --- the latter is better conditioned,
    ##                                       e.g., for (x = -1, alpha = 0.95, beta = 0.6)
    ## the former is better for  dstable(-122717558, alpha = 1.8, beta = 0.3, pm = 1)
    ## However, it can be that the maximum is at the boundary,  and
    ## g(.) > 1 everywhere or  g(.) < 1  everywhere  {in that case we could revert to optimize..}

    if((alpha >= 1 && !is.na(g. <- g(pi2)) && g. > 1) ||
       (alpha <	 1 && !is.na(g. <- g(pi2)) && g. < 1))
        g1.th2 <- g1( theta2 <- pi2..(1e-6) )
    else if((alpha <  1 && g(-theta0) > 1) ||
            (alpha >= 1 && g(-theta0) < 1))
        g1.th2 <- g1( theta2 <- .e.plus(-theta0, 1e-6) )
    else {
        ## when alpha ~=< 1 (0.998 e.g.),  g(x) is == 0 (numerically) on a wide range;
        ## uniroot is not good enough, and we should *increase* -theta0
        ## or decrease pi2 such that it can find the root:
        l.th <- -theta0
        u.th <- pi2
        if(alpha < 1) { ## g() is *in*creasing from 0 ..
	    while ((g.t <- g(.th <- (l.th + pi2)/2)) == 0) l.th <- .th
	    if(g.t == 1)# decrease upper limit {needed, e.g. for alpha = 1e-20}
                while ((g.t <- g(.th <- (l.th + u.th)/2)) == 1) u.th <- .th
	    if(abs(u.th - l.th) < 1e-13)# do not trust g()
                return(if(log)-Inf else 0)
	    if(verbose >= 2)
		cat(sprintf("\n -theta0=%g %s l.th=%g .. u.th=%g <= pi/2\n",
			    -theta0, if(-theta0 == l.th) "=" else "<",
			    l.th, u.th))
        }

        ur1 <- uniroot(function(th) g(th) - 1,
                       lower = l.th, upper = u.th, tol = .Machine$double.eps)
        ## consider using safeUroot() [ ~/R/Pkgs/copula/R/safeUroot.R ] !!
        ur2 <- tryCatch(uniroot(function(th) log(g(th)),
                                lower = l.th, upper = u.th, tol = .Machine$double.eps),
                        error=function(e)e)
	g.1 <- x.exp.m.x(ur1$f.root+1)
	g.2 <- if(inherits(ur2, "error")) -Inf else x.exp.m.x(exp(ur2$f.root))
        if(g.1 >= g.2) {
            theta2 <- ur1$root
            g1.th2 <- g.1 ## == g1(theta2)
        } else {
            theta2 <- ur2$root
            g1.th2 <- g.2
        }
    }
    ## now, because g1()'s peak (at th = theta2) may be extreme, we find two more intermediate values
    ## NB: Theoretically: Max = 0.3679 = g1(theta2)  ==> 1e-4 is a very small fraction of that
    ## to the left:
    eps <- 1e-4
    if((do1 <- g1.th2 > eps && g1(-theta0) < eps))
	th1 <- uniroot(function(th) g1(th) - eps, lower = -theta0, upper = theta2,
		       tol = tol)$root
    if((do4 <- g1.th2 > eps && g1(pi2) < eps))
	## to the right:
	th3 <- uniroot(function(th) g1(th) - eps, lower = theta2, upper = pi2,
		       tol = tol)$root

    if(do1) {
        r1 <- Int(-theta0, th1)
        r2 <- Int(         th1, theta2)
    } else {
        r1 <- 0
        r2 <- Int(-theta0,      theta2)
    }
    if(do4) {
        r3 <- Int(              theta2, th3)
        r4 <- Int(                      th3, pi2)
    } else {
        r3 <- Int(              theta2,      pi2)
        r4 <- 0
    }
    if(verbose)
	cat(sprintf("(%6.4g + %6.4g + %6.4g + %6.4g)= %g\n",
		    r1,r2,r3,r4, c2*(r1+r2+r3+r4)))
    if(log)
	log(c2)+ log(r1+r2+r3+r4)
    else
	c2*(r1+r2+r3+r4)
} ## {.fct1}


## ------------------------------------------------------------------------------

##' Auxiliary for dstable()  only used when alpha == 1 :
##' @param x  numeric *scalar*, >= 0
##' @param beta  0 < |beta| <= 1
##' @param tol
##' @param subdivisions
.fct2 <- function(x, beta, log, tol, subdivisions,
                  verbose = getOption("dstable.debug", default=FALSE))
{
    i2b <- 1/(2*beta)
    p2b <- pi*i2b # = pi/(2 beta)
    ea <- -p2b*x
    if(is.infinite(ea)) return(if(log)-Inf else 0)

    ##' g() is strictly monotone;
    ##' g(u) := original_g(u*pi/2)
    ##'  for beta > 0: increasing from g(-1) = 0   to  g(+1) = Inf
    ##'  for beta < 0: decreasing from g(-1) = Inf to  g(+1) = 0
    ##t0 <- -sign(beta)*pi2# g(t0) == 0  mathematically, but not always numerically
    u0 <- -sign(beta)# g(u0) == 0  mathematically, but not always numerically
    g <- function(u) {
        r <- u
        r[i <- abs(u-u0) < 1e-10] <- 0
        u <- u[!i]
        th <- u*pi2
	h <- p2b+ th # == g'/beta where g' := pi/2 + beta*th = pi/2* (1 + beta*u)
	r[!i] <- (h/p2b) * exp(ea + h*tanpi2(u)) / cospi2(u)
        r
    }

    ## Function to Integrate; u is a non-sorted vector!
    g2 <- function(u) {
	## g2 = g(.) exp(-g(.))
	x.exp.m.x( g(u) )
    }

    ## We know that the maximum of g2(.) is = exp(-1) = 0.3679  "at" g(.) == 1
    ## find that by uniroot :
    ur <- uniroot(function(u) g(u) - 1, lower = -1, upper = 1, tol = tol)
    u2 <- ur$root

    r1 <- .integrate2(g2, lower = -1, upper = u2,
		     subdivisions = subdivisions, rel.tol = tol, abs.tol = tol)
    r2 <- .integrate2(g2, lower = u2, upper = 1,
		     subdivisions = subdivisions, rel.tol = tol, abs.tol = tol)
    if(verbose) {
        cc <- pi2*abs(i2b)
	cat(sprintf(".fct2(%.11g, %.6g,..): c*sum(r1+r2)= %.11g*(%6.4g + %6.4g)= %g\n",
		    x,beta, cc, r1, r2, cc*(r1+r2)))
    }
    if(log)
	log(pi2) + log(abs(i2b)) + log(r1 + r2)
    else
	pi2*abs(i2b)*(r1 + r2)
}## {.fct2}

### ------------------------------------------------------------------------------


pstable <- function(q, alpha, beta, gamma = 1, delta = 0, pm = 0,
                    lower.tail = TRUE, log.p = FALSE,
		    tol = 64*.Machine$double.eps, subdivisions = 1000)
{
    ## A function implemented by Diethelm Wuertz

    ## Description:
    ##	 Returns probability for stable DF

    x <- q
    ## Parameter Check:
    ## NB: (gamma, delta) can be *vector*s (vectorized along x)
    stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
	      -1 <= beta, beta	<= 1, length(beta) == 1,
	      0 <= gamma, length(pm) == 1, pm %in% 0:2,
	      tol > 0, subdivisions > 0)

    ## Parameterizations:
    if (pm == 1) {
	delta <- delta + beta*gamma * .om(gamma,alpha)
    } else if (pm == 2) {
	delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
	gamma <- alpha^(-1/alpha) * gamma
    } ## else pm == 0

    ## Shift and Scale:
    x <- (x - delta) / gamma

    ## Return directly
    ## ------  first, special cases:
    if (alpha == 2) {
	pnorm(x, mean = 0, sd = sqrt(2), lower.tail=lower.tail, log.p=log.p)
    } else if (alpha == 1 && beta == 0) {
	pcauchy(x, lower.tail=lower.tail, log.p=log.p)
    } else {

        retValue <- function(F, useLower) { ## (vectorized in F)
            if(useLower) {
                if(log.p) log(F) else F
            } else { ## upper: 1 - F
                if(log.p) log1p(-F) else 1 - F
            }
        }
	## General Case
	if (alpha != 1) { ## 0 < alpha < 2	&  |beta| <= 1	from above
	    tanpa2 <- tan(pi2*alpha)
	    zeta <- -beta * tanpa2
	    theta0 <- min(max(-pi2, atan(-zeta) / alpha), pi2)
            if(finSupp <- (abs(beta) == 1 && alpha < 1)) {
                ## has *finite* support    [zeta, Inf)  if beta ==  1
                ##                         (-Inf, zeta] if beta == -1
            }
	    ## Loop over all x values:
	    vapply(x, function(z) {
		if(finSupp) {
                    if(beta == 1 && z <= zeta)
                        return(retValue(0., useLower=lower.tail))
                    else if(beta == -1 && z >= zeta)
                        return(retValue(1., useLower=lower.tail))
                    ## else .. one of the cases below
                }
		if(abs(z - zeta) < 2 * .Machine$double.eps) {
		    ## FIXME? same problem as dstable
		    r <- if(lower.tail) (1/2- theta0/pi) else 1/2+ theta0/pi
		    if(log.p) log(r) else r
		} else {
		    useLower <-
			((z > zeta && lower.tail) ||
			 (z < zeta && !lower.tail))
		    ## FIXME: for alpha > 1 -- the following computes F1 = 1 -c3*r(x)
		    ## and suffers from cancellation when 1-F1 is used below:
		    giveI <- !useLower && alpha > 1 # if TRUE, .FCT1() return 1-F
		    .F1 <- .FCT1(z, zeta, alpha=alpha, theta0=theta0,
				 giveI = giveI,
				 tol = tol, subdivisions = subdivisions)
		    if(giveI)
			if(log.p) log(.F1) else .F1
		    else retValue(.F1, useLower=useLower)
		}
	    }, 0.)
	}
	## Special Case alpha == 1	and  -1 <= beta <= 1 (but not = 0) :
	else { ## (alpha == 1)	and	 0 < |beta| <= 1  from above
	    useL <-
		if(beta >= 0)
		    lower.tail
		else {
		    beta <- -beta
		    x <- -x
		    !lower.tail
		}
	    if(giveI <- !useL && !log.p)
		useL <- TRUE
	    ## Loop over all x values:
	    retValue(vapply(x, function(z)
			    .FCT2(z, beta = beta, tol=tol, subdivisions=subdivisions,
				  giveI = giveI),
			    0.),
		     useLower = useL)
	}
    }
}## {pstable}

## ------------------------------------------------------------------------------

##' Auxiliary for pstable()  (for alpha != 1)
.FCT1 <- function(x, zeta, alpha, theta0, giveI, tol, subdivisions,
                  verbose = getOption("pstable.debug", default=FALSE))
{
    if(is.infinite(x))
	return(if(giveI) 0 else 1)
    stopifnot(is.finite(zeta))
    x.m.zet <- abs(x - zeta)
    ##-------->>>  identically as in .fct1() for dstable() above: <<<-----------
    ## FIXME: also provide "very small alpha" case, as in .fct1()
    if(x < zeta) theta0 <- -theta0

    a_1 <- alpha - 1
    cat0 <- cos(at0 <- alpha*theta0)

    g <- function(th) {
	r <- th
	## g(-pi/2) or g(pi/2) could become  NaN --> work around
	i.bnd <- abs(pi2 -sign(a_1)*th) < 64*.Machine$double.eps
	r[i.bnd] <- 0
	th <- th[io <- !i.bnd]
	att <- at0 + alpha*th ## = alpha*(theta0 + theta)
	r[io] <- (cat0 * cos(th) * (x.m.zet/sin(att))^alpha)^(1/a_1) * cos(att-th)
	r
    }

    if(verbose) cat(sprintf(".FCT1(%9g, %10g, th0=%.10g, %s..): ",
			    x,zeta, theta0, if(giveI)"giveI=TRUE," else ""))

    ## as g() is montone, the integrand  exp(-g(.)) is too ==> maximum is at the boundary
    ## however, integration can be inaccuracte when g(.) quickly jumps from Inf to 0
    ## _BUT_  empirically I find that good values l.th / u.th below are *INDEPENDENT* of x,
    l.th <- .e.plus(-theta0, 1e-6)
    if(alpha > 1 && g(l.th) == Inf) {
        ur <- uniroot(function(t) 1-2*(g(t)==Inf), lower=l.th, upper=pi2,
                      f.lower= -1, f.upper= 1, tol = 1e-8)
        l.th <- ur$root
        if(verbose) cat(sprintf(" g(-th0 +1e-6)=Inf: unirt(%d it) -> l.th=%.10g ",
                                ur$iter, l.th))
    }
    u.th <- .e.minus(pi2, 1e-6)
    if(alpha < 1 && g(u.th) == Inf) {
        ur <- uniroot(function(t) 1-2*(g(t)==Inf), lower=l.th, upper=u.th,
                      f.upper= -1, tol = 1e-8)
        u.th <- ur$root
        if(verbose) cat(sprintf(" g(pi/2 -1e-6)=Inf: unirt(%d it) -> u.th=%.10g ",
                                ur$iter, u.th))
    }
    r <- .integrate2(function(th) exp(-g(th)),
                     lower = l.th, upper = u.th, subdivisions = subdivisions,
                     rel.tol = tol, abs.tol = tol)

    if(verbose) cat(sprintf("--> Int r= %.11g\n", r))
    if(giveI) { ## { ==> alpha > 1 ==> c1 = 1; c3 = -1/pi}
	## return (1 - F) = 1 - (1 -1/pi * r) = r/pi :
	r/pi
    } else {
	c1 <- if(alpha < 1) 1/2 - theta0/pi else 1
	c3 <- sign(1-alpha)/pi
	## FIXME: for alpha > 1, F = 1 - |.|*r(x)
        ##    <==> cancellation iff we eventually want 1 - F() [-> 'lower.tail']
	c1 + c3* r
    }
} ## {.FCT1}

## ------------------------------------------------------------------------------

##' Auxiliary for pstable()  only used when alpha == 1 :
##' @param x numeric *scalar*
##' @param beta  >= 0 here
##' @param tol
##' @param subdivisions
.FCT2 <- function(x, beta, tol, subdivisions, giveI = FALSE,
                  verbose = getOption("pstable.debug", default=FALSE))
{
    i2b <- 1/(2*beta)
    p2b <- pi*i2b # = pi/(2 beta)
    ea <- -p2b*x
    if(is.infinite(ea))
	return(R.D.Lval(if(ea < 0) ## == -Inf  ==> g(.) == 0	==> G2(.) == 1
			1 else 0,  ## == +Inf  ==> g(.) == Inf	==> G2(.) == 0
			lower.tail= !giveI))

    ##' g() is strictly monotone;
    ##' g(u) := original_g(u*pi/2)
    ##'	 for beta > 0: increasing from g(-1) = 0   to  g(+1) = Inf
    ##'	 for beta < 0: decreasing from g(-1) = Inf to  g(+1) = 0
    ## original_g :
    ## g <- function(th) {
    ##     h <- p2b+ th # == g'/beta where g' := pi/2 + beta*th
    ##     (h/p2b) * exp(ea + h*tan(th)) / cos(th)
    ## }
    ##t0 <- -pi2# g(t0) == 0  mathematically, but not always numerically
    u0 <- -1 # g(u0) == 0  mathematically, but not always numerically
    g <- function(u) {
        r <- u
        r[i <- abs(u-u0) < 1e-10] <- 0
        u <- u[!i]
        th <- u*pi2
	h <- p2b+ th # == g'/beta where g' := pi/2 + beta*th = pi/2* (1 + beta*u)
	r[!i] <- (h/p2b) * exp(ea + h*tanpi2(u)) / cospi2(u)
        r
    }

    if(verbose)
	cat(sprintf(".FCT2(%.11g, %.6g, %s..): ",
		    x,beta, if(giveI) "giveI=TRUE," else ""))

    ## g(-u0) == +Inf {at other end}, mathematically ==> exp(-g(.)) == 0
    ## in the outer tails, the numerical integration can be inaccurate,
    ## because g(.) jumps from 0 to Inf,  but is 0 almost always
    ##   <==> g1(.) = exp(-g(.)) jumps from 1 to 0 and is 1 almost everywhere
    ##  ---> the integration "does not see the 0" and returns too large..
    u. <- 1
    if(g(uu <- .e.minus(u., 1e-6)) == Inf) {
        ur <- uniroot(function(t) 1-2*(g(t)==Inf), lower=-1, upper= uu,
                      f.lower= +1, f.upper= -1, tol = 1e-8)
        u. <- ur$root
        if(verbose) cat(sprintf(" g(%g)=Inf: unirt(%d it) -> u.=%.10g",
                                uu, ur$iter, u.))
    }

    ##' G2(.) = exp(-g(.)) is strictly monotone .. no need for 'theta2' !
    G2 <- if(giveI) function(u) expm1(-g(u)) else function(u) exp(-g(u))
    r <- .integrate2(G2, lower = -1, upper = u.,
                     subdivisions = subdivisions, rel.tol = tol, abs.tol = tol) / 2
    if(verbose) cat(sprintf("--> Int r= %.11g\n", r))
    if(giveI) -r else r
}## {.FCT2}

### ------------------------------------------------------------------------------

## -- utilities  (==^== Macros in R's  src/nmath/dpq.h ) :
R.D.Lval <- function(p, lower.tail) if(lower.tail) p else (1 - p) #   p
R.D.Cval <- function(p, lower.tail) if(lower.tail) (1 - p) else p # 1 - p
## R.D.qIv <- function(p, log.p)  if(log.p) exp(p) else p       # p  in qF(p,..)

##' == R.D.Lval(R.D.qIv(p))  "==="  p  in qF !
R.DT.qIv <- function(p, lower.tail, log.p) {
    if(log.p) if(lower.tail) exp(p) else - expm1(p)
    else R.D.Lval(p, lower.tail)
}

##' == R.D.Cval(R.D.qIv(p))  "===" (1 - p) in qF
R.DT.CIv <- function(p, lower.tail, log.p) {
    if(log.p) if(lower.tail) -expm1(p) else exp(p)
    else R.D.Cval(p, lower.tail)
}

qstable <- function(p, alpha, beta, gamma = 1, delta = 0, pm = 0,
                    lower.tail = TRUE, log.p = FALSE,
                    tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0,
                    integ.tol = 1e-7, subdivisions = 200)
{
    ## A function implemented by Diethelm Wuertz

    ## Description:
    ##	 Returns quantiles for stable DF

    ## Parameter Check:
    ## NB: (gamma, delta) can be *vector*s (vectorized along x)
    stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
	      -1 <= beta, beta	<= 1, length(beta) == 1,
	      0 <= gamma, length(pm) == 1, pm %in% 0:2,
	      tol > 0, subdivisions > 0)

    ## Parameterizations:
    if (pm == 1) {
	delta <- delta + beta*gamma * .om(gamma,alpha)
    } else if (pm == 2) {
	delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
	gamma <- alpha^(-1/alpha) * gamma
    } ## else pm == 0

    result <-
	## Special Cases:
	if (alpha == 2)
	    qnorm(p, mean = 0, sd = sqrt(2), lower.tail=lower.tail, log.p=log.p)
	else if (alpha == 1 && beta == 0)
	    qcauchy(p, lower.tail=lower.tail, log.p=log.p)
	else { ## -------------- 0 < alpha < 2 ---------------
            .froot <- function(x, p) {
                pstable(q = x, alpha=alpha, beta=beta, pm = 0,
                        lower.tail=lower.tail, log.p=log.p,
                        tol=integ.tol, subdivisions=subdivisions) - p
            }
            ## for approximate interval:
            .qN <- function(p) qnorm  (p, mean = 0, sd = sqrt(2),
                                     lower.tail=lower.tail, log.p=log.p)
            .qC <- function(p) qcauchy(p, lower.tail=lower.tail, log.p=log.p)

            ## Calculate:
            qst1 <- function(pp) {

		## 1) Find narrow interval  [xmin, xmax]  -----------------------
		##    NB: will deal with a too narrow interval later
		p0 <- R.DT.qIv(pp, lower.tail=lower.tail, log.p=log.p)
		left <- p0 < 0.5
		if (beta < 0) {
		    xmin <- -R.DT.CIv(pp, lower.tail=lower.tail, log.p=log.p)/p0
		    xmax <- if (left) .qN(pp) else .qC(pp)
		}
		else if (beta > 0 ) {
		    xmin <- if (left) .qC(pp) else .qN(pp)
		    xmax <- p0/R.DT.CIv(pp, lower.tail=lower.tail, log.p=log.p)
		}
		else { ## (beta == 0)
		    xmin <- if (left) .qC(pp) else .qN(pp)
		    xmax <- if (left) .qN(pp) else .qC(pp)
		}
		if(xmin >= xmax) { # fixup interval such that xmin < xmax
		    fdx <- if(xmin == xmax) .01*max(1e-7, abs(xmin)) else 1.01*(xmin-xmax)
		    xmin <- xmin - fdx
		    xmax <- xmax + fdx
		    stopifnot(xmin < xmax)
		}

                ## 2) root-finding  pstable(..) = p  inside the interval: -------
		dx <- 1
		repeat {
		    root <- .unirootNA(.froot, interval = c(xmin, xmax), p = pp,
				       extendInt = if(lower.tail) "upX" else "downX",
				       tol=tol, maxiter=maxiter, trace=trace)
		    if(!is.na(root))
			break
		    xmin <- xmin- dx
		    xmax <- xmax+ dx
		    if(xmin == -Inf && xmax == +Inf)
			stop("could not find an interval for x where pstable(x,*) - p changes sign")
		    dx <- dx * 2
		}
		root
	    }
	    vapply(p, qst1, 0.)
        }

    ## Result:
    result * gamma + delta
}

## ------------------------------------------------------------------------------


rstable <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 0, test=0)
{
    ## Description:
    ##	 Returns random variates for stable DF

    ## slightly amended along  copula::rstable1

    ## Parameter Check:
    ## NB: (gamma, delta) can be *vector*s (vectorized along x)
    stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
	      -1 <= beta, beta	<= 1, length(beta) == 1,
	      0 <= gamma, length(pm) == 1, pm %in% 0:2)

    ## Parameterizations:
    if (pm == 1) {
	delta <- delta + beta*gamma * .om(gamma,alpha)
    } else if (pm == 2) {
	delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
	gamma <- alpha^(-1/alpha) * gamma
    } ## else pm == 0

    ## Calculate uniform and exponential distributed random numbers:
    
    uni1 <- runif(n)
    uni2 <- runif(n)
    cau1 <- rcauchy(n)
    if(test==1) {
        uni1 <- 0.2
        uni2 <- 0.3
        cau1 <- 0.2
    }
    
    theta <- pi * (uni1-1/2)
    w <- -log(uni2)

    result <-
        ## If alpha is equal 1 then:
        if (alpha == 1 & beta == 0) {
            cau1
            ## Otherwise, if alpha is different from 1:
        } else {
            ## FIXME: learn from nacopula::rstable1R()
	    b.tan.pa <- beta*tan(pi2*alpha)
	    theta0 <- min(max(-pi2, atan(b.tan.pa) / alpha), pi2)
            c <- (1+b.tan.pa^2)^(1/(2*alpha))
	    a.tht <- alpha*(theta+theta0)
            r <- ( c*sin(a.tht)/
                  (cos(theta))^(1/alpha) ) *
                      (cos(theta-a.tht)/w)^((1-alpha)/alpha)
            ## Use Parametrization 0:
            r - b.tan.pa
        
        }

    ## Result:
    result * gamma + delta
}


## ------------------------------------------------------------------------------


##' Numerically Integrate -- basically the same as R's	integrate()
##' --------------------- main difference: no errors, but warnings
.integrate2 <- function(f, lower, upper, ..., subdivisions, rel.tol, abs.tol,
			stop.on.error = FALSE)
{
    ri <- integrate(f, lower, upper, ..., subdivisions=subdivisions,
		    rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)
    if((msg <- ri[["message"]]) != "OK")
	warning(msg) ## NB: "roundoff error ..." happens many times
    ri[["value"]]
}
