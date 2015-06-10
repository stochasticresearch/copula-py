## This is from 'fBasics',  but so small we will not import
## nor export, but just use it ...

## Now needs R >= 3.1.0 with its new argument(s) 'extendInt' etc
.unirootNA <-
    function(f, interval, ...,
             lower = min(interval), upper = max(interval),
             f.lower = f(lower, ...), f.upper = f(upper, ...),
             extendInt = c("no", "yes", "downX", "upX"),
             check.conv = FALSE,
             tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
{
    # Arguments:
    #   see 'uniroot'

    # Value:
    #   Returns the x value of f where the root is located. If
    #   no root exists,  NA will be returned instead. In that case,
    #   the function doesn't terminate with an error  as
    #   the standard function uniroot().

    # Example:
    #   .unirootNA(sin, c(1, 2)); .unirootNA(sin, c(-1, 1))

    # If there is no Root:
    if(is.na(f.lower) || is.na(f.upper) || f.lower * f.upper > 0)
        return(NA)
    ## else there is one :
    uniroot(f, interval = interval, ...,
            lower=lower, upper=upper, f.lower=f.lower, f.upper=f.upper,
	    extendInt=extendInt, check.conv=check.conv,
	    tol=tol, maxiter=maxiter, trace=trace)$root
}

## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_STABLEDIST_CHECK_EXTRA")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
