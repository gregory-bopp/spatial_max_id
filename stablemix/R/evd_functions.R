#' R wrapper for C function to calculate GEV cdf
#' @param q matrix of quantiles at which to evaluate GEV
#' @param loc matrix of GEV location parameters
#' @param scale matrix of GEV scale parameters
#' @param shape matrix of GEV shape parameters
#' @param lower.tail (logical) should lower tail be returned?
#' @param log.p (logical) should probabilites be log-transformed?
#' @return GEV cumulative probabilities for upper or lower tail
#' @export
pevdM <-
  function (q,
            loc,
            scale,
            shape,
            lower.tail = TRUE,
            log.p = FALSE)
  {
    if (any(scale <= 0))
      stop("pevd: invalid scale argument.  Must be > 0.")
    q <- (q - loc) / scale
    pevdC(q, loc, scale, shape, lower.tail, log.p)
  }

#' R wrapper for C function to calculate GEV quantile function
#' @param p matrix of probabilities for which to calculate quantiles
#' @param loc matrix of GEV location parameters
#' @param scale matrix of GEV scale parameters
#' @param shape matrix of GEV shape parameters
#' @param lower.tail (logical) should lower tail be returned?
#' @return GEV quantiles
#' @export
qevdM <- function (p, loc, scale, shape)
{
  if (any(scale <= 0))
    stop("qevd: invalid scale argument.  Must be > 0.")
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=
      1)
    stop("qevd: invalid p argument.  Must have 0 < p < 1.")
  qevdC(p, loc, scale, shape)
}
