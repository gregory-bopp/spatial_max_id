#' Convert GEV observations to log-stable mixture scale, but leave NA alone
#' @param y GEV scale observation matrix (each column corresponds to a different
#' location)
#' @param mu GEV location parameter matrix
#' @param sigma GEV scale parameter matrix
#' @param xi GEV shape parameter matrix
#' @param qfuns Quantile functions for log-stable mixture distribution. Should
#' be equal in length to the number of columns of \code{y} matrix.
#' @export
#' @return Marginally transformed matrix of y observations (to log-stable
#' mixture scale).
gev2lsm_w_na <- function (y, mu, sigma, xi, qfuns)
{
  if (max(length(mu), length(sigma), length(xi)) == 1) {
    ps <- extRemes::pevd(y, mu, sigma, xi)
  }
  else {
    ps <- pevdM(y, mu, sigma, xi)
  }
  if (is.matrix(ps)) {
    for (i in 1:ncol(ps)) {
      keep <- !is.na(ps[, i])
      ps[keep, i] <- qfuns[[i]](ps[keep, i])
    }
  }
  else {
    ps <- qfuns(ps)
  }
  return(ps)
}


#' Calculate observation model log-likelihood excluding NAs
#' @param lz matrix nrep x nloc. log-stable mixture observations
#' @param y matrix nrep x nloc. log-stable mixture observations transformed to
#' GEV scale with location, scale, and shape given by mu, sigma, and xi
#' @param lA (n by number of obs_coord) matrix of 1/alpha-norm scaled
#'         kernel functions at each observation location. Transformed to log
#'         scale.
#' @param lK (number of obs_coord by number of knot_coord) matrix of
#'         log-kernel functions at each observation location.
#' @param alpha Positive stable index
#' @param theta Hougaard exponential tilting parameter
#' @param mu GEV location parameter matrix
#' @param sigma GEV scale parameter matrix
#' @param xi GEV shape parameter matrix
#' @param tol tolerance bounding alpha away from zero and one
#' @export
#' @return observation model log-likelihood
llikM_w_na <-
  function(lz,
           y,
           lA,
           lK,
           alpha,
           theta,
           mu,
           sigma,
           xi,
           tol = 0.02) {
    if (alpha < tol | alpha > 1 - tol | theta < 0 | any(sigma <= 0)) {
      return(-Inf)
    }
    else{
      ll <- sum(
        extRemes::devd(
          exp(lz),
          loc = exp(lA * alpha),
          scale = alpha * exp(lA * alpha),
          shape = alpha,
          type = "GEV",
          log = TRUE
        ) + lz,
        na.rm = T
      ) +
        sum(extRemes::devd(
          y,
          loc = mu,
          scale = sigma,
          shape = xi,
          type = "GEV",
          log = TRUE
        ),
        na.rm = T) -
        sum(dlstabmixM(lz, alpha, theta, lK, logp = TRUE), na.rm = T)
      if (is.na(ll)) {
        return(-Inf)
      }
      else{
        return(ll)
      }
    }
  }

#' Calculate log-scores
#' @param yhold GEV scale data on holdout set
#' @param alpha_smps posterior samples of alpha (positive stable index)
#' @param theta_smps posterior samples of theta (Exponentail tilting par.)
#' @param lB_smps posterior samples of log-transformed basis scaling factors
#' (each row is a different MCMC sample, and the columns contain the C style
#' column-major ordering when transforming a matrix to a vector of lB as defined
#' in for example rstabmix.)
#' @param lK_smps posterior samples of log-transformed basis functions at
#' holdout locations. See \code{mcmc_lnorm_basis} for dimensions of array.
#' @param mu_smps posterior samples of GEV location parameters
#' holdout locations. Rows: mcmc samples. Columns: holdout locations
#' @param sigma_smps posterior samples of GEV scale parameters
#' holdout locations. Rows: mcmc samples. Columns: holdout locations
#' @param xi_smps posterior samples of GEV shape parameters
#' holdout locations. Rows: mcmc samples. Columns: holdout locations
#' @export
#' @return log-scores for holdout set.
calc_log_scores <-
  function(yhold,
           alpha_smps,
           theta_smps,
           lB_smps,
           lK_smps,
           mu_smps,
           sigma_smps,
           xi_smps) {
    nloc <- ncol(yhold)
    nyr <- nrow(yhold)
    L <- ncol(lB_smps) / nyr
    nsmps <- length(alpha_smps)

    lscores <- rep(NA, length = nsmps)
    for (i in 1:nsmps) {
      mu <- matrix(mu_smps[i, ],
                   nrow = nyr,
                   ncol = nloc,
                   byrow = T)
      sigma <-
        matrix(sigma_smps[i, ],
               nrow = nyr,
               ncol = nloc,
               byrow = T)
      xi <- matrix(xi_smps[i, ],
                   nrow = nyr,
                   ncol = nloc,
                   byrow = T)
      lB <- matrix(lB_smps[i, ], nrow = L, ncol = nyr)
      qfuns <-
        create_qlstabmix(alpha_smps[i], theta_smps[i], lK_smps[i, , ])
      lz <- gev2lsm_w_na(yhold, mu, sigma, xi, qfuns)
      lA <- mklA(lK_smps[i, , ], lB, alpha_smps[i])
      lscores[i] <-
        llikM_w_na(lz,
                   yhold,
                   lA,
                   lK_smps[i, , ],
                   alpha_smps[i],
                   theta_smps[i],
                   mu,
                   sigma,
                   xi)
    }
    return(lscores)
  }
