#' Add prefix to a variable name (helper)
#' @description If variable name is missing, titled "(Intercept)"
#' @param x variable name
#' @param prefix string to add before variable name with underscore
#'
#' @return prefixed string
#' @export
make_vnames <- function(x, prefix) {
  x <- gsub("[(]Intercept[)]", "", x)
  x <- ifelse(x == "", prefix, paste0(prefix, "_", x))
  return(x)
}

#' Log likelihood for observation model
#' @param lz (n-rep by nloc) matrix of log-transformed stable-mixture
#' observations
#' @param y (n-rep by nloc) matrix of stable-mixture observations on GEV scale
#' @param lA (n by number of obs_coord) matrix of 1/alpha-norm scaled
#'         kernel functions at each observation location.
#' @param lK (nloc by nbasis) matrix of log-transformed sum constrained kernel
#' functions
#' @param alpha stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param theta stable theta parameter (exponential tilting parameter). This
#' parameter controls the tail dependence class. \code{theta} = 0 gives
#' asymptotic dependence, and \code{theta} > 0 gives asymptotic independence.
#' @param mu GEV location matrix
#' @param sigma GEV scale matrix
#' @param xi GEV shape matrix
#' @param tol alpha tolerance to bound away from zero and one
#' @export
llikM <-
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
        ) + lz
      ) +
        sum(extRemes::devd(
          y,
          loc = mu,
          scale = sigma,
          shape = xi,
          type = "GEV",
          log = TRUE
        )) -
        sum(dlstabmixM(lz, alpha, theta, lK, logp = TRUE))
      if (is.na(ll)) {
        return(-Inf)
      }
      else{
        return(ll)
      }
    }
  }


#' Full-conditional for alpha parameter
#' @inheritParams llikM
#' @param lB matrix of log-basis scaling factor matrix (e.g. from Hougaard
#' distribution) dimension: n basis x n independent replicates
#' @return log-full conditional for alpha
#' @export
lpst_alpha <- function(alpha,
                       lz,
                       y,
                       lA,
                       lB,
                       lK,
                       theta,
                       mu,
                       sigma,
                       xi,
                       tol = 0.02) {
  if (alpha < tol | alpha > 1 - tol) {
    return(-Inf)
  }
  llhexpstab(lB, alpha, theta) +
    llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi) +
    log(1 / (1 - 2 * tol))
}

#' Full-conditional for theta parameter
#' @inheritParams llikM
#' @param lB matrix of log-basis scaling factor matrix (e.g. from Hougaard
#' distribution) dimension: n basis x n independent replicates
#' @return log-full conditional for theta
#' @export
lpst_theta <-
  function(theta, lz, y, lB, lA, lK, alpha, mu, sigma, xi) {
    if (theta < 0) {
      return(-Inf)
    }
    llhexpstab(lB, alpha, theta) +
      dnorm(theta, sd = 10, log = T) +
      llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi)
  }

#' Log-full conditional for random scaling factors on log-scale
#' @inheritParams llikM
#' @param lB matrix of log-basis scaling factor matrix (e.g. from Hougaard
#' distribution) dimension: n basis x n independent replicates
#' @return Log-full conditional for random scaling factors on log-scale
#' @export
lpst_lB <- function(lB, lz, y, lA, lK, alpha,
                    theta, mu, sigma, xi) {
  llhexpstab(lB, alpha, theta) +
    llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi)
}

#' Log-full conditional for kernel bandwidth parameter in Gaussian density case
#' @inheritParams llikM
#' @param tau kernel bandwidth
#' @param tau_rng range to constrain bandwidth to
#' @return log-full conditional for kernel bandwidth parameter
#' @export
lpst_tau <- function(tau,
                     lz,
                     y,
                     lA,
                     lK,
                     alpha,
                     theta,
                     mu,
                     sigma,
                     xi,
                     tau_rng) {
  if ((tau <= 0) | (tau < tau_rng[1]) | (tau > tau_rng[2])) {
    return(-Inf)
  }
  else{
    return(
      llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi) +
        dnorm(tau, sd = tau_rng[2], log = TRUE)
    )
  }
}

#' Log-prior function for log-Gaussian process kernels
#' @param lK log-Gaussian process kernels which sum to 1 at each location after
#' exponentiation. Matrix dimensions: nloc x nbasis
#' @param s observation location matrix (nloc x 2)
#' @param gvar GP variance
#' @param gscl GP scale/range
#' @return log-prior for kernel functions.
#' @export
llKbasis <- function(lK, s, gvar, gscl) {
  G <-
    lK[, -1] - matrix(lK[, 1], ncol = ncol(lK) - 1, nrow = nrow(lK))
  model <- RandomFields::RMexp(var = gvar, scale = gscl)
  RandomFields::RFlikelihood(model, x = s, data = G)$loglikelihood
}
#' Log-full conditional for log-GP basis variance parameter
#' @inheritParams llKbasis
#' @export
lpst_gvar <- function(gvar, s, lK, gscl) {
  if (gvar > 0) {
    llKbasis(lK, s, gvar, gscl) +
      dnorm(gvar, mean = 0, sd = 10, log = T)
  }
  else{
    return(-Inf)
  }
}

#' Log-full conditional for log-GP basis variance parameter
#' @inheritParams llKbasis
#' @param rng_max standard deviation for prior
#' @export
lpst_gscl <- function(gscl, s, lK, gvar, rng_max) {
  if (gscl > 0) {
    llKbasis(lK, s,  gvar, gscl) +
      dnorm(gscl,
            mean = 0,
            sd = rng_max,
            log = T)
  }
  else{
    return(-Inf)
  }
}

#' Log-full conditional for log-Gaussian basis functions
#' @inheritParams llikM
#' @param s observation location matrix (nloc x 2)
#' @param gvar GP variance
#' @param gscl GP scale/range
#' @export
lpst_lK <- function(lK,
                    s,
                    lz,
                    y,
                    lA,
                    alpha,
                    theta,
                    mu,
                    sigma,
                    xi,
                    gvar,
                    gscl) {
  llKbasis(lK, s, gvar, gscl) +
    llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi)
}

#' GP Log-likelihood wrapper (RandomFields)
#' @description RandomFields wrapper to calculate log-likelihood of GP with
#' exponential covariance function
#' @param gp Matrix of Gaussian process realizations (see RandomFields for
#' format)
#' @param s nloc x 2 matrix of observation location coordinates
#' @param gvar GP variance
#' @param gscl GP scale/range
#' @return log-likelihood of GP observations
#' @export
llgp <- function(gp, s, gvar, gscl) {
  model <- RandomFields::RMexp(var = gvar, scale = gscl)
  RandomFields::RFlikelihood(model, x = s, data = gp)$loglikelihood
}

#' Log-full conditional for marginal GEV parameters with GP priors
#' @inheritParams llikM
#' @param gp GEV parameter Gaussian process matrix
#' @param mpar_gvar variance parameter for marginal GEV GP
#' @param mpar_gscl scale/range parameter for marginal GEV GP
#' @param s nloc x 2 matrix of observation location coordinates
#' @param rng_max optionally define marginal parameter priors based on range of
#' observation coordinates.
#' @export
lpst_mpar_gp <- function(gp,
                         mpar_gvar,
                         mpar_gscl,
                         s,
                         lz,
                         y,
                         lA,
                         lK,
                         alpha,
                         theta,
                         mu,
                         sigma,
                         xi,
                         rng_max = NULL) {
  if ((mpar_gscl > 0) & (mpar_gvar > 0) & all(sigma > 0)) {
    if (is.null(rng_max)) {
      llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi) +
        llgp(gp, s, mpar_gvar, mpar_gscl) +
        dnorm(mpar_gvar, sd = 1, log = TRUE)
    }
    else{
      llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi) +
        llgp(gp, s, mpar_gvar, mpar_gscl) +
        dnorm(mpar_gvar, sd = 10, log = TRUE) +
        dnorm(mpar_gscl, sd = rng_max, log = TRUE)
    }
  }
  else{
    -Inf
  }
}

#' Log-full conditional for GP exp. covariance parameters
#' @param gp matrix of Gaussian process realizations (See RandomFields for form)
#' @param s nloc x 2 matrix of observation location coordinates
#' @param gvar GP variance
#' @param gscl GP scale/range
#' @param rng_max optionally define prior for range parameter in terms of the
#' range of the observation coordinates.
#' @export
lpst_gp_cov_par <- function(gp, s, gvar, gscl, rng_max = NULL) {
  if (is.null(rng_max)) {
    model <- RandomFields::RMexp(var = gvar, scale = gscl)
    RandomFields::RFlikelihood(model, x = s, data = gp)$loglikelihood +
      dnorm(gvar, sd = 10, log = TRUE)
  }
  else{
    model <- RandomFields::RMexp(var = gvar, scale = gscl)
    RandomFields::RFlikelihood(model, x = s, data = gp)$loglikelihood +
      dnorm(gvar, sd = 10, log = TRUE) +
      dnorm(gscl, sd = rng_max, log = TRUE)
  }
}

#' Log-full conditional for marginal GEV parameter GP mean function parameters
#' @param beta GP mean function coefficients in X beta
#' @inheritParams llikM
#' @param is.xi (logical) is modeled GP for the shape parameter?
#' @export
lpst_mpar_mean <-
  function(beta,
           lz,
           y,
           lA,
           lK,
           alpha,
           theta,
           mu,
           sigma,
           xi,
           is.xi = FALSE) {
    if (all(sigma > 0)) {
      if (!is.xi) {
        llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi) +
          sum(dnorm(beta, sd = 10, log = TRUE))
      }
      else{
        llikM(lz, y, lA, lK, alpha, theta, mu, sigma, xi) +
          sum(dnorm(beta, sd = 10, log = TRUE))
      }
    }
    else{
      -Inf
    }
  }


# Assuming the smoothness parameter of whittle-matern is nu =1/2
# Takes scale (gvar) and range (gscl) and returns micro ergodic par
#' @export
to_microe <- function(gvar, gscl) {
  return(gvar / gscl)
}

# Takes micro ergodic (microe) and range (gscl) and returns scale (gvar)
#' @export
from_microe <- function(microe, gscl) {
  return(microe * gscl)
}


#' Metropolis update for random scaling factors (supplied on log-scale)
#' @param lB_cur Current mcmc samples of log-scaling factor for a single year
#' @param lzsub log-stable mixture observations for the corresponding year to
#' the supplied lB_cur vector
#' @param ysub GEV scale observations after marginal transformations of lzsub
#' @param lK_cur (nloc by nbasis) matrix of log-kernel functions for current
#' mcmc iteration
#' @param alpha_cur stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param theta_cur stable theta parameter (exponential tilting parameter). This
#' parameter controls the tail dependence class. \code{theta} = 0 gives
#' asymptotic dependence, and \code{theta} > 0 gives asymptotic independence.
#' @param mu_cur GEV location matrix for current mcmc iteration
#' @param sigma_cur GEV scale matrix for current mcmc iteration
#' @param xi_cur GEV shape matrix for current mcmc iteration
#' @param pos index in adaptive variance updates matrix
#' @param next_pos next index in adaptive variance updates matrix
#' @param acpt_lB acceptance matrix of (1: accept) and (0: reject) for lB
#' scaling factors. Used for determining adaptive proposal updates.
#' @param lB_prop_var proposal variance vector for the lB scaling factors

#' @export
update_lB <-
  function(lB_cur,
           lzsub,
           ysub,
           lK_cur,
           alpha_cur,
           theta_cur,
           mu_cur,
           sigma_cur,
           xi_cur,
           pos,
           next_pos,
           acpt_lB,
           lB_prop_var) {
    lA_cur <-  mklA(lK_cur, lB_cur, alpha_cur)
    lB_prop <- lB_cur
    for (j in 1:length(lB_cur)) {
      lB_prop[j] <- rnorm_prop(lB_cur[j], tune_var = lB_prop_var[pos, j])
      lA_prop <- mklA(lK_cur, lB_prop, alpha_cur)
      if (log(runif(1)) <
          lpst_lB(
            lB_prop[j],
            lzsub,
            ysub,
            lA_prop,
            lK_cur,
            alpha_cur,
            theta_cur,
            mu_cur,
            sigma_cur,
            xi_cur
          ) -
          lpst_lB(
            lB_cur[j],
            lzsub,
            ysub,
            lA_cur,
            lK_cur,
            alpha_cur,
            theta_cur,
            mu_cur,
            sigma_cur,
            xi_cur
          ))
      {
        lB_cur[j] <- lB_prop[j]
        lA_cur <- lA_prop
        acpt_lB[pos, j] <- 1
      }
      else
      {
        lB_prop[j] <- lB_cur[j]
        acpt_lB[pos, j] <- 0
      }
    }
    return(list(lB_cur = lB_cur, acpt_lB = acpt_lB))
  }

#' Random walk Normal proposal
#' @param cur current mcmc sample to center proposal at
#' @param tune_var tuning variance for mcmc proposal
#' @return random walk Normal proposal
#' @export
rnorm_prop <- function(cur, tune_var) {
  return(rnorm(1, mean = cur, sd = sqrt(tune_var)))
}

#' Spring function for automatic proposal variance tuning update
#' @param cur_var current proposal variance
#' @param acpt_rt current acceptance rate (between 0 and 1)
#' @param opt_rt target acceptance rate (between 0 and 1)
#' @param gamma1 spring tuning parameter
#' @return new proposal variance
#' @export
update_var <- function(cur_var, acpt_rt, opt_rt, gamma1) {
  exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt))
}

#' Block Cholesky de-noising proposals for lK (see paper appendix for details)
#' @param lK_sub_cur Current values of log-Gaussian process basis functions at
#' a subset of the observation locations
#' @param Csub Covariance matrix of the Gaussian process at subset of obs
#' locations in \code{lK_sub_cur}
#' @param lK_prop_var_sub proposal variance for random walk
#' @param nbasis_per_rep (not currently used)
#' @return lK proposal on region subset
#' @export
chol_prop_lK <-
  function(lK_sub_cur,
           Csub,
           lK_prop_var_sub,
           nbasis_per_rep) {
    Lt <- chol(Csub)
    white <- t(solve(Lt)) %*% lK_sub_cur
    white_shifted <-
      white + rnorm(prod(dim(lK_sub_cur)), mean = 0, sd = lK_prop_var_sub)
    lK_prop <- t(Lt) %*% white_shifted
    lK_prop <- t(apply(lK_prop, 1, function(x) {
      x - log(sum(exp(x)))
    }))
    return(lK_prop)
  }
