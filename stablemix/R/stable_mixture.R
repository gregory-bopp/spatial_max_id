#' Calcuate pairwise Euclidean distances between two matrices with d-dim columns
#' @description Calculate pairwise Euclidean distances between all rows of
#' a matrix A and all rows of another matrix B. The output matrix will be of dimension
#' @param A Matrix 1 (nrow = m)
#' @param B Matrix 2 (nrow = n)
#' @return An m by n matrix of distances between the rows of A and rows of B,
#'  where m = \code{nrow(x)} and n = \code{nrow(y)}
#' @examples
#' A1 <- matrix(1:6, ncol = 2)
#' B1 <- matrix(1:8, ncol = 2)
#' calc_dist(A1, B1)
#'
#' @export
calc_dist <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  result = matrix(ncol = nrow(B), nrow = nrow(A))
  for (i in 1:nrow(A))
    for (j in 1:nrow(B))
      result[i, j] = sqrt(sum((A[i,] - B[j,]) ^ 2))

  result
}

#' Construct kernel basis matrix
#' @description Generate a set of kernel basis functions at the observation locations
#' in the matrix obs_coord such that the kernels sum to one at each observation
#' location. Kernels can either be two-dimensional Gaussian density functions
#' centered at knot_coord, or a realization of a log-Gaussian process. The
#' output is returned on the log-scale.
#' @param obs_coord  matrix of observation coordinates. The rows correspond
#'   to different locations, and the columns correspond to the dimensions of
#'   those coordinates.
#' @param nbasis number of kernel (basis) functions to generate. Kernels will
#' sum to one across all nbasis at each location.
#' @param type (string) type of kernel function to create. Should be one of
#' either c("smith", "br"). For fixed, Gaussian density
#' kernels, use "smith". For a log-Gaussian process kernels, use "br",
#' which stands for Brown-Resnick.
#' @param knot_coord matrix of knot locations. The rows correspond to
#'   different locations. Needed only if \code{type} = 'smith'.
#' @param kern_bw Kernel bandwidth (sd) for 'smith' kernel bases
#' @param gvar Log-Gaussian process variance term ('brown-resnick'). Needed
#' only if \code{type} = 'br'
#' @param gscl Log-Gaussian process scale/range term ('brown-resnick'). Needed
#' only if \code{type} = 'br'
#' @return Matrix of log-kernel functions. Each row corresponds to a different
#' \code{obs_coord} and each column contains a single basis function, so that
#' the dimension is nrow(obs_coord) x nbasis.
#' @export
#' @examples
#' # Not run
#' # Fixed, Gaussian density kernels
#' L <- 10        # Number of basis functions
#' n <- 100       # Number of spatial locations
#' v <- seq(0,1, l = L)
#' s <- seq(0.1, 0.9, l = n)
#' K <- make_kernel(s, L, "smith", knot_coord = v, kern_bw = 0.1)
#' matplot(s, exp(K), type = "l")
#'
#' # log-Gaussian process kernels
#' set.seed(2017)
#' L <- 10        # Number of basis functions
#' n <- 100       # Number of spatial locations
#' s <- seq(0.1, 0.9, l = n)
#' K <- make_kernel(s, L, "br", gvar = 1, gscl = 1)
#' matplot(s, exp(K), type = "l")
make_kernel <-
  function(obs_coord,
           nbasis,
           type = c("smith", "br"),
           knot_coord = NULL,
           kern_bw = NULL,
           gvar = NULL,
           gscl = NULL) {
    type <- match.arg(type)
    if (type == "smith") {
      dists <- calc_dist(obs_coord, knot_coord)
      K <- dnorm(dists, sd = kern_bw, log = TRUE)
      K <-
        K - matrix(log(apply(K, 1, function(x) {
          sum(exp(x))
        })),
        nrow = nrow(dists),
        ncol = ncol(dists))
      return(K)
    }
    else{
      if (is.matrix(obs_coord)) {
        x <- obs_coord[, 1]
        y <- obs_coord[, 2]
        nobs <- nrow(obs_coord)
      }
      else{
        x <- obs_coord
        y <- NULL
        nobs <- length(x)
      }
      model <- RandomFields::RMexp(var = gvar, scale = gscl)
      RandomFields::RFoptions(spConform = FALSE)
      z.gaus <-
        RandomFields::RFsimulate(
          n = nbasis - 1,
          model = model,
          x = x,
          y = y
        )
      K <- matrix(0, nrow = nobs, nbasis)
      K[,-1] <- z.gaus
      K <-
        K - matrix(log(apply(K, 1, function(x) {
          sum(exp(x))
        })), nrow = nobs, ncol = nbasis)
      return(K)
    }
  }


#' Simulate from stable-mixture model
#' @description Simulate n realizations from stable-mixture model
#' @export
#' @param n number of independent replicates to generate (e.g., years)
#' @param obs_coord  matrix of observation coordinates. The rows correspond
#'   to different locations, and the columns correspond to the dimensions of
#'   those coordinates.
#' @param nbasis number of kernel (basis) functions to generate. Kernels will
#' sum to one across all nbasis at each location.
#' @param type (string) type of kernel function to create. Should be one of
#' either c("smith", "br"). For fixed, Gaussian density
#' kernels, use "smith". For a log-Gaussian process kernels, use "br",
#' which stands for Brown-Resnick.
#' @param knot_coord matrix of knot locations. The rows correspond to
#'   different locations. Needed only if \code{type} = 'smith'.
#' @param kern_bw Kernel bandwidth (sd) for 'smith' kernel bases
#' @param gvar Log-Gaussian process variance term ('brown-resnick'). Needed
#' only if \code{type} = 'br'
#' @param gscl Log-Gaussian process scale/range term ('brown-resnick'). Needed
#' only if \code{type} = 'br'
#' @param alpha stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param delta stable delta parameter (scale): \code{delta} > 0.
#' Elsewhere, this is fixed equal to \code{alpha}.
#' @param theta stable theta parameter (exponential tilting parameter). This
#' parameter controls the tail dependence class. \code{theta} = 0 gives
#' asymptotic dependence, and \code{theta} > 0 gives asymptotic independence.
#' @param return_all (logical) TRUE: just return stable mixture process,
#' FALSE: return components
#' @param return_logZ (logical) should return stable mixture realizations be
#' returned on the log-scale?
#' @return list containing:
#'         lZ: (n by number of obs_coord) matrix of (log) stable mixture
#'         process observations.
#'         lA: (n by number of obs_coord) matrix of 1/alpha-norm scaled
#'         kernel functions at each observation location.
#'         lK: (number of obs_coord by number of knot_coord) matrix of
#'         log-kernel functions at each observation location.
#'         lB: (number of bases by n) stable random effects (i.e. kernel scaling
#'         terms).
#' @examples
#' n <- 3                    # Number of independent replicates (years)
#' L <-30                    # Number of random effects per year
#' alpha <- 0.1             # Stable alpha
#' theta <- 0.001            # Stable exp tilt
#' tau <- 0.5                  # Gauss kernel bw (sd)
#' v <- seq(0, 10, l = L)    # Kernel knots
#' s <-seq(3,8, l = 1000)     # Spatial locations
#' Z <- rstabmix(n, obs_coord = s, knot_coord = v,
#'               nbasis = L, alpha = alpha, delta = alpha,
#'               theta = theta, kern_bw = tau)
#' matplot(s, exp(t(Z)), type = "l")
rstabmix <-
  function(n,
           obs_coord,
           knot_coord,
           nbasis,
           alpha,
           delta,
           theta,
           kern_bw = NULL,
           gvar = NULL,
           gscl = NULL,
           type = c("smith", "br"),
           return_all = FALSE,
           return_logZ = TRUE) {
    if (is.array(obs_coord)) {
      nloc <- nrow(obs_coord)
    }
    else{
      if (!is.vector(obs_coord)) {
        stop("obs_coord must be either a vector or a matrix.")
      }
      nloc <- length(obs_coord)
    }
    type <- match.arg(type)
    lB <- log(matrix(
      rhexpstab(nbasis * n, alpha, delta, theta),
      nrow = nbasis,
      ncol = n
    ))
    lK <-
      make_kernel(
        obs_coord,
        nbasis = nbasis,
        type = type,
        knot_coord = knot_coord,
        kern_bw = kern_bw,
        gscl = gscl,
        gvar = gvar
      )
    lA <- mklA(lK, lB, alpha)

    Z <- matrix(
      extRemes::revd(
        n * nloc,
        loc = exp(lA * alpha),
        scale = alpha * exp(lA * alpha),
        shape = alpha,
        type = "GEV"

      ),
      nrow = nrow(lA)
    )
    if (return_logZ) {
      Z <- log(Z)
    }
    if (!return_all) {
      return(Z)
    }
    else{
      if (return_logZ) {
        return(list(
          lZ = Z,
          lA = lA,
          lK = lK,
          lB = lB
        ))
      }
      if (!return_logZ) {
        return(list(
          Z = Z,
          lA = lA,
          lK = lK,
          lB = lB
        ))
      }
    }
  }


#' Stable mixture CDF
#'
#' @param z vector of quantiles
#' @param alpha stable alpha parameter
#' @param theta stable theta parameter (exponential decay)
#' @param basis Matrix of basis functions at observation locations
#' @param kern_bw Kernel bandwidth parameter (if using Gaussian density basis)
#' @param obs_coord Matrix of observation coordinates
#' @param knot_coord Matrix of knot coordinates (if using Gaussian density basis)
#' @param logp (logical) Should the log-cdf be returned?
#'
#' @export
pstabmix <-
  function(z,
           alpha,
           theta,
           basis = NULL,
           kern_bw = NULL,
           obs_coord = NULL,
           knot_coord = NULL,
           logp = FALSE) {
    delta <- alpha
    if (is.null(basis)) {
      basis <- make_kernel(
        obs_coord,
        nbasis = nrow(obs_coord),
        type = "smith",
        knot_coord = knot_coord
      )
    }
    result <-
      -(delta / alpha) * sum(((theta + (basis / z) ^ (1 / alpha)) ^ alpha - theta ^
                                alpha))
    if (logp) {
      return(result)
    }
    else{
      return(exp(result))
    }
  }

#' Stable mixture quantile function
#' @param p vector of probabilities
#' @param alpha stable alpha parameter
#' @param theta stable theta parameter (exponential decay)
#' @param basis Matrix of basis functions at observation locations
#' @examples
#' # Not run
#' n <- 100
#' L <- 10
#' tau <- 1
#' v <- seq(0,1, l = L)
#' s <- seq(0.1, 0.9, l = n)
#' K <- make_kernel(s, L, "smith", knot_coord = v, kern_bw = tau) # Random basis functions
#' matplot(s, K, type = "l")
#' basis <- K[1,]
#' qout <- qstabmix(0.5, 0.5, 0.1, 0, basis)
#' pstabmix(qout$root, 0.5, 0.1, 0, basis)
#' @export
qstabmix <- function(p, alpha, theta, basis) {
  fdiff <- function(x, alpha, theta, basis)
    pstabmix(x, alpha, theta, basis, logp = TRUE) - log(p)
  upper <- find_upper(p, alpha, theta, basis)
  lower <- find_lower(p, alpha, theta, basis)
  uniroot(
    fdiff,
    c(lower, upper),
    alpha = alpha,
    theta = theta,
    basis = basis,
    tol = .Machine$double.eps ^ 0.9
  )
}

#' Calculate quantile function for stable-mixture observation model
#' @description Calculate quantile function by numerically inverting the CDF.
#' Since the marginal distributions are assumed to be the same across years but
#' vary by spatial location, storing the quantile functions can be dramatically
#' faster than performing root-finding for each observation when the number of
#' years is large.
#' @inheritParams rstabmix
#' @param lK matrix of log-basis functions evaluated at observation locations.
#' nrow: number of observation locations. ncol: number of basis functions.
#' @param grd_len number of design points when defining grid for numerical
#' inversion of CDF. Grids of 100 or more have been consistent with root-finding
#' inversion for stable-mixture cdfs.
#' @return
#' R function that takes a single argument (probability 0 < p < 1) and returns
#' corresponding stable-mixture observation model quantile function
#' @export
qstabmix_spline <- function(alpha, theta, lK, grd_len = 100) {
  if (theta == 0) {
    return(function(p) {
      -log(-log(p))
    })
  }
  else{
    zl <- find_lower(0.001, alpha, theta, lK)
    zu <- find_upper(0.999, alpha, theta, lK)
    qs <- seq(zl, zu, l = grd_len)
    ps <- plstabmix(qs, alpha, theta, lK)
    splinefun(ps, qs, method = "monoH.FC")
  }
}

##' Transform marginal distribution from log-stablemix to gev
#'
#' @param lz log-stable mixture observation matrix (rows: independent replicates
#' columns: observation locations)
#' @inheritParams rstabmix
#' @param lK log-kernel basis matrix (rows: locations, columns: different
#' kernel functions)
#' @param mu GEV location vector (one value per observation location)
#' @param sigma GEV scale vector (one value per observation location)
#' @param xi GEV shape vector (one value per observation location)
#'
#' @return Transformed matrix of observations on GEV scale
#' @export
#' @examples
#' alpha <- 0.9                # positive stable index
#' theta <- 0.01               # tilting parameter
#' tau <- 0.1                  # kernel bandwidth
#' n <- 30                     # number of replicates (years)
#' nloc <- 10^2                # number of spatial locations
#' L <- 25                     # number of random effects per year
#' s <- as.matrix(expand.grid(seq(0.1,0.9, l = sqrt(nloc)),
#'                            seq(0.1,0.9, l = sqrt(nloc))))
#' v <- as.matrix(expand.grid(seq(0, 1, l = sqrt(L)),
#'                            seq(0, 1, l = sqrt(L))))
#' lZ <- rstabmix(n, s, v, nbasis = L, alpha = alpha,
#'                delta = alpha, theta = theta,
#'                kern_bw = tau, type = "smith",
#'                gvar = gvar, gscl =gscl)
#' mu <- sigma <- matrix(1, nrow = n, ncol = nloc)
#' xi <- matrix(0, nrow = n, ncol = nloc)
#' y <- lsm2gev(lZ, alpha, theta, lK, mu, sigma, xi)
lsm2gev <- function(lz, alpha, theta, lK, mu, sigma, xi) {
  ps <- plstabmixM(lz, alpha, theta, lK, logp = FALSE)
  if (max(length(mu), length(sigma), length(xi)) == 1) {
    extRemes::qevd(ps, mu, sigma, xi)
  }
  else{
    qevdM(ps, mu, sigma, xi)
  }
}

#' Create list of quantile functions
#' @description Create list of quantile functions for log-stablemix data model
#' elements of the list correspond to a single observation location. List is in
#' same order as rows of log-basis matrix \code{lK}
#' @inheritParams rstabmix
#' @param lK log-kernel basis matrix (rows: locations, columns: different
#' kernel functions)
#'
#' @return List of log-stable mixture observation model quantile functions
#' @export
create_qlstabmix <- function(alpha, theta, lK) {
  if (is.matrix(lK)) {
    return <- apply(lK, MARGIN = 1, function(x, alpha, theta) {
      qstabmix_spline(alpha, theta, x)
    }, alpha = alpha, theta = theta)
  }
  else{
    qstabmix_spline(alpha, theta, lK)
  }
}

#' Transform marginal distribution from log-stablemix to gev
#' @description Transform from gev to log-stablemix.
#' Need to run create_qlstabmix first.
#' @param y (matrix) gev scale data to transform. Columns: locations,
#' rows: replicates.
#' @param mu GEV location parameter
#' @param sigma GEV scale parameter
#' @param xi GEV shape parameter
#' @param qfuns (list of log-stablemix quantile functions) for each column of y.
#' See \code{create_qlstabmix}
#'
#' @return Matrix of same dimension as \code{y} of transformed GEV data to
#' log-stable mixture scale
#' @export
gev2lsm <- function(y, mu, sigma, xi, qfuns) {
  if (max(length(mu), length(sigma), length(xi)) == 1) {
    ps <- extRemes::pevd(y, mu, sigma, xi)
  }
  else{
    ps <- pevdM(y, mu, sigma, xi)
  }
  if (is.matrix(ps)) {
    for (i in 1:ncol(ps)) {
      ps[, i] <- qfuns[[i]](ps[, i])
    }
  }
  else{
    ps <- qfuns(ps)
  }
  return(ps)
}

#' Wrapper for create_qlstabmix and gev2lsm
#' @description  All-in-one transformation from GEV to log-stabmix
# when the columns of y have different bases
#' @inheritParams gev2lsm
#' @param alpha stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param theta stable theta parameter (exponential tilting parameter). This
#' parameter controls the tail dependence class. \code{theta} = 0 gives
#' @return Matrix of same dimension as \code{y} of transformed GEV data to
#' log-stable mixture scale
#' @export
gev2lsm_all <- function(y, alpha, theta, lK, mu, sigma, xi) {
  qfuns <- create_qlstabmix(alpha, theta, lK)
  gev2lsm(y, mu, sigma, xi, qfuns)
}

#' Simulate Gaussian process
#' @description Simulate a Gaussian process (GP)with exponential covariance
#' function using RandomFields pacakge.
#' @param s (nloc x 2 dimensional matrix) of locations at which to simulate GP
#' @param gvar GP variance parameter
#' @param gscl GP scale/range parameter
#' @export
#' @examples
#' nloc = 100
#' s <- as.matrix(expand.grid(seq(0.1,0.9, l = sqrt(nloc)),
#'                            seq(0.1,0.9, l = sqrt(nloc))))
#' Z <- rgp(s, 1, 1)
rgp <- function(s, gvar, gscl) {
  model <- RandomFields::RMexp(var = gvar, scale = gscl)
  RandomFields::RFoptions(spConform = FALSE)
  RandomFields::RFsimulate(
    n = 1,
    model = model,
    x = s[, 1],
    y = s[, 2]
  )
}

#' Simulate marginal GEV parameter GP processes
#' @description Simulate a Gaussian process with exponential covariance
#' and known mean that may depend on covariates. Used for generating
#' marginal GEV parameters.
#' @param formula R formula for GP mean function using covariates in \code{data}
#' @param data data.frame containing covariates specified in formula
#' @param betas vector of coefficients with which to linearly combine covariates
#' @param nrep number of data replicates of stable mixture process with common
#' margins across time.
#' @param s nloc x 2 dimensional matrix of observation locations
#' @param cov_par list of GP covariance function parameters. Must contain
#' named gvar and gscl elements. See \code{rgp}.
#' @param log.link (logical) should a log link function be used. Useful for GEV
#' scale parameter, which must be positive.
#' @param return_all (logical) should the mean-zero residual process (to the
#' known mean) be returned in addition to the GP with non-zero mean?
#'
#' @export
sample_marpar_gp <-
  function(formula,
           data,
           betas,
           nrep,
           s,
           cov_par = NULL,
           log.link = FALSE,
           return_all = FALSE) {
    nloc <- nrow(s)
    X <- model.matrix(formula, data = data)
    par <- matrix(X %*% betas,
                  nrow = nrep,
                  ncol = nloc,
                  byrow = TRUE)
    if (!is.null(cov_par)) {
      gp <- rgp(s, cov_par[["gvar"]], cov_par[["gscl"]])
      par <- par + matrix(gp,
                          nrow = nrep,
                          ncol  = nloc,
                          byrow = TRUE)
    }
    if (log.link) {
      par <- exp(par)
    }
    if (return_all) {
      return(list(cent_gp = gp, shift_gp = par))
    }
    else{
      return(par)
    }
  }


#' Subtract known mean from Gaussian process to get residual (mean zero)
#' process.
#'
#' @param mpar Matrix containing gaussian process copies in all rows
#' @param X design matrix of known covariates for GP mean
#' @param betas linear coefficients for X matrix
#' @param nrep number of rows in mpar
#' @param log.link (logical) was a log-link used to construct the mpar GP
#'
#' @export
get_gp <- function(mpar, X, betas, nrep, log.link = FALSE) {
  if (log.link) {
    mpar <- log(mpar)
  }
  c(mpar[1, ] - matrix(X %*% betas, nrow = nrep, byrow = TRUE)[1, ])
}

#' Add mean function to mean-zero GP
#' @description Combine mean function with a mean-zero Gaussian process
#' @param X matrix of known covariates used in mean function
#'
#' @param beta linear coefficients for X matrix
#' @param gp mean zero gaussian process
#' @param nrep number of data replicates of stable mixture process that this GP
#' will later interface with
#' @param nloc number of spatial locations
#' @param log.link (logical) should a log link be used?
#'
#' @export
calc_mpar <- function(X, beta, gp, nrep, nloc, log.link = FALSE) {
  par <-
    matrix(X %*% beta,
           nrow = nrep,
           ncol  = nloc,
           byrow = TRUE) + matrix(gp,
                                  nrow = nrep,
                                  ncol  = nloc,
                                  byrow = TRUE)
  if (log.link) {
    par <- exp(par)
  }
  return(par)
}


#' Make conditional draws from stable-mixture dist'n
#' @description Make a draw from stable-mixture conditional on the basis
#' functions, scaling factors, etc.
#' @param n number of replicates to simulate
#'
#' @param obs_coord ncoord x 2 dimensional matrix of observation coordinates
#' @param alpha stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param lK ncoord x nbasis dimensional matrix of log-transformed normalized
#' kernels at each observation location. This should come from mcmc posterior
#' draws.
#' @param lB nreplicate x nyear matrix of log-transformed basis scaling factors.
#' This should come from mcmc posterior draws.
#' @param sub_index Indices of observation locations if only a subset of
#' locations are desired
#' @param return_logZ (logical) should draws of stable-mixture variates be
#' returned on the log-scale?
#'
#' @export
rpostpred <-
  function(n,
           obs_coord,
           alpha,
           lK,
           lB,
           sub_index = NULL,
           return_logZ = TRUE) {
    if (is.array(obs_coord)) {
      nloc <- nrow(obs_coord)
    }
    else{
      if (!is.vector(obs_coord)) {
        stop("obs_coord must be either a vector or a matrix.")
      }
      nloc <- length(obs_coord)
    }
    lA <- mklA(lK, lB, alpha)
    if (!is.null(sub_index)) {
      lA <- lA[sub_index]
      Z <- extRemes::revd(
        length(lA),
        loc = exp(lA * alpha),
        scale = alpha * exp(lA * alpha),
        shape = alpha,
        type = "GEV"
      )
    }
    else{
      Z <- matrix(
        extRemes::revd(
          n * nloc,
          loc = exp(lA * alpha),
          scale = alpha * exp(lA * alpha),
          shape = alpha,
          type = "GEV"
        ),
        nrow = nrow(lA)
      )
    }
    if (return_logZ) {
      Z <- log(Z)
    }
    return(Z)
  }

#' log-stable mixture observation model quantile estimation
#' @description Numerical inversion of log-stablemix observation model CDF using
#'  root finding for a single location
#' @param p desired quantile (0 < p < 1)
#' @param alpha stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param theta stable theta parameter (exponential tilting parameter). This
#' parameter controls the tail dependence class. \code{theta} = 0 gives
#' asymptotic dependence, and \code{theta} > 0 gives asymptotic independence.
#' @param basis Vector of sum=1 constrained basis functions at single
#' observation location.
#'@export
qlstabmix_rf <- function(p, alpha, theta, basis) {
  fdiff <- function(x, alpha, theta, basis)
    plstabmix(x, alpha, theta, basis, logp = TRUE) - log(p)
  upper <- find_upper(p, alpha,  theta, basis)
  lower <- find_lower(p, alpha,  theta, basis)
  uniroot(
    fdiff,
    c(lower, upper),
    alpha = alpha,
    theta = theta,
    basis = basis,
    tol = .Machine$double.eps ^ 0.9
  )
}

#' Simulate from stable-mixture model with a fixed seed
#' @description Simulate n realizations from stable-mixture model using a fixed
#' seed in order to visualize the effects of various parameters
#' @export
#' @param n number of independent replicates to generate (e.g., years)
#' @param obs_coord  matrix of observation coordinates. The rows correspond
#'   to different locations, and the columns correspond to the dimensions of
#'   those coordinates.
#' @param nbasis number of kernel (basis) functions to generate. Kernels will
#' sum to one across all nbasis at each location.
#' @param type (string) type of kernel function to create. Should be one of
#' either c("smith", "br"). For fixed, Gaussian density
#' kernels, use "smith". For a log-Gaussian process kernels, use "br",
#' which stands for Brown-Resnick.
#' @param knot_coord matrix of knot locations. The rows correspond to
#'   different locations. Needed only if \code{type} = 'smith'.
#' @param kern_bw Kernel bandwidth (sd) for 'smith' kernel bases
#' @param gvar Log-Gaussian process variance term ('brown-resnick'). Needed
#' only if \code{type} = 'br'
#' @param gscl Log-Gaussian process scale/range term ('brown-resnick'). Needed
#' only if \code{type} = 'br'
#' @param alpha stable alpha parameter: 0 < \code{alpha} < 1. Controls the
#' tradeoff between all signal (alpha = 0) and all nugget (alpha = 1).
#' @param delta stable delta parameter (scale): \code{delta} > 0.
#' Elsewhere, this is fixed equal to \code{alpha}.
#' @param theta stable theta parameter (exponential tilting parameter). This
#' parameter controls the tail dependence class. \code{theta} = 0 gives
#' asymptotic dependence, and \code{theta} > 0 gives asymptotic independence.
#' @param return_all (logical) TRUE: just return stable mixture process,
#' FALSE: return components
#' @param return_logZ (logical) should return stable mixture realizations be
#' returned on the log-scale?
#' @return list containing:
#'         lZ: (n by number of obs_coord) matrix of (log) stable mixture
#'         process observations.
#'         lA: (n by number of obs_coord) matrix of 1/alpha-norm scaled
#'         kernel functions at each observation location.
#'         lK: (number of obs_coord by number of knot_coord) matrix of
#'         log-kernel functions at each observation location.
#'         lB: (number of bases by n) stable random effects (i.e. kernel scaling
#'         terms).
#' @export
rstabmix_fix_seed <-
  function (n,
            obs_coord,
            knot_coord,
            nbasis,
            alpha,
            delta,
            theta,
            kern_bw = NULL,
            gvar = NULL,
            gscl = NULL,
            type = c("smith",
                     "br"),
            return_all = FALSE,
            return_logZ = TRUE,
            seed = 1)
  {
    if (is.array(obs_coord)) {
      nloc <- nrow(obs_coord)
    }
    else {
      if (!is.vector(obs_coord)) {
        stop("obs_coord must be either a vector or a matrix.")
      }
      nloc <- length(obs_coord)
    }
    type <- match.arg(type)
    set.seed(seed)
    lB <- log(matrix(
      rhexpstab(nbasis * n, alpha, delta, theta),
      nrow = nbasis,
      ncol = n
    ))
    set.seed(seed)
    lB[] <- lB[order(lB)]
    lB[] <-
      lB[sample(1:length(order(lB)), size = length(order(lB)), rep = F)]
    set.seed(seed)
    RandomFields::RFoptions(seed = seed)
    lK <- make_kernel(
      obs_coord,
      nbasis = nbasis,
      type = type,
      knot_coord = knot_coord,
      kern_bw = kern_bw,
      gscl = gscl,
      gvar = gvar
    )
    lA <- mklA(lK, lB, alpha)
    Z <- matrix(
      extRemes::revd(
        n * nloc,
        loc = exp(lA * alpha),
        scale = alpha * exp(lA * alpha),
        shape = alpha,
        type = "GEV"
      ),
      nrow = nrow(lA)
    )
    if (return_logZ) {
      Z <- log(Z)
    }
    if (!return_all) {
      return(Z)
    }
    else {
      if (return_logZ) {
        return(list(
          lZ = Z,
          lA = lA,
          lK = lK,
          lB = lB
        ))
      }
      if (!return_logZ) {
        return(list(
          Z = Z,
          lA = lA,
          lK = lK,
          lB = lB
        ))
      }
    }
  }
