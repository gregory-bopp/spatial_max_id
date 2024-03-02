#' See Stephenson (2009) equation (2)
#'@export
afunR <- function(x, alpha) {
  ((sin(alpha * x) / sin(x)) ^ (1 / (1 - alpha))) * (sin((1 - alpha) * x) /
                                                       sin(alpha * x))
}

#' See Stephenson (2009) equation (2) integrand
#'@export
dhintR <- function(u, s, alpha) {
  (alpha / (1 - alpha)) * (1 / s) ^ (1 / (1 - alpha)) * afunR(pi * u, alpha) *
    exp(-(1 / s) ^ (alpha / (1 - alpha)) * afunR(pi * u, alpha))
}

#' See Stephenson (2009) equation (2) integral for single observation
#'@export
dps_oneR <- function(x, alpha, return_log = TRUE) {
  val <- tryCatch(
    integrate(
      dhintR,
      lower =  0,
      upper = 1,
      s = x,
      alpha = alpha
    )$value,
    error = function(e)
      return(0)
  )
  if (return_log == TRUE) {
    return(log(val))
  }
  else{
    return(val)
  }
}

#' Vectorized version of dps_oneR
#' @export
dpsR <- Vectorize(dps_oneR, vectorize.args = "x")


#' Gamma function, which is used in parameterization of positive stable
#' distribution in the \code{stabledist} package
#' @export
gamf <- function(alpha) {
  iu <- complex(real = 0, imaginary = 1)
  return(abs(1 - iu * tan(pi * alpha / 2)) ^ (-1 / alpha))
}

#' Random positive stable distribution simulation wrapper of stabledist
#' @param n number of realizations to simulate
#' @param alpha positive stable index (0 < alpha < 1)
#' @export
#' @return realizations from positive stable distribution
#' @examples
#' a <- rps(100, 0.5)
rps <- function(n, alpha) {
  stabledist::rstable(
    n,
    alpha = alpha,
    beta = 1,
    gamma = gamf(alpha),
    delta = 0,
    pm = 1
  )
}

#' Log-likelihood of positive stable distribution
#' @param x positive stable observation vector
#' @inheritParams rps
#' @param return_log (logical) should log-likelihood be returned?
#' @param tol numeric tolerance for alpha parameter (away from zero and 1)
#' @return log-likelihood of positive stable distribution
#' @export
#' @examples
#' a <- rps(100, 0.5)
#' lps(a, 0.5)
lps <- function(x,
                alpha,
                return_log = T,
                tol = 1e-2) {
  if ((alpha > tol) & (alpha < 1 - tol)) {
    l <- sum(dps(x, alpha))
  }
  else{
    l <- -Inf
  }
  if (return_log == TRUE) {
    return(l)
  }
  else{
    return(exp(l))
  }
}

#' R implementation of exponentially tilted positive stable (log-)likelihood
#' @param x exponentially tilted positive stable observation vector
#' @param alpha positive stable index (0 < alpha < 1)
#' @param delta Hougaard scale parameter
#' @param theta Hougaard exponentially tilting parameter
#' @param return_log (logical) should log-likelihood be returned
#' @param tol tolerance for alpha parameter (bounding away from 0 and 1)
#' @return exponentially tilted positive stable log-likelihood
#' @export
lhexpstabR <-
  function(x,
           alpha,
           delta,
           theta,
           return_log = TRUE,
           tol = 1e-2) {
    if (alpha > 1 - tol | alpha < tol | delta <= 0 | theta < 0) {
      l <- -Inf
    }
    else{
      n <- length(x)
      beta <- (delta / alpha) ^ (1 / alpha)
      l <- lps(x / beta, alpha) - sum(theta * x) -
        n * (log(beta) - (beta * theta) ^ alpha)
    }
    if (return_log == TRUE) {
      return(l)
    }
    else{
      return(exp(l))
    }
  }

#' R implementation of density function for exponentially tilted PS(alpha)
#' @inheritParams lhexpstabR
#' @return exponentially tilted positive stable density
#' @export
dhexpstabR <-
  function(x,
           alpha,
           delta,
           theta,
           return_log = TRUE,
           tol = 1e-2) {
    if (alpha > 1 - tol | alpha < tol | delta <= 0 | theta < 0) {
      l <- -Inf
    }
    else{
      beta <- (delta / alpha) ^ (1 / alpha)
      l <- dps(x / beta, alpha) - (theta * x) -
        (log(beta) - (beta * theta) ^ alpha)
    }
    if (return_log == TRUE) {
      return(l)
    }
    else{
      return(exp(l))
    }
  }


#' sinc function
#' @description Apply the Method of Luc Devroye (2009) for exponentially
#'              tilted stable distributions
#' @param x sinc argument
#' @param tol tolerance to bound from zero
#' @examples
#' x <- seq(1e-20, 0.01, l= 1000)
#' s <- sapply(x, sinc)
#' plot(x, s, type = "l")
#' @export
sinc <- function(x, tol = 1e-10) {
  if (abs(x) > tol) {
    return(sin(x) / x)
  }
  else{
    return(1)
  }
}

#' B() Zolotarev function (see Devroye (2009) section 9)
#' @export
bzolo <- function(x, alpha) {
  sinc(x) / ((sinc(alpha * x) ^ alpha) * (sinc((1 - alpha) * x) ^ (1 - alpha)))
}

#' A() Zolotarev function (see Devroye (2009) section 9)
#' @export
azolo <- function(x, alpha) {
  (bzolo(x, alpha) / ((alpha ^ alpha) * (1 - alpha) ^ (1 - alpha))) ^ (1 /
                                                                         (alpha - 1))
}

#' Zolotarev density (see Devroye (2009) section 9)
#' @export
dzolo <- function(x, alpha, b) {
  if ((alpha <= 0) | (alpha > 1)) {
    stop("out of parameter space")
  }
  if ((x >= 0) & (x <= pi)) {
    C <- gamma(1 + b * alpha) * gamma(1 + b * (1 - alpha)) /
      (pi * gamma(1 + b))
    return(((bzolo(x, alpha) / ((alpha ^ alpha) * (1 - alpha) ^ (1 - alpha)
    )) ^ b) * C)
  }
  else{
    return(0)
  }
}

#' Simulate single Hougaard RV using Devroye (2009) in reparameterized form
#' @description Apply the Method of Luc Devroye (2009) for exponentially
#' tilted stable distributions from the paper
#' `Random Variate Generation for Exponentially and
#' Polynomially Tilted Stable Distributions'. This method applies only to
#' the theta > 0 case. Here parameterized using the Devroye (2009) lambda
#' parameter. See section 9 of Devroye (2009). This algorithm is useful for
#' large theta.
#' @inheritParams rps
#' @param lambda = theta * (alpha/delta) ^ (1 / alpha)
#' @return single realization from Hougaard dist'n
#' @export
rexpstab_one <- function(alpha, lambda) {
  if ((alpha <= 0) | (alpha > 1))
    stop("alpha needs to be in 0 and 1")
  if (lambda < 0)
    stop("lambda needs to be >= 0")
  gamma <- (lambda ^ alpha) * alpha * (1 - alpha)
  xi <- ((2 + sqrt(pi / 2)) * sqrt(2 * gamma) + 1) / pi
  psi <- (1 / pi) * exp(-(gamma * pi ^ 2) / 8) *
    (2 + sqrt(pi / 2)) * sqrt(gamma * pi)
  w1 <- xi * sqrt(pi / (2 * gamma))
  w2 <- 2 * psi * sqrt(pi)
  w3 <- xi * pi
  b <- (1 - alpha) / alpha
  while (TRUE) {
    U <- pi + 1
    Z <- 2
    while ((U >= pi) || (Z > 1)) {
      V <- runif(1)
      Wp <- runif(1)
      if (gamma >= 1) {
        if (V < w1 / (w1 + w2)) {
          U <- abs(rnorm(1)) / sqrt(gamma)
          if (U >= pi)
            next
        }
        else{
          U <- pi * (1 - Wp ^ 2)
        }
      }
      else{
        if (V < w3 / (w3 + w2)) {
          U <- pi * Wp
        }
        else{
          U <- pi * (1 - Wp ^ 2)
        }
      }
      W <- runif(1)
      zeta <- sqrt(bzolo(U, alpha))
      phi <- (sqrt(gamma) + alpha * zeta) ^ (1 / alpha)
      z <- phi / (phi - (sqrt(gamma)) ^ (1 / alpha))
      rho1 <-
        ifelse((U >= 0) & (gamma >= 1),
               xi * exp(-gamma * (U ^ 2) / 2), 0)
      rho2 <- ifelse((0 < U) & (U < pi), psi / sqrt(pi - U), 0)
      rho3 <- ifelse((0 <= U) & (U <= pi) & (gamma < 1), xi, 0)
      rho <-
        pi * exp(-(lambda ^ alpha) * (1 - (zeta ^ (-2)))) *
        (rho1 + rho2 + rho3) /
        ((1 + sqrt(pi / 2)) * sqrt(gamma) / zeta + z)
      Z <- W * rho
    }
    a <- azolo(U, alpha)
    m <- (b * lambda / a) ^ alpha
    d <- sqrt(m * alpha / a)
    a1 <- d * sqrt(pi / 2)
    a2 <- d
    a3 <- z / a
    s <- a1 + a2 + a3
    Vp <- runif(1)
    Np <- rnorm(1)
    Ep <- rexp(1)
    if (Vp < a1 / s) {
      X <- m - d * abs(Np)
    }
    else{
      if (Vp < (a1 + a2) / s) {
        X <- runif(1, min = m, max = m + d)
      }
      else{
        X <- m + d + Ep * a3
      }
    }
    E <- rexp(1)
    if ((X >= 0) &
        (a * (X - m) + lambda * (X ^ (-b) - m ^ (-b)) -
         ifelse(X < m, (Np ^ 2) / 2, 0) -
         ifelse(X > m + d, Ep, 0) <= E)) {
      return(1 / (X ^ b))
    }
  }
}

#' Simulate from Hougaard distribution
#' @param n number of realizations to simulate
#' @inheritParams rps
#' @param delta Hougaard scale parameter
#' @param theta Hougaard exponential tilting parameter
#' @export
#' @return n realizations from Hougaard(alpha, delta, theta)
#' @export
rhexpstab <- function(n, alpha, delta, theta) {
  if (theta > 0) {
    y <- 1:n
    for (i in 1:n) {
      th <- delta / alpha
      lambda <- theta * th ^ (1 / alpha)
      X <- rexpstab_one(alpha, lambda)
      y[i] <- X * th ^ (1 / alpha)
    }
    return(y)
  }
  else if (theta == 0) {
    rps(n, alpha) * (delta / alpha) ^ (1 / alpha)
  }
}

#' Monte Carlo estimate of Laplace transform of Hougaard distribution
#' @description used for a sanity check to compare against known Laplace transform
#' @export
laplace_mc <- function(t, alpha, delta, theta, n = 10000) {
  y <- rhexpstab(n, alpha, delta, theta)
  mean(exp(-y * t))
}

#' Exact Laplace transform of Hougaard distribution
#' @export
exact_lap <- function(t, alpha, delta, theta) {
  exp(-(delta / alpha) * ((theta + t) ^ alpha - theta ^ alpha))
}

#' Monte Carlo estimate of Hougaard mean (for theta > 0)
#' @export
mc_mean <- function(alpha, delta, theta, n = 10000) {
  y <- rhexpstab(n, alpha, delta, theta)
  mean(y)
}

#' Exact mean of Hougaard distribution (for theta > 0)
#' @export
exact_mean <- function(alpha, delta, theta) {
  delta * theta ^ (alpha - 1)
}
