#' Calculate theoretical chi_u in Gaussian density basis case
#' @param nq number of quantiles at which to calculate chi_u. This is only used
#' if \code{us} is not passed.
#' @param qlim Range of u over which chi_u should be calculated.
#' @param alpha Positive stable index
#' @param theta Hougaard exponential tilting parameter
#' @param obs_coord Observation coordinate matrix/vector defining two locations
#' for which chi_u should be calculated.
#' @param knot_coord Gaussian density basis function knot coordinates
#' @param kern_bw Gaussian density basis kernel bandwidth (sd)
#' @param lK 2 x nbasis matrix of log-transformed, normalized, log-Gaussian
#' density basis functions
#' @param us vector of u values to calculate chi_u for
#' @export
#' @return chi_u, chibar_u and u
fixed_chi_exact <- function(nq,
                            qlim,
                            alpha,
                            theta,
                            obs_coord,
                            knot_coord,
                            kern_bw,
                            lK = NULL,
                            us = NULL) {
  if (is.null(us)) {
    us <- seq(qlim[1], qlim[2], l = nq)
  }
  else{
    nq <- length(us)
  }
  if (is.null(lK)) {
    if (is.matrix(knot_coord)) {
      L <- nrow(knot_coord)
    }
    else{
      L <- length(knot_coord)
    }
    lK <- make_kernel(obs_coord,
                      L,
                      "smith",
                      knot_coord = knot_coord,
                      kern_bw = kern_bw)
  }
  else{
    L <- length(lK[1, ])
  }
  basis1 <- lK[1,]
  basis2 <- lK[2,]
  chi <- chibar <- rep(NA, nq)
  for (i in 1:nq) {
    if (theta == 0) {
      q1 <- q2 <- -log(-log(us[i]))
    }
    else{
      q1 <- qlstabmix_rf(us[i], alpha, theta, basis1)$root
      q2 <- qlstabmix_rf(us[i], alpha, theta, basis2)$root
    }
    if (theta == 0) {
      lbicdf <-
        -sum((
          exp(-q1 / alpha + basis1 / alpha) + exp(-q2 / alpha + basis2 / alpha)
        ) ^ alpha)
    }
    else{
      lbicdf <-  (L * theta ^ alpha -
                    sum((theta +
                           exp((1 / alpha) * (basis1 - q1)
                           )   +
                           exp((1 / alpha) * (basis2 - q2)
                           )) ^ alpha))
    }
    chi[i] <- 2 - (1 - exp(lbicdf)) / (1 - us[i])
    chibar[i] <-
      2 * log(1 - us[i]) / log(1 - 2 * us[i] + exp(lbicdf)) - 1
  }
  obj <- list()
  obj$chi <- chi
  obj$chibar <- chibar
  obj$us <- us
  return(obj)
}

#' Calculate conditional chi_u conditionally on basis functions
#' @inheritParams fixed_chi_exact
#' @param lK 2 x nbasis matrix of log-transformed, normalized, basis functions
#' @param gvar GP variance
#' @param gvar GP exponential cov. range/scale
#' @export
#' @return chi_u, chibar_u and u
lnorm_cond_chi_exact <- function(nq,
                                 qlim,
                                 alpha,
                                 theta,
                                 obs_coord,
                                 L,
                                 gvar,
                                 gscl,
                                 lK = NULL,
                                 us = NULL) {
  if (is.null(us)) {
    us <- seq(qlim[1], qlim[2], l = nq)
  }
  else{
    nq <- length(us)
  }
  if (is.null(lK)) {
    lK <- make_kernel(obs_coord,
                      L,
                      "br",
                      gvar = gvar,
                      gscl = gscl)
  }
  basis1 <- lK[1,]
  basis2 <- lK[2,]
  chi <- chibar <- rep(NA, nq)
  for (i in 1:nq) {
    q1 <- qlstabmix_rf(us[i], alpha, theta, basis1)$root
    q2 <- qlstabmix_rf(us[i], alpha, theta, basis2)$root
    lbicdf <-  (L * theta ^ alpha -
                  sum((theta +
                         exp((1 / alpha) * (basis1 - q1)
                         )   +
                         exp((1 / alpha) * (basis2 - q2)
                         )) ^ alpha))
    chi[i] <- 2 - (1 - exp(lbicdf)) / (1 - us[i])
    chibar[i] <-
      2 * log(1 - us[i]) / log(1 - 2 * us[i] + exp(lbicdf)) - 1
  }
  obj <- list()
  obj$chi <- chi
  obj$chibar <- chibar
  obj$us <- us
  return(obj)
}

#' Wrapper for fixed_chi_exact to calculate chi_u by iterating over several parameters
#' @inheritParams fixed_chi_exact
#' @param alphas vector of positive stable indices for which to calculate chi_u
#' @param thetas Hougaard exponential tilting parameters for which to calculate
#' chi_u
#' @return data.frames of chi_u and chibar_u for different parameter settings
#' @export
prepare_fixed_chis <-
  function(nq,
           qlim,
           obs_coord,
           knot_coord,
           kern_bw,
           alphas,
           thetas) {
    pars <- expand.grid(alphas, thetas)
    colnames(pars) <- c("alpha", "theta")
    chibarm <- chim <- matrix(NA, nrow = nrow(pars), ncol = nq)
    for (i in 1:nrow(pars)) {
      alpha <- pars[i, "alpha"]
      delta <- pars[i, "alpha"]
      theta <- pars[i, "theta"]
      chis <-
        fixed_chi_exact(nq, qlim, alpha, theta, obs_coord, knot_coord, kern_bw)
      chim[i,] <- chis$chi
      chibarm[i, ] <- chis$chibar
    }
    us <- seq(qlim[1], qlim[2], l = nq)
    dus <- data.frame(dummy = paste0("V", 1:nq), u = us)
    chiw <- cbind.data.frame(chim, pars)
    colnames(chiw)[1:nq] <- paste0("V", 1:nq)
    chilong <-
      tidyr::gather(chiw, dummy, chi, V1:V1000, factor_key = FALSE)
    chid <-
      dplyr::inner_join(chilong, dus, by = "dummy") %>% select(-dummy)


    chibarw <- cbind.data.frame(chibarm, pars)
    colnames(chibarw)[1:nq] <- paste0("V", 1:nq)
    chibarlong <-
      tidyr::gather(chibarw, dummy, chibar, V1:V1000, factor_key = FALSE)
    chibard <-
      dplyr::inner_join(chibarlong, dus, by = "dummy") %>% select(-dummy)
    return(list(chid = chid, chibard = chibard))
  }

#' Monte Carlo estimate of chi_u as a function of distance for log-GP basis
#' @description Wrapper for lnorm_cond_chi_exact to calculate chi_u by iterating
#' over several parameters. Monte Carlo is used to marginalize over basis
#' functions.
#' @inheritParams fixed_chi_exact
#' @param alphas vector of positive stable indices for which to calculate chi_u
#' @param thetas Hougaard exponential tilting parameters for which to calculate
#' chi_u
#' @param lK 2 x nbasis matrix of log-transformed, normalized, basis functions
#' @param gvar GP variance
#' @param gvar GP exponential cov. range/scale
#' @return data.frames of chi_u and chibar_u for different parameter settings
#' @param L number of basis functions
#' @param nmc Number of Monte Carlo replicates to use.
#' @return data.frames of chi_u and chibar_u for different parameter settings
#' @export
prepare_lnorm_chis <-
  function(nq,
           qlim,
           obs_coord,
           L,
           gvar,
           gscl,
           alphas,
           thetas,
           nmc = 100) {
    pars <- expand.grid(alphas, thetas)
    colnames(pars) <- c("alpha", "theta")
    chibarm <- chim <- matrix(NA, nrow = nrow(pars), ncol = nq)
    for (i in 1:nrow(pars)) {
      chis <- chibars <- matrix(NA, nrow = nmc, ncol = nq)
      alpha <- pars[i, "alpha"]
      delta <- pars[i, "alpha"]
      theta <- pars[i, "theta"]
      for (j in 1:nmc) {
        chicomb <-
          lnorm_cond_chi_exact(nq, qlim, alpha, theta, obs_coord, L, gvar, gscl)
        chis[j, ] <- chicomb$chi
        chibars[j, ] <- chicomb$chibar
      }
      chim[i,] <- apply(chis, 2, mean)
      chibarm[i, ] <- apply(chibars, 2, mean)
    }
    us <- seq(qlim[1], qlim[2], l = nq)
    dus <- data.frame(dummy = paste0("V", 1:nq), u = us)
    chiw <- cbind.data.frame(chim, pars)
    colnames(chiw)[1:nq] <- paste0("V", 1:nq)
    chilong <-
      tidyr::gather(chiw, dummy, chi, V1:V1000, factor_key = FALSE)
    chid <-
      dplyr::inner_join(chilong, dus, by = "dummy") %>% select(-dummy)


    chibarw <- cbind.data.frame(chibarm, pars)
    colnames(chibarw)[1:nq] <- paste0("V", 1:nq)
    chibarlong <-
      tidyr::gather(chibarw, dummy, chibar, V1:V1000, factor_key = FALSE)
    chibard <-
      dplyr::inner_join(chibarlong, dus, by = "dummy") %>% select(-dummy)
    return(list(chid = chid, chibard = chibard))
  }




#' Empirical estimate of joint exceedance probabilities
#' @description Calculate joint exceedance probability for one pair of vectors
#' @param x first vector of observations
#' @param y second vector of observations
#' @param p upper quantile for which to estimate joint exceedance probability
#' @return empirical estimate of joint exceedance probability.
#' @export
calc_jpexcd_one_pair <- function(x, y, p) {
  n <- length(x)
  r1 <- rank(x, na.last = "keep")
  r2 <- rank(y, na.last = "keep")
  sum((r1 / n > p) & (r2 / n > p), na.rm = T) / n
}

#' Calculate joint exceedance probability for all pairs listed in pair_ids
#' @param X Matrix where data replicates are stored in separate columns
#' @param pair_ids 2-column matrix containing column indicies of X matrix. The
#' function \code{calc_jpexcd_one_pair} will be applied to each pair of columns.
#' @inheritParams calc_jpexcd_one_pair
#' @return Mean joint exceedance probability estimate of exceeding pth quantile
#' as estimated by comparing all pairs of observation vectors listed in
#' \code{pair_ids} matrix
#' @export
calc_jpexcd <- function(X, pair_ids, p) {
  jpexcd <- 1:nrow(pair_ids)
  for (i in 1:nrow(pair_ids)) {
    xy <- X[, pair_ids[i, ]]
    cc <- complete.cases(xy)
    jpexcd[i] <- calc_jpexcd_one_pairC(xy[cc, 1], xy[cc, 2], p)
  }
  mean(jpexcd)
}

#' Calculate chi_u for various distances (binned)
#' @param X matrix of data (columns: locations, rows: years)
#' @param coords n-location x 2 dimensional matrix of observation coordinates,
#' corresponding to the columns of X.
#' @param nbin number of equally spaced bins to use
#' @param ps vector of quantiles over which to estimate chi_u
#' @param bins (optional) alternatively, can explicitly pass a vector of bin
#' boundaries
#' @return empirical estimate of chi_u for various binned distances
#' @export
calc_chi_u_by_dist <- function(X, coords, nbin, ps, bins = NULL) {
  D <- calc_dist(coords, coords)
  if (!is.null(bins)) {
    nbin <- length(bins) - 1
  }
  else{
    bins <- seq(min(D[D > 0]), max(D), l = nbin + 1)
  }
  D[!upper.tri(D)] <- NA
  h <- bins[-length(bins)] + diff(bins) / 2
  dfs <- list()
  for (j in 1:length(ps)) {
    p <- ps[j]
    bl <- list()
    for (i in 1:(length(bins) - 1)) {
      bl[[i]] <- which(matrix(
        between(c(D), bins[i], bins[i + 1]),
        nrow = nrow(D),
        ncol = ncol(D)
      ), arr.ind = T)
    }

    jpexcd <- h
    for (k in 1:length(bl)) {
      if (nrow(bl[[k]]) != 0) {
        jpexcd[k] <- calc_jpexcd(X, bl[[k]], p)
      }
      else{
        if (k == 1) {
          jpexcd <- NA
        }
        else{
          jpexcd[k] <-  jpexcd[k - 1]
        }
      }
    }
    chi <- jpexcd / (1 - p)
    dfs[[j]] <- data.frame(h = h, chi = chi, p = p)
  }
  if (length(ps) == 1) {
    return(dfs[[1]])
  }
  else{
    return(do.call(rbind.data.frame, dfs))
  }
}
