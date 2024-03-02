#' Posterior predictive estimat of chi_u as a function of distance
#' @param mh MCMC samples from \code{mcmc_lnorm_basis} or
#' \code{mcmc_fixed_basis} output.
#' @param coords n-location x 2 dimensional matrix of observation coordinates
#' @param nbin number of bins to use (unique distances) when estimating chi_u
#' @param ps vector of probabilities (u terms in chi_u)
#' @param type one of "smith" or "br" for Gaussian density or log-Gaussian
#' process basis respectively.
#' @param knot_coord n-knot x 2 matrix of knot coordinates (needed only if
#' \code{type} = "smith")
#' @param bins (optional) vector of distances defining bin boundaries. If this
#' is not supplied, bins are defined between shortest distances and greatest
#' distance between observation locations.
#' @param thin_int (optional) thinning interval to use.
#' @param n.core number of cores to parallelize calculation over.
#' @export
#' @return Posterior predictive estimate of chi_u for different bin distances,
#' and quantiles defined by \code{ps}.
#'@export
calc_chi_u_by_dist_mcmc <-
  function (mh,
            coords,
            nbin,
            ps,
            type = c("smith", "br"),
            knot_coord = NULL,
            bins = NULL,
            thin_int = 1,
            n.core = NULL)
  {
    type <- match.arg(type)
    np <- length(ps)
    if (type == "smith") {
      L <- nrow(knot_coord)
    }
    else{
      L <- dim(mh$smp_lK)[3]
    }
    D <- calc_dist(coords, coords)
    nmcmc <- nrow(mh$smp_par)
    if (!is.null(bins)) {
      nbin <- length(bins) - 1
    }
    else {
      bins <- seq(min(D[D > 0]), max(D), l = nbin + 1)
    }
    D[!upper.tri(D)] <- NA
    h <- bins[-length(bins)] + diff(bins) / 2
    dfs <- list()
    bl <- list()
    for (i in 1:(length(bins) - 1)) {
      bl[[i]] <- which(matrix(
        between(c(D), bins[i], bins[i +
                                      1]),
        nrow = nrow(D),
        ncol = ncol(D)
      ), arr.ind = T)
    }

    mcmcids <- which((1:nmcmc) %% thin_int == 0)
    if (!is.null(n.core)) {
      n.core = parallel::detectCores()
    }
    doMC::registerDoMC(cores = n.core)
    chidfall <- foreach(i = mcmcids) %dopar% {
      if (type == "smith") {
        lK <- make_kernel(coords,
                          L,
                          "smith",
                          knot_coord = knot_coord,
                          kern_bw = mh$smp_par[i, "tau"])
      }
      for (k in 1:length(bl)) {
        if (nrow(bl[[k]]) != 0) {
          chimat <- matrix(nrow = np, ncol = nrow(bl[[k]]))
          if (type == "smith") {
            for (j in 1:nrow(bl[[k]])) {
              chimat[, j] <-
                fixed_chi_exact(
                  NULL,
                  NULL,
                  mh$smp_par[i, "alpha"],
                  mh$smp_par[i, "theta"],
                  obs_coord = NULL,
                  knot_coord = NULL,
                  kern_bw = NULL,
                  lK = lK[bl[[k]][j, ], ],
                  us = ps
                )$chi
            }
          }
          else{
            for (j in 1:nrow(bl[[k]])) {
              chimat[, j] <-
                lnorm_cond_chi_exact(
                  NULL,
                  NULL,
                  mh$smp_par[i, "alpha"],
                  mh$smp_par[i, "theta"],
                  obs_coord = NULL,
                  L = L,
                  gvar = mh$smp_par[i, "gvar"],
                  gscl = mh$smp_par[i, "gscl"],
                  lK = mh$smp_lK[i, bl[[k]][j, ], ],
                  us = ps
                )$chi
            }
          }
          chi <- apply(chimat, 1, mean, na.rm = T)
        }
        else {
          chi <- rep(NA, np)
        }
        dfs[[k]] <- data.frame(
          h = h[k],
          chi = chi,
          p = ps,
          mcit = i
        )
      }
      do.call(rbind.data.frame, dfs)
    }
    return(chidfall)
  }
