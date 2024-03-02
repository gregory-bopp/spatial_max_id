#' Combine MCMC samples from log-Gaussian process basis model
#' @description Since cluster jobs sometimes have max wall times,
#' it may be necessary to break MCMC run into smaller jobs
#' of a few days at a time. This function chains together the output from a
#' collection of smaller jobs that were run sequentially.
#' @param file.dir Directory containing all MCMC datasets (and nothing else)
#' @param mar_gp_which vector of strings containing some subset of mu, sigma, xi
#' indicating which GEV parameters were modeled as Gaussian processes
#' @param burnin How many mcmc iterations should be discarded as burnin, if any?
#' @export
#' @return combined mcmc output list of same format as in stored R files
lnorm_combine_mcmc <-
  function(file.dir,
           mar_gp_which = c("mu", "sigma"),
           burnin = 100) {
    files <- list.files(file.dir, full.names = TRUE)
    nms <- list.files(file.dir, full.names = FALSE)
    ord <- order(as.numeric(unlist(strsplit(nms, ".Rdata"))))
    files <- files[ord]
    for (i in 1:length(files)) {
      load(files[i])
      if (i == 1) {
        smp_par <- out$smp_par
        smp_lB <- out$smp_lB
        smp_lK <- out$smp_lK
        if ("mu" %in% mar_gp_which) {
          smp_mu_gp <- out$smp_mu_gp
          smp_mu_gvar <- out$smp_mu_gvar
          smp_mu_gscl <- out$smp_mu_gscl
        }
        else{
          smp_mu_gp <- NULL
          smp_mu_gvar <- NULL
          smp_mu_gscl <- NULL
        }
        if ("sigma" %in% mar_gp_which) {
          smp_sigma_gp <- out$smp_sigma_gp
          smp_sigma_gvar <- out$smp_sigma_gvar
          smp_sigma_gscl <- out$smp_sigma_gscl
        }
        else{
          smp_sigma_gp <- NULL
          smp_sigma_gvar <- NULL
          smp_sigma_gscl <- NULL
        }
        if ("xi" %in% mar_gp_which) {
          smp_xi_gp <- out$smp_xi_gp
          smp_xi_gvar <- out$smp_xi_gvar
          smp_xi_gscl <- out$smp_xi_gscl
        }
        else{
          smp_xi_gp <- NULL
          smp_xi_gvar <- NULL
          smp_xi_gscl <- NULL
        }
      }
      else{
        smp_par <- rbind(smp_par, out$smp_par)
        smp_lB <- rbind(smp_lB, out$smp_lB)
        smp_lK <- abind::abind(smp_lK, out$smp_lK, along = 1)
        if ("mu" %in% mar_gp_which) {
          smp_mu_gp <- rbind(smp_mu_gp, out$smp_mu_gp)
          smp_mu_gvar <- c(smp_mu_gvar, out$smp_mu_gvar)
          smp_mu_gscl <- c(smp_mu_gscl, out$smp_mu_gscl)
        }
        if ("sigma" %in% mar_gp_which) {
          smp_sigma_gp <- rbind(smp_sigma_gp, out$smp_sigma_gp)
          smp_sigma_gvar <- c(smp_sigma_gvar, out$smp_sigma_gvar)
          smp_sigma_gscl <- c(smp_sigma_gscl, out$smp_sigma_gscl)
        }
        if ("xi" %in% mar_gp_which) {
          smp_xi_gp <- rbind(smp_xi_gp, out$smp_xi_gp)
          smp_xi_gvar <- c(smp_xi_gvar, out$smp_xi_gvar)
          smp_xi_gscl <- c(smp_xi_gscl, out$smp_xi_gscl)
        }
      }
    }
    nmc <- nrow(smp_par)
    if ("mu" %in% mar_gp_which) {
      nloc <- ncol(smp_mu_gp)
      smp_mu <- matrix(nrow =  nmc, ncol = nloc)
    }
    if ("sigma" %in% mar_gp_which) {
      nloc <- ncol(smp_sigma_gp)
      smp_sigma <- matrix(nrow =  nmc, ncol = ncol(smp_sigma_gp))
    }
    if ("xi" %in% mar_gp_which) {
      nloc <- ncol(smp_xi_gp)
      smp_xi <- matrix(nrow =  nmc, ncol = ncol(smp_xi_gp))
    }
    if (any(mar_gp_which %in% c("mu", "sigma", "xi"))) {
      for (i in 1:nmc) {
        if ("mu" %in% mar_gp_which) {
          smp_mu[i, ] <- calc_mpar(matrix(1, nrow = nloc, ncol = 1),
                                   smp_par[i, "mu"],
                                   smp_mu_gp[i, ],
                                   1,
                                   nloc,
                                   log.link = FALSE)
        }
        else{
          smp_mu <- NULL
        }
        if ("sigma" %in% mar_gp_which) {
          smp_sigma[i, ] <- calc_mpar(matrix(1, nrow = nloc, ncol = 1),
                                      smp_par[i, "sigma"],
                                      smp_sigma_gp[i, ],
                                      1,
                                      nloc,
                                      log.link = TRUE)
        }
        else{
          smp_sigma <- NULL
        }

        if ("xi" %in% mar_gp_which) {
          smp_xi[i, ] <- calc_mpar(matrix(1, nrow = nloc, ncol = 1),
                                   smp_par[i, "xi"],
                                   smp_xi_gp[i, ],
                                   1,
                                   nloc,
                                   log.link = FALSE)
        }
        else{
          smp_xi <- NULL
        }
      }
    }
    else{
      smp_mu <- smp_sigma <- smp_xi <- NULL
    }

    return(
      list(
        smp_par = smp_par[-c(1:burnin), ],
        smp_lB = smp_lB[-c(1:burnin), ],
        smp_lK = smp_lK[-c(1:burnin), , ],
        smp_mu_gp = smp_mu_gp[-c(1:burnin), ],
        smp_mu_gvar = smp_mu_gvar[-c(1:burnin)],
        smp_mu_gscl = smp_mu_gscl[-c(1:burnin)],
        smp_sigma_gp = smp_sigma_gp[-c(1:burnin), ],
        smp_sigma_gvar = smp_sigma_gvar[-c(1:burnin)],
        smp_sigma_gscl = smp_sigma_gscl[-c(1:burnin)],
        smp_xi_gp = smp_xi_gp[-c(1:burnin), ],
        smp_xi_gvar = smp_xi_gvar[-c(1:burnin)],
        smp_xi_gscl = smp_xi_gscl[-c(1:burnin)],
        smp_mu = smp_mu[-c(1:burnin), ],
        smp_sigma = smp_sigma[-c(1:burnin), ],
        smp_xi = smp_xi[-c(1:burnin), ]
      )
    )
  }


#' Combine MCMC samples from Gaussian density basis model
#' @description Since cluster jobs typically can only run for a few days at a
#' time, for long jobs, it may be necessary to break MCMC run into smaller jobs
#' of a few days at a time. This function chains together the output from a
#' collection of smaller jobs that were run sequentially.
#' @inheritParams lnorm_combine_mcmc
#' @export
fixed_combine_mcmc <-
  function(file.dir,
           mar_gp_which = c("mu", "sigma"),
           burnin = 100) {
    files <- list.files(file.dir, full.names = TRUE)
    nms <- list.files(file.dir, full.names = FALSE)
    ord <- order(as.numeric(unlist(strsplit(nms, ".Rdata"))))
    files <- files[ord]
    for (i in 1:length(files)) {
      load(files[i])
      if (i == 1) {
        smp_par <- out$smp_par
        smp_lB <- out$smp_lB
        if ("mu" %in% mar_gp_which) {
          smp_mu_gp <- out$smp_mu_gp
          smp_mu_gvar <- out$smp_mu_gvar
          smp_mu_gscl <- out$smp_mu_gscl
        }
        else{
          smp_mu_gp <- NULL
          smp_mu_gvar <- NULL
          smp_mu_gscl <- NULL
        }
        if ("sigma" %in% mar_gp_which) {
          smp_sigma_gp <- out$smp_sigma_gp
          smp_sigma_gvar <- out$smp_sigma_gvar
          smp_sigma_gscl <- out$smp_sigma_gscl
        }
        else{
          smp_sigma_gp <- NULL
          smp_sigma_gvar <- NULL
          smp_sigma_gscl <- NULL
        }
        if ("xi" %in% mar_gp_which) {
          smp_xi_gp <- out$smp_xi_gp
          smp_xi_gvar <- out$smp_xi_gvar
          smp_xi_gscl <- out$smp_xi_gscl
        }
        else{
          smp_xi_gp <- NULL
          smp_xi_gvar <- NULL
          smp_xi_gscl <- NULL
        }
      }
      else{
        smp_par <- rbind(smp_par, out$smp_par)
        smp_lB <- rbind(smp_lB, out$smp_lB)
        if ("mu" %in% mar_gp_which) {
          smp_mu_gp <- rbind(smp_mu_gp, out$smp_mu_gp)
          smp_mu_gvar <- c(smp_mu_gvar, out$smp_mu_gvar)
          smp_mu_gscl <- c(smp_mu_gscl, out$smp_mu_gscl)
        }
        if ("sigma" %in% mar_gp_which) {
          smp_sigma_gp <- rbind(smp_sigma_gp, out$smp_sigma_gp)
          smp_sigma_gvar <- c(smp_sigma_gvar, out$smp_sigma_gvar)
          smp_sigma_gscl <- c(smp_sigma_gscl, out$smp_sigma_gscl)
        }
        if ("xi" %in% mar_gp_which) {
          smp_xi_gp <- rbind(smp_xi_gp, out$smp_xi_gp)
          smp_xi_gvar <- c(smp_xi_gvar, out$smp_xi_gvar)
          smp_xi_gscl <- c(smp_xi_gscl, out$smp_xi_gscl)
        }
      }
    }

    nmc <- nrow(smp_par)
    if ("mu" %in% mar_gp_which) {
      nloc <- ncol(smp_mu_gp)
      smp_mu <- matrix(nrow =  nmc, ncol = nloc)
    }
    if ("sigma" %in% mar_gp_which) {
      nloc <- ncol(smp_sigma_gp)
      smp_sigma <- matrix(nrow =  nmc, ncol = ncol(smp_sigma_gp))
    }
    if ("xi" %in% mar_gp_which) {
      nloc <- ncol(smp_xi_gp)
      smp_xi <- matrix(nrow =  nmc, ncol = ncol(smp_xi_gp))
    }
    if (any(mar_gp_which %in% c("mu", "sigma", "xi"))) {
      for (i in 1:nmc) {
        if ("mu" %in% mar_gp_which) {
          smp_mu[i, ] <- calc_mpar(matrix(1, nrow = nloc, ncol = 1),
                                   smp_par[i, "mu"],
                                   smp_mu_gp[i, ],
                                   1,
                                   nloc,
                                   log.link = FALSE)
        }
        else{
          smp_mu <- NULL
        }
        if ("sigma" %in% mar_gp_which) {
          smp_sigma[i, ] <- calc_mpar(matrix(1, nrow = nloc, ncol = 1),
                                      smp_par[i, "sigma"],
                                      smp_sigma_gp[i, ],
                                      1,
                                      nloc,
                                      log.link = TRUE)
        }
        else{
          smp_sigma <- NULL
        }

        if ("xi" %in% mar_gp_which) {
          smp_xi[i, ] <- calc_mpar(matrix(1, nrow = nloc, ncol = 1),
                                   smp_par[i, "xi"],
                                   smp_xi_gp[i, ],
                                   1,
                                   nloc,
                                   log.link = FALSE)
        }
        else{
          smp_xi <- NULL
        }
      }
    }
    else{
      smp_mu <- smp_sigma <- smp_xi <- NULL
    }

    return(
      list(
        smp_par = smp_par[-c(1:burnin), ],
        smp_lB = smp_lB[-c(1:burnin), ],
        smp_mu_gp = smp_mu_gp[-c(1:burnin), ],
        smp_mu_gvar = smp_mu_gvar[-c(1:burnin)],
        smp_mu_gscl = smp_mu_gscl[-c(1:burnin)],
        smp_sigma_gp = smp_sigma_gp[-c(1:burnin), ],
        smp_sigma_gvar = smp_sigma_gvar[-c(1:burnin)],
        smp_sigma_gscl = smp_sigma_gscl[-c(1:burnin)],
        smp_xi_gp = smp_xi_gp[-c(1:burnin), ],
        smp_xi_gvar = smp_xi_gvar[-c(1:burnin)],
        smp_xi_gscl = smp_xi_gscl[-c(1:burnin)],
        smp_mu = smp_mu[-c(1:burnin), ],
        smp_sigma = smp_sigma[-c(1:burnin), ],
        smp_xi = smp_xi[-c(1:burnin), ]
      )
    )
  }
