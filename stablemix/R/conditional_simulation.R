#' Conditionally simulate Gaussian process basis functions
#' @param lK normalized log-Gaussian process basis function matrix on which to
#' condition observed at \code{sold} and transformed to log-scale
#' @param snew n-new-location x 2 dimensional coordinate matrix of locations at
#' which to conditionally simulate lK.
#' @param sold n-obs-location x 2 dimensional coordinate matrix of locations at
#' which samples of lK are to be conditioned on
#' @param gvar Gaussian process exponential covariance function variance
#' @param gscl Gassian process exponential covariance function scale/range
#' @export
cond_sim_lK <- function(lK,
                        snew,
                        sold,
                        gvar,
                        gscl) {
  nnold <- nrow(sold)
  nnew <- nrow(snew)
  L <- ncol(lK)
  G <-
    lK[,-1] - matrix(lK[, 1], ncol = L - 1, nrow = nnold) 
  Gimp <- matrix(NA_real_, nrow = nnew, ncol = L - 1)
  model <- RandomFields::RMexp(var = gvar, scale = gscl)
  RandomFields::RFoptions(spConform = FALSE)
  for (i in 1:(L - 1)) {
    Gimp[, i] <- RandomFields::RFsimulate(
      n = 1,
      model = model,
      x = snew[, 1],
      y = snew[, 2],
      data = data.frame(x = sold[, 1], y = sold[, 2], G[, i])
    )
  }
  lKimp <- matrix(0, nrow = nnew, L)
  lKimp[,-1] <- Gimp
  lKimp <-
    lKimp - matrix(log(apply(lKimp, 1, function(x) {
      sum(exp(x))
    })), nrow = nnew, ncol = L)
}


#' Conditionally simulate GEV marginal parameter Gaussian processes
#' @param mpar_gp Marginal parameter Gaussian process observations on which
#' to condition
#' @param nrep number of replicates of stable-mixture observations
#' @param betas linear coefficients in, e.g., loc = X beta
#' @param formula R formula describing form of covariates in GP mean
#' @param data_new data.frame containing covariate values at \code{snew}
#' locations
#' @param log.link (logical) should a log link be used for GEV parameter? This
#' should be set to TRUE when parameter is constrained to be positive. For
#' example in the GEV scale case.
#' @inheritParams cond_sim_lK
#' @export
cond_sim_mpar <-
  function(mpar_gp,
           snew,
           sold,
           gvar,
           gscl,
           nrep = NULL,
           betas = NULL,
           formula = NULL,
           data_new = NULL,
           log.link = FALSE) {
    model <- RandomFields::RMexp(var = gvar, scale = gscl)
    RandomFields::RFoptions(spConform = FALSE)
    mparimp_gp <- RandomFields::RFsimulate(
      n = 1,
      model = model,
      x = snew[, 1],
      y = snew[, 2],
      data = data.frame(x = sold[, 1], y = sold[, 2], mpar_gp)
    )
    if (!is.null(formula)) {
      nnew <- nrow(snew)
      X <- model.matrix(formula, data = data_new)
      mparimp_mean <-
        matrix(X %*% betas,
               nrow = nrep,
               ncol = nnew,
               byrow = TRUE)
      mparimp_gp <-
        matrix(mparimp_gp,
               nrow = nrep,
               ncol = nnew,
               byrow = TRUE) + mparimp_mean
    }
    if (log.link == TRUE) {
      mparimp_gp <- exp(mparimp_gp)
    }
    return(mparimp_gp)
  }


#' Wrapper for cond_sim_lK and cond_sim_mpar to make conditional draws of latent processes
#' @param out MCMC sample output from mcmc_lnorm_basis
#' @param obs_coords nloc x 2 dimensional matrix of observation location coordinates
#' @param pred_coords n-prediction-location x 2 dimensional matrix of
#' prediction location coordinates at which samples of the latent processes
#' are to be conditionally simulated
#' @param thin_int (optional) thinning interval - if only a subset of MCMC
#' samples should be used (e.g. for memory constraints/testing)
#' @param sim_y (logical) should conditional draws of the stable-mixture
#' response be simulated at \code{pred_coords}?
#' @export
lnorm_cond_sim_all <-
  function(out,
           obs_coords,
           pred_coords,
           thin_int = 1,
           sim_y = TRUE) {
    # Inside function ---------------------------------------------------------
    npred <- nrow(pred_coords)
    nmcmc <- nrow(out$smp_par)
    n_smp <-
      ifelse((thin_int == 1),
             yes = nmcmc,
             no = floor(nmcmc / thin_int))

    L <- dim(out$smp_lK)[3]
    nrep <- dim(out$smp_lB)[2] / dim(out$smp_lK)[3]
    dfnew <- data.frame(x = rep(1, npred * nrep))
    lK_smps <- array(NA_real_, dim = c(n_smp, npred, L))
    mu_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)
    sigma_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)
    xi_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)

    if (sim_y) {
      y_smps <- array(NA_real_, dim = c(n_smp, nrep, npred))
    }

    smp_loc <- 1
    for (i in 1:nmcmc) {
      if (i %% thin_int == 0) {
        alpha_cur <- out$smp_par[i, "alpha"]
        theta_cur <- out$smp_par[i, "theta"]
        # lK ----------------------------------------------------------------------
        lKimp <-
          cond_sim_lK(out$smp_lK[i, , ],
                      pred_coords,
                      obs_coords,
                      out$smp_par[i, "gvar"],
                      out$smp_par[i, "gscl"])

        # Mu ----------------------------------------------------------------------
        muimp <-
          cond_sim_mpar(
            out$smp_mu_gp[i, ],
            pred_coords,
            obs_coords,
            out$smp_mu_gvar[i],
            out$smp_mu_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "mu"],
            formula = ~ 1,
            data_new = dfnew
          )

        # Sigma ----------------------------------------------------------------------
        sigmaimp <-
          cond_sim_mpar(
            out$smp_sigma_gp[i, ],
            pred_coords,
            obs_coords,
            out$smp_sigma_gvar[i],
            out$smp_sigma_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "sigma"],
            formula = ~ 1,
            data_new = dfnew,
            log.link = TRUE
          )


        # Xi ----------------------------------------------------------------------
        if (is.null(out$smp_xi_gp)) {
          xiimp <- matrix(out$smp_par[i, "xi"], nrow = nrep, ncol = npred)
        }
        else{
          xiimp <- cond_sim_mpar(
            out$smp_xi_gp[i, ],
            pred_coords,
            obs_coords,
            out$smp_xi_gvar[i],
            out$smp_xi_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "xi"],
            formula = ~ 1,
            data_new = dfnew
          )
        }
        # Response ----------------------------------------------------------------
        # n is nrep
        if (sim_y) {
          lzimp <-
            rpostpred(nrep,
                      pred_coords,
                      alpha_cur,
                      lKimp,
                      matrix(out$smp_lB[i, ], nrow = L))
          yimp <-
            lsm2gev(lzimp,
                    alpha_cur,
                    theta_cur,
                    lKimp,
                    muimp,
                    sigmaimp,
                    xiimp)
          y_smps[smp_loc, , ] <- yimp
        }
        lK_smps[smp_loc, , ] <- lKimp
        mu_smps[smp_loc, ] <- muimp[1, ]
        sigma_smps[smp_loc, ] <- sigmaimp[1, ]
        xi_smps[smp_loc, ] <- xiimp[1, ]
        smp_loc <- smp_loc + 1
      }
    }
    if (sim_y) {
      return(list(
        y = y_smps,
        lK = lK_smps,
        mu = mu_smps,
        sigma = sigma_smps,
        xi = xi_smps
      ))
    }
    else{
      return(list(
        lK = lK_smps,
        mu = mu_smps,
        sigma = sigma_smps,
        xi = xi_smps
      ))
    }
  }


#' Wrapper for cond_sim_mpar to make conditional draws of latent processes in the Gaussian density basis case
#' @param out MCMC sample output from mcmc_fixed_basis
#' @param obs_coords nloc x 2 dimensional matrix of observation location coordinates
#' @param pred_coords n-prediction-location x 2 dimensional matrix of
#' prediction location coordinates at which samples of the latent processes
#' are to be conditionally simulated
#' @param thin_int (optional) thinning interval - if only a subset of MCMC
#' samples should be used (e.g. for memory constraints/testing)
#' @param sim_y (logical) should conditional draws of the stable-mixture
#' response be simulated at \code{pred_coords}?
#' @export
fixed_cond_sim_all <-
  function(out,
           obs_coords,
           pred_coords,
           knot_coords,
           thin_int = 1,
           sim_y = TRUE) {
    # Inside function ---------------------------------------------------------
    npred <- nrow(pred_coords)
    nmcmc <- nrow(out$smp_par)
    n_smp <-
      ifelse((thin_int == 1),
             yes = nmcmc,
             no = floor(nmcmc / thin_int))
    L <- nrow(knot_coords)
    nrep <- ncol(out$smp_lB) / L
    dfnew <- data.frame(x = rep(1, npred * nrep))

    lK_smps <- array(NA_real_, dim = c(n_smp, npred, L))
    mu_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)
    sigma_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)
    xi_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)

    if (sim_y) {
      y_smps <- array(NA_real_, dim = c(n_smp, nrep, npred))
    }

    smp_loc <- 1
    for (i in 1:nmcmc) {
      if (i %% thin_int == 0) {
        alpha_cur <- out$smp_par[i, "alpha"]
        theta_cur <- out$smp_par[i, "theta"]
        # lK ----------------------------------------------------------------------
        lKimp <-
          make_kernel(
            pred_coords,
            nbasis =  L,
            type = "smith",
            knot_coord = knot_coords,
            kern_bw = out$smp_par[i, "tau"]
          )
        # Mu ----------------------------------------------------------------------
        muimp <-
          cond_sim_mpar(
            out$smp_mu_gp[i, ],
            pred_coords,
            obs_coords,
            out$smp_mu_gvar[i],
            out$smp_mu_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "mu"],
            formula = ~ 1,
            data_new = dfnew
          )

        # Sigma ----------------------------------------------------------------------
        sigmaimp <-
          cond_sim_mpar(
            out$smp_sigma_gp[i, ],
            pred_coords,
            obs_coords,
            out$smp_sigma_gvar[i],
            out$smp_sigma_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "sigma"],
            formula = ~ 1,
            data_new = dfnew,
            log.link = TRUE
          )

        # Xi ----------------------------------------------------------------------
        if (is.null(out$smp_xi_gp)) {
          xiimp <- matrix(out$smp_par[i, "xi"], nrow = nrep, ncol = npred)
        }
        else{
          xiimp <- cond_sim_mpar(
            out$smp_xi_gp[i, ],
            pred_coords,
            obs_coords,
            out$smp_xi_gvar[i],
            out$smp_xi_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "xi"],
            formula = ~ 1,
            data_new = dfnew
          )
        }
        # Response ----------------------------------------------------------------
        if (sim_y) {
          lzimp <-
            rpostpred(nrep,
                      pred_coords,
                      alpha_cur,
                      lKimp,
                      matrix(out$smp_lB[i, ], nrow = L))
          yimp <-
            lsm2gev(lzimp,
                    alpha_cur,
                    theta_cur,
                    lKimp,
                    muimp,
                    sigmaimp,
                    xiimp)
          y_smps[smp_loc, , ] <- yimp
        }
        lK_smps[smp_loc, , ] <- lKimp
        mu_smps[smp_loc, ] <- muimp[1, ]
        sigma_smps[smp_loc, ] <- sigmaimp[1, ]
        xi_smps[smp_loc, ] <- xiimp[1, ]
        smp_loc <- smp_loc + 1
      }
    }
    if (sim_y) {
      return(list(
        y = y_smps,
        lK = lK_smps,
        mu = mu_smps,
        sigma = sigma_smps,
        xi = xi_smps
      ))
    }
    else{
      return(list(
        lK = lK_smps,
        mu = mu_smps,
        sigma = sigma_smps,
        xi = xi_smps
      ))
    }
  }


#' Conditional simulation for calculation of chi for Gaussian density basis
#' Wrapper for cond_sim_lK and cond_sim_mpar to make conditional draws of latent processes
#' @param out MCMC sample output from mcmc_fixed_basis
#' @param obs_coords nloc x 2 dimensional matrix of observation location coordinates
#' @param pred_coords n-prediction-location x 2 dimensional matrix of
#' prediction location coordinates at which samples of the latent processes
#' are to be conditionally simulated
#' @param thin_int (optional) thinning interval - if only a subset of MCMC
#' samples should be used (e.g. for memory constraints/testing)
#' @param sim_y (logical) should conditional draws of the stable-mixture
#' response be simulated at \code{pred_coords}?
#' @param knot_coords knot locations of Gaussian density basis
#' @param npred number of predictive draws to make
#' @export
fixed_cond_sim_spatial_chi <-
  function (out,
            obs_coords,
            pred_coords,
            knot_coords,
            thin_int = 1,
            sim_y = TRUE,
            npred = 100)
  {
    pred_sub_coords_list <- list()
    nmcmc <- nrow(out$smp_par)
    n_smp <-
      ifelse((thin_int == 1),
             yes = nmcmc,
             no = floor(nmcmc / thin_int))
    L <- nrow(knot_coords)
    nrep <- ncol(out$smp_lB) / L
    dfnew <- data.frame(x = rep(1, npred * nrep))
    y_smps <- array(NA_real_, dim = c(n_smp, nrep, npred))

    smp_loc <- 1
    for (i in 1:nmcmc) {
      if (i %% thin_int == 0) {
        pred_sub_coords_list[[smp_loc]] <- pred_sub_coords <-
          pred_coords[sample(1:nrow(pred_coords),
                             size = npred,
                             replace = F),]
        alpha_cur <- out$smp_par[i, "alpha"]
        theta_cur <- out$smp_par[i, "theta"]
        lKimp <-
          make_kernel(
            pred_sub_coords,
            nbasis = L,
            type = "smith",
            knot_coord = knot_coords,
            kern_bw = out$smp_par[i,
                                  "tau"]
          )
        muimp <- cond_sim_mpar(
          out$smp_mu_gp[i,],
          pred_sub_coords,
          obs_coords,
          out$smp_mu_gvar[i],
          out$smp_mu_gscl[i],
          nrep = nrep,
          betas = out$smp_par[i, "mu"],
          formula = ~ 1,
          data_new = dfnew
        )
        sigmaimp <- cond_sim_mpar(
          out$smp_sigma_gp[i,],
          pred_sub_coords,
          obs_coords,
          out$smp_sigma_gvar[i],
          out$smp_sigma_gscl[i],
          nrep = nrep,
          betas = out$smp_par[i,
                              "sigma"],
          formula = ~ 1,
          data_new = dfnew,
          log.link = TRUE
        )
        if (is.null(out$smp_xi_gp)) {
          xiimp <- matrix(out$smp_par[i, "xi"], nrow = nrep,
                          ncol = npred)
        }
        else {
          xiimp <- cond_sim_mpar(
            out$smp_xi_gp[i,],
            pred_sub_coords,
            obs_coords,
            out$smp_xi_gvar[i],
            out$smp_xi_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "xi"],
            formula = ~ 1,
            data_new = dfnew
          )
        }
        lzimp <- rpostpred(nrep,
                           pred_sub_coords,
                           alpha_cur,
                           lKimp,
                           matrix(out$smp_lB[i,], nrow = L))
        yimp <- lsm2gev(lzimp, alpha_cur, theta_cur,
                        lKimp, muimp, sigmaimp, xiimp)
        y_smps[smp_loc, ,] <- yimp
        smp_loc <- smp_loc + 1
      }
    }
    return(list(y = y_smps, pred_sub_coords_list = pred_sub_coords_list))
  }


#' Wrapper for cond_sim_lK and cond_sim_mpar to make conditional draws of latent processes
#' @param out MCMC sample output from mcmc_lnorm_basis
#' @param obs_coords nloc x 2 dimensional matrix of observation location coordinates
#' @param pred_coords n-prediction-location x 2 dimensional matrix of
#' prediction location coordinates at which samples of the latent processes
#' are to be conditionally simulated
#' @param thin_int (optional) thinning interval - if only a subset of MCMC
#' samples should be used (e.g. for memory constraints/testing)
#' @param sim_y (logical) should conditional draws of the stable-mixture
#' response be simulated at \code{pred_coords}?
#' @param npred number of predictive draws to make
#' @export
lnorm_cond_sim_spatial_chi <- function(out,
                                       obs_coords,
                                       pred_coords,
                                       thin_int = 1,
                                       sim_y = TRUE,
                                       npred = 100) {
  pred_sub_coords_list <- list()
  nmcmc <- nrow(out$smp_par)
  n_smp <-
    ifelse((thin_int == 1),
           yes = nmcmc,
           no = floor(nmcmc / thin_int))
  L <- dim(out$smp_lK)[3]
  nrep <- dim(out$smp_lB)[2] / dim(out$smp_lK)[3]
  dfnew <- data.frame(x = rep(1, npred * nrep))

  y_smps <- array(NA_real_, dim = c(n_smp, nrep, npred))

  smp_loc <- 1
  for (i in 1:nmcmc) {
    if (i %% thin_int == 0) {
      pred_sub_coords_list[[smp_loc]] <- pred_sub_coords <-
        pred_coords[sample(1:nrow(pred_coords),
                           size = npred,
                           replace = F),]
      alpha_cur <- out$smp_par[i, "alpha"]
      theta_cur <- out$smp_par[i, "theta"]
      # lK ----------------------------------------------------------------------
      lKimp <-
        cond_sim_lK(out$smp_lK[i, , ],
                    pred_sub_coords,
                    obs_coords,
                    out$smp_par[i, "gvar"],
                    out$smp_par[i, "gscl"])

      # Mu ----------------------------------------------------------------------
      muimp <-
        cond_sim_mpar(
          out$smp_mu_gp[i, ],
          pred_sub_coords,
          obs_coords,
          out$smp_mu_gvar[i],
          out$smp_mu_gscl[i],
          nrep = nrep,
          betas = out$smp_par[i, "mu"],
          formula = ~ 1,
          data_new = dfnew
        )

      # Sigma ----------------------------------------------------------------------
      sigmaimp <-
        cond_sim_mpar(
          out$smp_sigma_gp[i, ],
          pred_sub_coords,
          obs_coords,
          out$smp_sigma_gvar[i],
          out$smp_sigma_gscl[i],
          nrep = nrep,
          betas = out$smp_par[i, "sigma"],
          formula = ~ 1,
          data_new = dfnew,
          log.link = TRUE
        )


      # Xi ----------------------------------------------------------------------
      if (is.null(out$smp_xi_gp)) {
        xiimp <- matrix(out$smp_par[i, "xi"], nrow = nrep, ncol = npred)
      }
      else{
        xiimp <-
          cond_sim_mpar(
            out$smp_xi_gp[i, ],
            pred_sub_coords,
            obs_coords,
            out$smp_xi_gvar[i],
            out$smp_xi_gscl[i],
            nrep = nrep,
            betas = out$smp_par[i, "xi"],
            formula = ~ 1,
            data_new = dfnew
          )
      }
      # Response ----------------------------------------------------------------
      lzimp <-
        rpostpred(nrep,
                  pred_sub_coords,
                  alpha_cur,
                  lKimp,
                  matrix(out$smp_lB[i, ], nrow = L))
      yimp <-
        lsm2gev(lzimp, alpha_cur, theta_cur, lKimp, muimp, sigmaimp, xiimp)
      y_smps[smp_loc, , ] <- yimp
      smp_loc <- smp_loc + 1
    }
  }
  return(list(y = y_smps, pred_sub_coords_list = pred_sub_coords_list))
}

