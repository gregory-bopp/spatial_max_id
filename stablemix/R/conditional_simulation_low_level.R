#' Simulate GEV marginal parameter
#' @param snew observation locations
#' @param gvar GP variance
#' @param gscl GP range
#' @param nrep number of replicates
#' @param betas linear coefficients in, e.g., loc = X beta
#' @param formula R formula describing form of covariates in GP mean
#' @param data_new data.frame containing covariate values at \code{snew}
#' locations
#' @param log.link (logical) should a log link be used for GEV parameter? This
#' should be set to TRUE when parameter is constrained to be positive. For
#' example in the GEV scale case.
#' @export
sim_mpar <- function(snew,
                     gvar,
                     gscl,
                     nrep = NULL,
                     betas = NULL,
                     formula = NULL,
                     data_new = NULL,
                     log.link = FALSE){
  model <- RandomFields::RMexp(var = gvar, scale = gscl)
  RandomFields::RFoptions(spConform = FALSE)
  mparimp_gp <- RandomFields::RFsimulate(
    n = 1,
    model = model,
    x = snew[,1],
    y = snew[,2]
  )
  if(!is.null(formula)){
    nnew <- nrow(snew)
    X <- model.matrix(formula, data = data_new)
    mparimp_mean <- matrix(X %*% betas, nrow = nrep, ncol = nnew, byrow = TRUE)
    mparimp_gp <- matrix(mparimp_gp, nrow = nrep, ncol = nnew, byrow = TRUE) + mparimp_mean
  }
  if(log.link == TRUE){
    mparimp_gp <- exp(mparimp_gp)
  }
  return(mparimp_gp)
}

#' Make conditional simulations of response conditional on low-level processes (log-normal basis)
#' @param out MCMC output object
#' @param pred_coords matrix of prediction coordinates
#' @param thin_int Thinning interval
#' @param sim_y (logical) Should the response be simulated at prediction locations
#' @export
lnorm_cond_sim_low_all <- function(out,
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
      lKimp <- make_kernel(
        pred_coords,
        nbasis = L,
        type = "br",
        gvar = out$smp_par[i, "gvar"],
        gscl = out$smp_par[i, "gscl"]
      )

      # Mu ----------------------------------------------------------------------
      muimp <-  sim_mpar(
        pred_coords,
        gvar = out$smp_mu_gvar[i],
        gscl = out$smp_mu_gscl[i],
        nrep = nrep,
        betas = out$smp_par[i, "mu"],
        formula = ~ 1,
        data_new = dfnew
      )

      # Sigma ----------------------------------------------------------------------
      sigmaimp <-  sim_mpar(
        pred_coords,
        gvar = out$smp_sigma_gvar[i],
        gscl = out$smp_sigma_gscl[i],
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
        xiimp <- sim_mpar(
          pred_coords,
          gvar = out$smp_xi_gvar[i],
          gscl = out$smp_xi_gscl[i],
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

#' Make conditional simulations of response conditional on low-level processes  (density basis)
#' @param out MCMC output object
#' @param pred_coords matrix of prediction coordinates
#' @param knot_coords matrix of knot coordinates for Gaussian density basis
#' @param thin_int Thinning interval
#' @param sim_y (logical) Should the response be simulated at prediction locations
#' @export
fixed_cond_sim_low_all <- function(out,
                                   pred_coords,
                                   knot_coords,
                                   thin_int = 1,
                                   sim_y = TRUE){
  # Inside function ---------------------------------------------------------
  npred <- nrow(pred_coords)
  nmcmc <- nrow(out$smp_par)
  n_smp <- ifelse((thin_int == 1), yes = nmcmc, no = floor(nmcmc/thin_int))
  L <- nrow(knot_coords)
  nrep <- ncol(out$smp_lB)/L
  dfnew <- data.frame(x = rep(1, npred*nrep))

  lK_smps <- array(NA_real_, dim = c(n_smp, npred, L))
  mu_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)
  sigma_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)
  xi_smps <- matrix(NA_real_, nrow = n_smp, ncol = npred)

  if(sim_y){
    y_smps <- array(NA_real_, dim = c(n_smp, nrep, npred))
  }

  smp_loc <- 1
  for(i in 1:nmcmc){
    if (i%%thin_int == 0) {
      alpha_cur <- out$smp_par[i,"alpha"]
      theta_cur <- out$smp_par[i,"theta"]
      # lK ----------------------------------------------------------------------
      lKimp <- make_kernel(pred_coords,
                           nbasis =  L,
                           type = "smith",
                           knot_coord = knot_coords,
                           kern_bw = out$smp_par[i,"tau"])

      # Mu ----------------------------------------------------------------------
      muimp <-  sim_mpar(pred_coords,
                         gvar = out$smp_mu_gvar[i],
                         gscl = out$smp_mu_gscl[i],
                         nrep = nrep,
                         betas = out$smp_par[i,"mu"],
                         formula = ~ 1,
                         data_new = dfnew)

      # Sigma ----------------------------------------------------------------------
      sigmaimp <-  sim_mpar(pred_coords,
                            gvar = out$smp_sigma_gvar[i],
                            gscl = out$smp_sigma_gscl[i],
                            nrep = nrep,
                            betas = out$smp_par[i,"sigma"],
                            formula = ~ 1,
                            data_new = dfnew,
                            log.link = TRUE)

      # Xi ----------------------------------------------------------------------
      if(is.null(out$smp_xi_gp)){
        xiimp <-matrix(out$smp_par[i,"xi"], nrow = nrep, ncol = npred)
      }
      else{
        xiimp <- sim_mpar(pred_coords,
                          gvar = out$smp_xi_gvar[i],
                          gscl = out$smp_xi_gscl[i],
                          nrep = nrep,
                          betas = out$smp_par[i,"xi"],
                          formula = ~ 1,
                          data_new = dfnew)
      }
      # Response ----------------------------------------------------------------
      # n is nrep
      if(sim_y){
        lzimp <- rpostpred(nrep, pred_coords, alpha_cur, lKimp, matrix(out$smp_lB[i,], nrow = L))
        yimp <- lsm2gev(lzimp, alpha_cur, theta_cur, lKimp, muimp, sigmaimp, xiimp)
        y_smps[smp_loc,,] <- yimp
      }
      lK_smps[smp_loc,,] <- lKimp
      mu_smps[smp_loc,] <- muimp[1,]
      sigma_smps[smp_loc,] <- sigmaimp[1,]
      xi_smps[smp_loc,] <- xiimp[1,]
      smp_loc <- smp_loc + 1
    }
  }
  if(sim_y){
    return(list(y = y_smps, lK = lK_smps, mu = mu_smps, sigma = sigma_smps, xi = xi_smps))
  }
  else{
    return(list(lK = lK_smps, mu = mu_smps, sigma = sigma_smps, xi = xi_smps))
  }
}
