#' Metropolis-Hastings algorithm for log-Gaussian process basis stable mixture model
#' @description Main Metropolis-Hastings loop to generate posterior samples for
#' stable-mixture model when using log-Gaussian process kernel function basis.
#' @param n_iter number of MCMC iterations/posterior samples to generate
#' @param obs_coord n-location x 2 matrix of observation location coordinates
#' @param y n-replicate x n-location matrix of observations with GEV margins
#' to which a stable-mixture model will be fit.
#' @param lB n-basis x n-replicate matrix of initial values for log-transformed
#' random kernel scaling factors (assumed to follow a Hougaard distribution).
#' @param lK n-location x n-basis matrix of initial values for normalized
#' (sum to 1) log-Gaussian process prior kernel functions, which are transformed
#' to the log-scale after normalization.
#' @param init named list of initial values for model parameters. Should contain
#' alpha, theta, tau (kernel bw), loc (mean function betas for GEV location),
#' scale (mean function betas for log(GEV scale)), shape.
#' @param tune_var named list of tuning proposal variance parameters containing
#' alpha, theta, tau (kernel bandwidth) loc, scale, shape (each of which should
#' match the beta coefficient vector dimension in e.g. loc = X %*% beta),
#' lB, and gps (which is itslef a list with names gvar, gscl, and gp) for GEV
#' marginal parameter Gaussian processes.
#' @param data data.frame containing covariates with which to model GEV
#' marginal parameter means.
#' @param mar_gp_par Named list containing initial values for GEV marginal
#' parameter Gaussian process hyper parameter initial values. E.g.,
#' list(mu = list(gvar = 0.01, gscl = 1/3),
#' sigma = list(gvar = 0.01, gscl = 1/3),
#' xi =list(gvar = 0.01, gscl = 1/3))
#' @param mar_gp_init (Optional) list of initial values for GEV location,
#' scale, and shape GP residual processes (i.e. with mean zero).
#' @param mar_gp_which vector containing a subset of c("mu","sigma","xi")
#' indicating which GEV marginal parameters should be modeled as spatial
#' Gaussian processes, and which should be held constant in space. Note: if
#' "sigma" is specified, it will be modeled as a log-Gaussian process to
#' satisfy the positivity constraint.
#' @param loc R formula for (GP mean function of, depending) GEV location
#' parameter. Any specified covariates should be present in 'data'
#' @param scale R formula for (GP mean function of, depending) GEV scale
#' parameter. Any specified covariates should be present in 'data'
#' @param shape R formula for (GP mean function of, depending) GEV shape
#' parameter. Any specified covariates should be present in the 'data'
#' @param opt_rt Target acceptance rate for sampler (applies to all parameters)
#' Should be a number between 0 and 1.
#' @param select_lat If only a subset of the latent parameters should be
#' returned (e.g. for the first 5 locations), pass a vector of indices (e.g. 1:5)
#' @param thin_int Thinning interval if there are memory constraints for
#' mcmc sample storage.
#' @param tau_tol (optional) kernel bandwidth tolerance to bound bandwidth
#' away from zero, which might cause kernel degeneracy and numerical instability.
#' @param alpha_tol (optional) tolerance to bound alpha parameter away from 0
#' and 1.
#' @param clusters (optional) list of indices for homogeneous subregions of
#' observation locations (should match row numbers of observation location
#' matrix). This blocking scheme is used to improve the efficiency of the
#' Gaussian process updates.
#' @param infer_general_theta (not currently implemented)
#' @param parallelize (logical) If multiple cores are available, can parallelize
#' over the independent replicates of the process when updating the random
#' basis scaling factors belonging to different replicates (e.g., years)
#' @export
#' @return list containing
#'   smp_par: posterior samples of model parameters
#'   smp_lB: posterior samples of basis scaling factors (on log scale)
#'   smp_lK: posterior samples of normalized log-Gaussian process basis
#'           functions (on log scale)
#'   smp_mu_gp: location parameter (mean-zero) residual gaussian process
#'   smp_sigma_gp: log(scale parameter) (mean-zero) residual gaussian process
#'   smp_xi_gp: shape parameter (mean-zero) residual gaussian process
#'   smp_mu_gvar: GEV location parameter GP variance parameter posterior samples
#'   smp_sigma_gvar: GEV log(scale parameter) GP variance parameter posterior samples
#'   smp_xi_gvar: GEV shape parameter GP variance parameter posterior samples
#'   smp_mu_gscl: GEV location parameter GP location parameter posterior samples
#'   smp_sigma_gvar: GEV log(scale parameter) GP location parameter posterior samples
#"   smp_xi_gvar: GEV shape parameter GP location parameter posterior samples
#'   acpt_rts: matrix of acceptance rates for sampled parameters
#'   acpt_rts_loc: matrix of acceptance rates for GEV location parameter GPs
#'   acpt_rts_scale: matrix of acceptance rates for GEV log(scale) parameter GPs
#'   acpt_rts_shape: matrix of acceptance rates for GEV shape parameter GPs
#'   acpt_rts_lB: matrix of acceptance rates for random basis scaling factors
#'   acpt_rts_lK: matrix of acceptance rates for random basis function terms
#'   prop_var: proposal variances for model parameter updates
#'   loc_prop_var: proposal variances for GEV location parameter
#'   scale_prop_var: proposal variances for GEV log(scale) parameter
#'   shape_prop_var: proposal variances for GEV shape parameter
#'   clusters: observation location cluster block definitions for efficient
#'   mcmc sampling
mcmc_lnorm_basis <-
  function(n_iter,
           obs_coord,
           y,
           lB,
           lK,
           init,
           tune_var,
           data,
           mar_gp_par,
           mar_gp_init = NULL,
           mar_gp_which = NULL,
           loc = ~ 1,
           scale = ~ 1,
           shape = ~ 1,
           opt_rt = 0.40,
           select_lat = NULL,
           thin_int = 1,
           alpha_tol = 0.07,
           clusters = NULL,
           infer_general_theta = FALSE,
           parallelize = FALSE) {

  if(!is.matrix(y)){
    stop("y must either be a matrix with independent
         replicates in each row. Use drop = FALSE if using
         only a single replicate.")
  }
  if(init["theta"] ==0){
    fix_theta0 = TRUE
  }
  else{
    fix_theta0 = FALSE
  }

  Xloc <- model.matrix(loc, data = data)
  Xscale <- model.matrix(scale, data = data)
  Xshape <- model.matrix(shape, data = data)
  locnames <- make_vnames(colnames(Xloc), "mu")
  scalenames <- make_vnames(colnames(Xscale), "sigma")
  shapenames <- make_vnames(colnames(Xshape), "xi")

  ploc <- ncol(Xloc)
  pscale <- ncol(Xscale)
  pshape <- ncol(Xshape)

  nrep <- nrow(y)
  nobs <- length(y)
  nbasis <- length(lB)
  nbasis_per_rep <- nbasis/nrep
  D <- calc_dist(obs_coord, obs_coord)
  nloc <- nrow(D)
  obs_rng <- max(D)

  win_len <- min(2*nobs, n_iter, 50)            # Window for acceptance rates and prop vars
  parnames <- c("alpha", "theta","gvar","gscl")

  if(is.null(clusters)){
    km <- kmeans(obs_coord, floor(nrow(obs_coord)/20))
    clusters <- list(mu = km$cluster, sigma= km$cluster, xi = km$cluster, lK = km$cluster)
  }
  nclus <- lapply(clusters, function(x)length(unique(x)))



  # Initialize parameters and latent variables
  alpha_cur <- init[["alpha"]];  theta_cur <- init[["theta"]];
  gvar_cur <- init[["gvar"]]; gscl_cur <- init[["gscl"]]
  locpar_cur <- init[["loc"]]; scalepar_cur <- init[["scale"]]; shapepar_cur <- init[["shape"]]

  if(is.null(mar_gp_init)){
    mu_cur <- matrix(Xloc%*%locpar_cur, nrow = nrep, ncol = nloc, byrow = TRUE)
    sigma_cur <- matrix(exp(Xscale%*%scalepar_cur), nrow = nrep, ncol = nloc, byrow = TRUE)
    xi_cur <- matrix(Xshape%*%shapepar_cur, nrow = nrep, ncol = nloc, byrow = TRUE)

    # Marginal parameter Gaussian processes
    mu_gp_cur <- mu_gp_prop <- rep(0, nloc)
    sigma_gp_cur <- sigma_gp_prop <- rep(0, nloc)
    xi_gp_cur <- xi_gp_prop <- rep(0, nloc)
  }
  else{
    mu_cur <-calc_mpar(Xloc, locpar_cur, mar_gp_init[["mu"]], nrep, nloc)
    sigma_cur <-calc_mpar(Xscale, scalepar_cur, mar_gp_init[["sigma"]], nrep, nloc, log.link = TRUE)
    xi_cur <-calc_mpar(Xshape, shapepar_cur, mar_gp_init[["xi"]], nrep, nloc)

    mu_gp_cur <- mu_gp_prop  <- mar_gp_init[["mu"]]
    sigma_gp_cur <- sigma_gp_prop <- mar_gp_init[["sigma"]]
    xi_gp_cur <- xi_gp_prop <- mar_gp_init[["xi"]]
  }

  lB_prop <- lB_cur <- lB
  lK_prop <- lK_cur <- lK

  # Store parameter and latent variable samples
  n_smp <- ifelse((thin_int ==1), yes = n_iter,
                  no = floor(n_iter/thin_int))      # Number of samples retained after thinning
  smp_alpha <- smp_theta <- smp_gvar <- smp_gscl <- rep(NA_real_, n_smp)
  smp_locpar <- matrix(NA_real_, nrow = n_smp, ncol = ploc)
  smp_scalepar <- matrix(NA_real_, nrow = n_smp, ncol = pscale)
  smp_shapepar <- matrix(NA_real_, nrow = n_smp, ncol = pshape)

  if(!is.null(select_lat)){
    smp_lB <- matrix(NA_real_, nrow = n_smp, ncol = length(select_lat))
    smp_lK <- array(NA_real_, dim = c(n_smp, nloc, length(select_lat)))
  }
  else{
    smp_lB <- matrix(NA_real_, nrow = n_smp, ncol = nbasis)
    smp_lK <- array(NA_real_, dim = c(n_smp, nloc, nbasis_per_rep))
  }

  # Automatic Tuning of Proposal variance
  c0 <- 10;  c1 <- 0.8;  tune_k <- 3

  # Proposal Variances
  prop_var <- matrix(NA_real_, nrow = win_len, ncol = length(parnames),
                     dimnames = list(NULL, parnames))
  prop_var[1,] <- c(tune_var[["alpha"]], tune_var[["theta"]],  tune_var[["gvar"]], tune_var[["gscl"]])
  loc_prop_var <- matrix(NA_real_, nrow = win_len, ncol = ploc)
  loc_prop_var[1,] <- tune_var[["loc"]]
  scale_prop_var <- matrix(NA_real_, nrow = win_len, ncol = pscale)
  scale_prop_var[1,] <- tune_var[["scale"]]
  shape_prop_var <- matrix(NA_real_, nrow = win_len, ncol = pshape)
  shape_prop_var[1,] <- tune_var[["shape"]]
  mgvar_prop_var <- mgscl_prop_var <- matrix(NA_real_, win_len, ncol = 3,
                                             dimnames = list(NULL, c("mu","sigma","xi")))

  lB_prop_var <- matrix(NA_real_, nrow = win_len, ncol =nbasis)
  lB_prop_var[1,] <- tune_var[["lB"]]
  lK_prop_var <- matrix(NA_real_, nrow = win_len, ncol = nclus[["lK"]])
  lK_prop_var[1,] <- tune_var[["lK"]]

  # Acceptance matrix for parameters
  acpt <- matrix(NA_real_, nrow = win_len, ncol = length(parnames), dimnames = list(NULL, parnames))
  acpt_loc <- matrix(NA_real_, nrow = win_len, ncol = ploc)
  acpt_scale <- matrix(NA_real_, nrow = win_len, ncol = pscale)
  acpt_shape <- matrix(NA_real_, nrow = win_len, ncol = pshape)
  acpt_mgvar <- acpt_mgscl <- matrix(NA_real_, nrow = win_len, ncol = 3,
                                     dimnames = list(NULL, c("mu","sigma","xi")))
  acpt_lB <- matrix(NA_real_, nrow = win_len, ncol = nbasis)
  acpt_lK <- matrix(NA_real_, nrow = win_len, ncol = nclus[["lK"]])

  # Make Basis functions and A terms
  lA_cur <- mklA(lK_cur, lB_cur, alpha_cur)
  qfuns_cur <- create_qlstabmix(alpha_cur, theta_cur, lK_cur)

  if("mu"%in%mar_gp_which){
    mu_gvar_cur <- mar_gp_par[["mu"]][["gvar"]]
    mu_gscl_cur <- mar_gp_par[["mu"]][["gscl"]]
    smp_mu_gp <- matrix(NA_real_, nrow = n_smp, ncol = nloc)
    smp_mu_gvar <- smp_mu_gscl <- rep(NA_real_, n_smp)
    mu_gp_prop_var  <- matrix(NA_real_, nrow = win_len, ncol =nloc)
    mu_gp_prop_var[1,] <-  tune_var[["gps"]][["gp"]]
    mgvar_prop_var[1,"mu"] <- tune_var[["gps"]][["gvar"]]
    mgscl_prop_var[1,"mu"] <- tune_var[["gps"]][["gscl"]]
    acpt_gp_mu <-  matrix(NA_real_, nrow = win_len, ncol = nclus[["mu"]])
  }
  if("sigma"%in%mar_gp_which){
    sigma_gvar_cur <- mar_gp_par[["sigma"]][["gvar"]]
    sigma_gscl_cur <- mar_gp_par[["sigma"]][["gscl"]]
    smp_sigma_gp <- matrix(NA_real_, nrow = n_smp, ncol = nloc)
    smp_sigma_gvar <- smp_sigma_gscl <- rep(NA_real_, n_smp)
    sigma_gp_prop_var <- matrix(NA_real_, nrow = win_len, ncol =nloc)
    sigma_gp_prop_var[1,] <-  tune_var[["gps"]][["gp"]]
    mgvar_prop_var[1,"sigma"] <- tune_var[["gps"]][["gvar"]]
    mgscl_prop_var[1,"sigma"] <- tune_var[["gps"]][["gscl"]]
    acpt_gp_sigma <- matrix(NA_real_, nrow = win_len, ncol = nclus[["sigma"]])
  }
  if("xi"%in%mar_gp_which){
    xi_gvar_cur <- mar_gp_par[["xi"]][["gvar"]]
    xi_gscl_cur <- mar_gp_par[["xi"]][["gscl"]]
    smp_xi_gp <- matrix(NA_real_, nrow = n_smp, ncol = nloc)
    smp_xi_gvar <- smp_xi_gscl <- rep(NA_real_, n_smp)
    xi_gp_prop_var <- matrix(NA_real_, nrow = win_len, ncol =nloc)
    xi_gp_prop_var[1,] <-  tune_var[["gps"]][["gp"]]
    mgvar_prop_var[1,"xi"] <- tune_var[["gps"]][["gvar"]]
    mgscl_prop_var[1,"xi"] <- tune_var[["gps"]][["gscl"]]
    acpt_gp_xi <- matrix(NA_real_, nrow = win_len, ncol = nclus[["xi"]])
  }
  # Check for missing values (draw from posterior predictive)
  miss.ind <- is.na(y)
  n.miss <- sum(miss.ind)
  if(n.miss ==0){
    lz_cur <- gev2lsm(y, mu_cur, sigma_cur, xi_cur, qfuns_cur)
  }

  smp_loc <- 1                                      # index of retained samples
  for(i in 1:n_iter){
    # Draws from posterior predictive distribution
    if(n.miss>0){
      lzp <- rpostpred(nrep, obs_coord, alpha_cur, lK_cur, lB_cur)
      yp <- lsm2gev(lzp, alpha_cur, theta_cur, lK_cur, mu_cur, sigma_cur, xi_cur)
      y[miss.ind] <- yp[miss.ind]
      lz_cur <-gev2lsm(y, mu_cur, sigma_cur, xi_cur, qfuns_cur)
    }

    pos <- (i-1)%%win_len + 1
    nxt_pos <- i%%win_len + 1
    gamma1 <- c0 / (i + tune_k)^(c1)                # Automatic tuning constant
    # Update alpha
    alpha_prop <- rnorm_prop(alpha_cur, tune_var = prop_var[pos, "alpha"])
    if(alpha_prop > alpha_tol){
      qfuns_prop <- create_qlstabmix(alpha_prop, theta_cur, lK_cur)
      lz_prop <- gev2lsm(y, mu_cur, sigma_cur, xi_cur, qfuns_prop)
      lA_prop <- mklA(lK_cur, lB_cur, alpha_prop)
      if((!any(is.na(lz_prop)))&
         (log(runif(1)) <
         lpst_alpha(alpha_prop, lz_prop, y, lA_prop, lB_cur, lK_cur,
                    theta_cur, mu_cur, sigma_cur, xi_cur, tol = alpha_tol)-
         lpst_alpha(alpha_cur, lz_cur, y, lA_cur, lB_cur, lK_cur,
                    theta_cur, mu_cur, sigma_cur, xi_cur, tol = alpha_tol)
      )){
        alpha_cur <- alpha_prop
        qfuns_cur <- qfuns_prop
        lz_cur <- lz_prop
        lA_cur <- lA_prop
        acpt[pos, "alpha"] <- 1
      }
      else{
        acpt[pos, "alpha"] <- 0
      }
    }
    else{
      acpt[pos, "alpha"] <- 0
    }
    prop_var[nxt_pos, "alpha"] <- update_var(prop_var[pos, "alpha"], acpt_rt = mean(acpt[,"alpha"], na.rm = TRUE),
                                             opt_rt = opt_rt, gamma1 = gamma1)

    # Update theta
    if(!fix_theta0){
      theta_prop <- exp(rnorm_prop(log(theta_cur), tune_var = prop_var[pos, "theta"]))
      if(theta_prop > 0){
        qfuns_prop <- create_qlstabmix(alpha_cur, theta_prop, lK_cur)
        lz_prop <- gev2lsm(y, mu_cur, sigma_cur, xi_cur, qfuns_prop)
        if((!any(is.na(lz_prop)))&
           (log(runif(1)) <
            lpst_theta(theta_prop, lz_prop, y, lB_cur, lA_cur, lK_cur, alpha_cur, mu_cur, sigma_cur, xi_cur) -
            lpst_theta(theta_cur, lz_cur, y, lB_cur, lA_cur, lK_cur, alpha_cur, mu_cur, sigma_cur, xi_cur)
           )){
          theta_cur <- theta_prop
          qfuns_cur <- qfuns_prop
          lz_cur <- lz_prop
          acpt[pos, "theta"] <- 1
        }
        else{
          acpt[pos, "theta"] <- 0
        }
      }
      else{
        acpt[pos, "theta"] <- 0
      }
      prop_var[nxt_pos, "theta"] <- update_var(prop_var[pos, "theta"], acpt_rt = mean(acpt[,"theta"], na.rm = TRUE),
                                               opt_rt = opt_rt, gamma1 = gamma1)
    }

    # Update gvar and gscl
    gvar_prop <- rnorm_prop(gvar_cur, tune_var = prop_var[pos, "gvar"])
    if(log(runif(1)) <
       lpst_gvar(gvar_prop, obs_coord, lK_cur, gscl_cur) -
       lpst_gvar(gvar_cur, obs_coord, lK_cur, gscl_cur)
    ){
      gvar_cur <- gvar_prop
      acpt[pos, "gvar"] <- 1
    }
    else{
      acpt[pos, "gvar"] <- 0
    }
    prop_var[nxt_pos, "gvar"] <- update_var(prop_var[pos, "gvar"], acpt_rt = mean(acpt[,"gvar"], na.rm = TRUE),
                                          opt_rt = opt_rt, gamma1 = gamma1)

    # Update gscl
    gscl_prop <- rnorm_prop(gscl_cur, tune_var = prop_var[pos, "gscl"])
    if(gscl_prop < obs_rng){
      if(log(runif(1)) <
         lpst_gscl(gscl_prop, obs_coord, lK_cur, gvar_cur, obs_rng) -
         lpst_gscl(gscl_cur, obs_coord, lK_cur, gvar_cur, obs_rng)
      ){
        gscl_cur <- gscl_prop
        acpt[pos, "gscl"] <- 1
      }
      else{
        acpt[pos, "gscl"] <- 0
      }
    }
    else{
      acpt[pos, "gscl"] <- 0
    }
    prop_var[nxt_pos, "gscl"] <- update_var(prop_var[pos, "gscl"], acpt_rt = mean(acpt[,"gscl"], na.rm = TRUE),
                                            opt_rt = opt_rt, gamma1 = gamma1)

    # Update mu
    # betas
    locpar_prop <- locpar_cur
    for(j in 1:ploc){
      locpar_prop[j] <- rnorm_prop(locpar_cur[j], tune_var = loc_prop_var[pos, j])
      mu_prop <- calc_mpar(Xloc, locpar_prop, mu_gp_cur, nrep, nloc)
      lz_prop <- gev2lsm(y, mu_prop, sigma_cur, xi_cur, qfuns_cur)
      if((!any(is.na(lz_prop)))&
         (log(runif(1)) <
         lpst_mpar_mean(locpar_prop[j], lz_prop, y, lA_cur, lK_cur, alpha_cur, theta_cur, mu_prop, sigma_cur, xi_cur) -
         lpst_mpar_mean(locpar_cur[j], lz_cur, y, lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_cur)
      )){
        mu_cur <- mu_prop
        locpar_cur[j] <- locpar_prop[j]
        lz_cur <- lz_prop
        acpt_loc[pos, j] <- 1
      }
      else{
        locpar_prop[j] <- locpar_cur[j]
        acpt_loc[pos, j] <- 0
      }
      loc_prop_var[nxt_pos, j] <- update_var(loc_prop_var[pos, j], acpt_rt = mean(acpt_loc[,j], na.rm = TRUE),
                                             opt_rt = opt_rt, gamma1 = gamma1)
    }
    if("mu"%in%mar_gp_which){
      # gaussian process error
      cmod_cur <- RandomFields::RMexp(var = mu_gvar_cur, scale = mu_gscl_cur)
      ll_cur <- lpst_mpar_gp(mu_gp_cur, mu_gvar_cur, mu_gscl_cur,obs_coord, lz_cur, y,
                             lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_cur)
      Cmat <- RandomFields::RFcovmatrix(cmod_cur, obs_coord)
      for(j in 1:nclus[["mu"]]){
        mu_gp_prop <- mu_gp_cur
        sub <- (clusters[["mu"]] ==j)
        if(any(sub)){
          Csub <- Cmat[sub, sub]
          Lt <- chol(Csub)
          white <- mu_gp_cur[sub]%*%solve(Lt)
          white_shifted <- white + rnorm(sum(sub), mean = 0, sd = sqrt(mu_gp_prop_var[pos, j]))
          mu_gp_prop[sub] <- white_shifted%*%Lt
          mu_prop <- calc_mpar(Xloc, locpar_cur, mu_gp_prop, nrep, nloc)
          lz_prop <- gev2lsm(y, mu_prop, sigma_cur, xi_cur, qfuns_cur)
          ll_prop <- lpst_mpar_gp(mu_gp_prop, mu_gvar_cur, mu_gscl_cur,obs_coord,lz_prop, y,
                                  lA_cur, lK_cur, alpha_cur, theta_cur, mu_prop, sigma_cur, xi_cur)
          if((!any(is.na(lz_prop)))&
             (log(runif(1)) <ll_prop - ll_cur)){
            mu_cur <- mu_prop
            mu_gp_cur <- mu_gp_prop
            lz_cur <- lz_prop
            acpt_gp_mu[pos, j] <- 1
            ll_cur <- ll_prop
          }
          else{
            acpt_gp_mu[pos, j] <- 0
          }
          mu_gp_prop_var[nxt_pos, j] <- update_var(mu_gp_prop_var[pos, j], acpt_rt = mean(acpt_gp_mu[,j], na.rm = TRUE),
                                                   opt_rt = opt_rt, gamma1 = gamma1)
        }
      }
      # gvar
      mgvar_prop <- rnorm_prop(mu_gvar_cur, tune_var = mgvar_prop_var[pos,"mu"])
      if(mgvar_prop > 0){
        if(log(runif(1)) <
           lpst_gp_cov_par(mu_gp_cur, obs_coord, mgvar_prop, mu_gscl_cur) -
           lpst_gp_cov_par(mu_gp_cur, obs_coord, mu_gvar_cur, mu_gscl_cur)
        ){
          mu_gvar_cur <- mgvar_prop
          acpt_mgvar[pos,"mu"] <- 1
        }
        else{
          acpt_mgvar[pos,"mu"]  <- 0
        }
      }
      else{
        acpt_mgvar[pos,"mu"]<- 0
      }
      mgvar_prop_var[nxt_pos,"mu"] <- update_var(mgvar_prop_var[pos,"mu"], acpt_rt = mean(acpt_mgvar[,"mu"], na.rm = TRUE),
                                                 opt_rt = opt_rt, gamma1 = gamma1)
      # gscl
      mgscl_prop <- rnorm_prop(mu_gscl_cur, tune_var = mgscl_prop_var[pos,"mu"])
      if((mgscl_prop > 0)&(mgscl_prop<obs_rng)){
        if(log(runif(1)) <
           lpst_gp_cov_par(mu_gp_cur, obs_coord, mu_gvar_cur, mgscl_prop, rng_max = obs_rng) -
           lpst_gp_cov_par(mu_gp_cur, obs_coord, mu_gvar_cur, mu_gscl_cur, rng_max = obs_rng)
        ){
          mu_gscl_cur <- mgscl_prop
          acpt_mgscl[pos,"mu"] <- 1
        }
        else{
          acpt_mgscl[pos,"mu"]  <- 0
        }
      }
      else{
        acpt_mgscl[pos,"mu"] <- 0
      }
      mgscl_prop_var[nxt_pos,"mu"] <- update_var(mgscl_prop_var[pos,"mu"], acpt_rt = mean(acpt_mgscl[,"mu"], na.rm = TRUE),
                                                 opt_rt = opt_rt, gamma1 = gamma1)

    }

    # Update sigma
    # betas
    scalepar_prop <- scalepar_cur
    for(j in 1:pscale){
      scalepar_prop[j] <- rnorm_prop(scalepar_cur[j], tune_var = scale_prop_var[pos, j])
      sigma_prop <- calc_mpar(Xscale, scalepar_prop, sigma_gp_cur, nrep, nloc, log.link = TRUE)
      if(all(sigma_prop > 0)){
        lz_prop <- gev2lsm(y, mu_cur, sigma_prop, xi_cur, qfuns_cur)
        if((!any(is.na(lz_prop)))&
           (log(runif(1)) <
           lpst_mpar_mean(scalepar_prop[j], lz_prop, y, lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_prop, xi_cur) -
           lpst_mpar_mean(scalepar_cur[j], lz_cur, y, lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_cur)
        )){
          sigma_cur <- sigma_prop
          scalepar_cur[j] <- scalepar_prop[j]
          lz_cur <- lz_prop
          acpt_scale[pos, j] <- 1
        }
        else{
          scalepar_prop[j] <- scalepar_cur[j]
          acpt_scale[pos, j] <- 0
        }
      }
      else{
        scalepar_prop[j] <- scalepar_cur[j]
        acpt_scale[pos, j] <- 0
      }
      scale_prop_var[nxt_pos, j] <- update_var(scale_prop_var[pos, j], acpt_rt = mean(acpt_scale[,j], na.rm = TRUE),
                                               opt_rt = opt_rt, gamma1 = gamma1)
    }
    if("sigma"%in%mar_gp_which){
      # Gaussian Process
      cmod_cur <- RandomFields::RMexp(var = sigma_gvar_cur, scale = sigma_gscl_cur)
      ll_cur <- lpst_mpar_gp(sigma_gp_cur, sigma_gvar_cur, sigma_gscl_cur,obs_coord,lz_cur, y,
                             lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_cur)
      Cmat <- RandomFields::RFcovmatrix(cmod_cur, obs_coord)
      for(j in 1:nclus[["sigma"]]){
        sigma_gp_prop <- sigma_gp_cur
        sub <- (clusters[["sigma"]] ==j)
        if(any(sub)){
          Csub <- Cmat[sub, sub]
          Lt <- chol(Csub)
          white <- sigma_gp_cur[sub]%*%solve(Lt)
          white_shifted <- white + rnorm(sum(sub), mean = 0, sd = sqrt(sigma_gp_prop_var[pos, j]))
          sigma_gp_prop[sub] <- white_shifted%*%Lt
          sigma_prop <- calc_mpar(Xscale, scalepar_cur, sigma_gp_prop, nrep, nloc, log.link = TRUE)
          lz_prop <- gev2lsm(y, mu_cur, sigma_prop, xi_cur, qfuns_cur)
          ll_prop <- lpst_mpar_gp(sigma_gp_prop, sigma_gvar_cur, sigma_gscl_cur,obs_coord,lz_prop, y,
                                  lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_prop, xi_cur)
          if((!any(is.na(lz_prop)))&
             (log(runif(1)) <ll_prop - ll_cur)){
            sigma_cur <- sigma_prop
            sigma_gp_cur <- sigma_gp_prop
            lz_cur <- lz_prop
            acpt_gp_sigma[pos, j] <- 1
            ll_cur <- ll_prop
          }
          else{
            acpt_gp_sigma[pos, j] <- 0
          }
          sigma_gp_prop_var[nxt_pos, j] <- update_var(sigma_gp_prop_var[pos, j], acpt_rt = mean(acpt_gp_sigma[,j], na.rm = TRUE),
                                                      opt_rt = opt_rt, gamma1 = gamma1)
        }
      }
      # gvar
      mgvar_prop <- rnorm_prop(sigma_gvar_cur, tune_var = mgvar_prop_var[pos,"sigma"])
      if(mgvar_prop > 0){
        if(log(runif(1)) <
           lpst_gp_cov_par(sigma_gp_cur, obs_coord, mgvar_prop, sigma_gscl_cur) -
           lpst_gp_cov_par(sigma_gp_cur, obs_coord, sigma_gvar_cur, sigma_gscl_cur)
        ){
          sigma_gvar_cur <- mgvar_prop
          acpt_mgvar[pos,"sigma"] <- 1
        }
        else{
          acpt_mgvar[pos,"sigma"] <- 0
        }
      }
      else{
        acpt_mgvar[pos,"sigma"]  <- 0
      }
      mgvar_prop_var[nxt_pos,"sigma"] <- update_var(mgvar_prop_var[pos,"sigma"], acpt_rt = mean(acpt_mgvar[,"sigma"], na.rm = TRUE),
                                                    opt_rt = opt_rt, gamma1 = gamma1)
      # gscl
      mgscl_prop <- rnorm_prop(sigma_gscl_cur, tune_var = mgscl_prop_var[pos,"sigma"])
      if((mgscl_prop > 0)&(mgscl_prop<obs_rng)){
        if(log(runif(1)) <
           lpst_gp_cov_par(sigma_gp_cur, obs_coord, sigma_gvar_cur, mgscl_prop, rng_max = obs_rng) -
           lpst_gp_cov_par(sigma_gp_cur, obs_coord, sigma_gvar_cur, sigma_gscl_cur, rng_max = obs_rng)
        ){
          sigma_gscl_cur <- mgscl_prop
          acpt_mgscl[pos,"sigma"] <- 1
        }
        else{
          acpt_mgscl[pos,"sigma"]  <- 0
        }
      }
      else{
        acpt_mgscl[pos,"sigma"]  <- 0
      }
      mgscl_prop_var[nxt_pos,"sigma"] <- update_var(mgscl_prop_var[pos,"sigma"], acpt_rt = mean(acpt_mgscl[,"sigma"], na.rm = TRUE),
                                                    opt_rt = opt_rt, gamma1 = gamma1)

    }
    # Update xi
    # betas
    shapepar_prop <- shapepar_cur
    for(j in 1:pshape){
      shapepar_prop[j] <- rnorm_prop(shapepar_cur[j], tune_var = shape_prop_var[pos, j])
      xi_prop <- calc_mpar(Xshape, shapepar_prop, xi_gp_cur, nrep, nloc)
      lz_prop <- gev2lsm(y, mu_cur, sigma_cur, xi_prop, qfuns_cur)
      if((!any(is.na(lz_prop)))&
         (log(runif(1)) <
         lpst_mpar_mean(shapepar_prop[j], lz_prop, y, lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_prop, is.xi = TRUE) -
         lpst_mpar_mean(shapepar_cur[j], lz_cur, y, lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_cur, is.xi = TRUE)
      )){
        xi_cur <- xi_prop
        shapepar_cur[j] <- shapepar_prop[j]
        lz_cur <- lz_prop
        acpt_shape[pos, j] <- 1
      }
      else{
        shapepar_prop[j] <- shapepar_cur[j]
        acpt_shape[pos, j] <- 0
      }
      shape_prop_var[nxt_pos, j] <- update_var(shape_prop_var[pos, j], acpt_rt = mean(acpt_shape[,j], na.rm = TRUE),
                                               opt_rt = opt_rt, gamma1 = gamma1)
    }
    if("xi"%in%mar_gp_which){
      # gaussian process error
      cmod_cur <- RandomFields::RMexp(var = xi_gvar_cur, scale = xi_gscl_cur)
      ll_cur <- lpst_mpar_gp(xi_gp_cur, xi_gvar_cur, xi_gscl_cur,obs_coord, lz_cur, y,
                             lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_cur)
      Cmat <- RandomFields::RFcovmatrix(cmod_cur, obs_coord)
      for(j in 1:nclus[["xi"]]){
        xi_gp_prop <- xi_gp_cur
        sub <- (clusters[["xi"]] ==j)
        if(any(sub)){
          Csub <- Cmat[sub, sub]
          Lt <- chol(Csub)
          white <- xi_gp_cur[sub]%*%solve(Lt)
          white_shifted <- white + rnorm(sum(sub), mean = 0, sd = sqrt(xi_gp_prop_var[pos, j]))
          xi_gp_prop[sub] <- white_shifted%*%Lt
          xi_prop <- calc_mpar(Xshape, shapepar_cur, xi_gp_prop, nrep, nloc)
          lz_prop <- gev2lsm(y, mu_cur, sigma_cur, xi_prop, qfuns_cur)
          ll_prop <- lpst_mpar_gp(xi_gp_prop, xi_gvar_cur, xi_gscl_cur,obs_coord,lz_prop, y,
                                  lA_cur, lK_cur, alpha_cur, theta_cur, mu_cur, sigma_cur, xi_prop)
          if((!any(is.na(lz_prop)))&
             (log(runif(1)) <ll_prop - ll_cur)){
            xi_cur <- xi_prop
            xi_gp_cur <- xi_gp_prop
            lz_cur <- lz_prop
            acpt_gp_xi[pos, j] <- 1
            ll_cur <- ll_prop
          }
          else{
            acpt_gp_xi[pos, j] <- 0
          }
          xi_gp_prop_var[nxt_pos, j] <- update_var(xi_gp_prop_var[pos, j], acpt_rt = mean(acpt_gp_xi[,j], na.rm = TRUE),
                                                   opt_rt = opt_rt, gamma1 = gamma1)
        }
      }
      # gvar
      mgvar_prop <- rnorm_prop(xi_gvar_cur, tune_var = mgvar_prop_var[pos,"xi"])
      if(mgvar_prop > 0){
        if(log(runif(1)) <
           lpst_gp_cov_par(xi_gp_cur, obs_coord, mgvar_prop, xi_gscl_cur) -
           lpst_gp_cov_par(xi_gp_cur, obs_coord, xi_gvar_cur, xi_gscl_cur)
        ){
          xi_gvar_cur <- mgvar_prop
          acpt_mgvar[pos,"xi"] <- 1
        }
        else{
          acpt_mgvar[pos,"xi"]  <- 0
        }
      }
      else{
        acpt_mgvar[pos,"xi"]  <- 0
      }
      mgvar_prop_var[nxt_pos,"xi"] <- update_var(mgvar_prop_var[pos,"xi"], acpt_rt = mean(acpt_mgvar[,"xi"], na.rm = TRUE),
                                                 opt_rt = opt_rt, gamma1 = gamma1)
      # gscl
      mgscl_prop <- rnorm_prop(xi_gscl_cur, tune_var = mgscl_prop_var[pos,"xi"])
      if((mgscl_prop > 0)&(mgscl_prop<obs_rng)){
        if(log(runif(1)) <
           lpst_gp_cov_par(xi_gp_cur, obs_coord, xi_gvar_cur, mgscl_prop, rng_max = obs_rng) -
           lpst_gp_cov_par(xi_gp_cur, obs_coord, xi_gvar_cur, xi_gscl_cur, rng_max = obs_rng)
        ){
          xi_gscl_cur <- mgscl_prop
          acpt_mgscl[pos,"xi"] <- 1
        }
        else{
          acpt_mgscl[pos,"xi"] <- 0
        }
      }
      else{
        acpt_mgscl[pos,"xi"] <- 0
      }
      mgscl_prop_var[nxt_pos,"xi"] <- update_var(mgscl_prop_var[pos,"xi"], acpt_rt = mean(acpt_mgscl[,"xi"], na.rm = TRUE),
                                                 opt_rt = opt_rt, gamma1 = gamma1)
    }

    # Update lK
    cmod_cur <- RandomFields::RMexp(var = gvar_cur, scale = gscl_cur)
    ll_cur <- lpst_lK(lK_cur,obs_coord,lz_cur,y, lA_cur, alpha_cur,
                      theta_cur, mu_cur, sigma_cur, xi_cur, gvar_cur, gscl_cur)
    qfuns_prop <- qfuns_cur
    Cmat <- RandomFields::RFcovmatrix(cmod_cur, obs_coord)
    for(j in 1:nclus[["lK"]]){
      lK_prop <- lK_cur
      sub <- (clusters[["lK"]] ==j)
      if(any(sub)){
        subids <- which(sub)
        Csub <- Cmat[sub, sub]
        lK_prop[sub,] <- chol_prop_lK(lK_cur[sub,, drop = FALSE], Csub, lK_prop_var[pos, j], nbasis_per_rep)
        lA_prop <- mklA(lK_prop, lB_cur, alpha_cur)
        for(k in 1:sum(sub)){
          qfuns_prop[subids[k]] <- unlist(create_qlstabmix(alpha_cur, theta_cur, lK_prop[subids[k],, drop = FALSE]))
          lz_prop[,subids[k]] <-  gev2lsm(y[,subids[k],drop = FALSE], mu_cur[,subids[k],drop = FALSE],
                                    sigma_cur[,subids[k],drop = FALSE], xi_cur[,subids[k],drop = FALSE], qfuns_prop[subids[k]])
        }
        ll_prop <- lpst_lK(lK_prop, obs_coord,lz_prop, y, lA_prop, alpha_cur,
                           theta_cur, mu_cur, sigma_cur, xi_cur, gvar_cur, gscl_cur)
        if((!any(is.na(lz_prop)))& (log(runif(1)) < ll_prop -   ll_cur) ){
          lK_cur <- lK_prop
          lA_cur <- lA_prop
          ll_cur <- ll_prop
          qfuns_cur <- qfuns_prop
          acpt_lK[pos, j] <- 1
          lz_cur[,sub] <- lz_prop[,sub]
        }
        else{
          lz_prop[,sub] <- lz_cur[,sub]
          qfuns_prop <- qfuns_cur
          acpt_lK[pos, j] <- 0
        }
        lK_prop_var[nxt_pos, j] <- update_var(lK_prop_var[pos, j],
                                              acpt_rt = mean(acpt_lK[,j], na.rm = TRUE),
                                              opt_rt = opt_rt, gamma1 = gamma1)
      }
    }

    # Update B
    if(parallelize){
      n.core = parallel::detectCores()
      doMC::registerDoMC(cores = n.core)
      lBobj <- foreach::foreach(j = 1:nrep) %dopar% {
        bgrp_st <- (j-1)*nbasis_per_rep + 1
        bgrp_end <- j*nbasis_per_rep
        update_lB(lB_cur[,j, drop = FALSE], lz_cur[j,, drop = FALSE],y[j,, drop = FALSE], lK_cur, alpha_cur,
                  theta_cur, mu_cur[j,,drop = FALSE], sigma_cur[j,,drop = FALSE], xi_cur[j,,drop = FALSE], pos, next_pos, acpt_lB[,bgrp_st:bgrp_end],
                  lB_prop_var[,bgrp_st:bgrp_end])
      }
    }
    else{
      lBobj <- foreach::foreach(j = 1:nrep) %do% {
        bgrp_st <- (j-1)*nbasis_per_rep + 1
        bgrp_end <- j*nbasis_per_rep
        update_lB(lB_cur[,j, drop = FALSE], lz_cur[j,, drop = FALSE],y[j,, drop = FALSE], lK_cur, alpha_cur,
                  theta_cur, mu_cur[j,,drop = FALSE], sigma_cur[j,,drop = FALSE], xi_cur[j,,drop = FALSE], pos, next_pos, acpt_lB[,bgrp_st:bgrp_end],
                  lB_prop_var[,bgrp_st:bgrp_end])
      }
    }
    for(j in 1:nrep){
      bgrp_st <- (j-1)*nbasis_per_rep + 1
      bgrp_end <- j*nbasis_per_rep
      lB_cur[,j] <- lBobj[[j]][["lB_cur"]]
      acpt_lB[,bgrp_st:bgrp_end] <- lBobj[[j]][["acpt_lB"]]
    }
    for(j in 1:nbasis){
      lB_prop_var[nxt_pos, j] <- update_var(lB_prop_var[pos, j],
                                            acpt_rt = mean(acpt_lB[,j], na.rm = TRUE),
                                            opt_rt = opt_rt, gamma1 = gamma1)
    }
    lA_cur <- mklA(lK_cur, lB_cur, alpha_cur)

    # Retain samples
    if(i%%thin_int ==0){
      smp_alpha[smp_loc] <- alpha_cur
      smp_theta[smp_loc] <- theta_cur
      smp_gscl[smp_loc] <- gscl_cur
      smp_gvar[smp_loc] <- gvar_cur
      smp_locpar[smp_loc,] <- locpar_cur
      smp_scalepar[smp_loc,] <- scalepar_cur
      smp_shapepar[smp_loc,] <- shapepar_cur
      if("mu"%in%mar_gp_which){
        smp_mu_gp[smp_loc,] <- mu_gp_cur
        smp_mu_gvar[smp_loc] <- mu_gvar_cur
        smp_mu_gscl[smp_loc] <- mu_gscl_cur
      }
      else{
        smp_mu_gp <- smp_mu_gvar <- smp_mu_gscl <- NULL
      }
      if("sigma"%in%mar_gp_which){
        smp_sigma_gp[smp_loc,] <- sigma_gp_cur
        smp_sigma_gvar[smp_loc] <- sigma_gvar_cur
        smp_sigma_gscl[smp_loc] <- sigma_gscl_cur
      }
      else{
        smp_sigma_gp <- smp_sigma_gvar <- smp_sigma_gscl <- NULL
      }
      if("xi"%in%mar_gp_which){
        smp_xi_gp[smp_loc,] <- xi_gp_cur
        smp_xi_gvar[smp_loc] <- xi_gvar_cur
        smp_xi_gscl[smp_loc] <- xi_gscl_cur
      }
      else{
        smp_xi_gp <- smp_xi_gvar <- smp_xi_gscl <- NULL
      }

      if(is.null(select_lat)){
        smp_lB[smp_loc,] <- c(lB_cur)
        smp_lK[smp_loc,,] <- lK_cur
      }
      else{
        for(j in 1:length(select_lat)){
          smp_lB[smp_loc,j] <- lB_cur[select_lat[j]]
          smp_lK[smp_loc,,j] <- lK_cur[,select_lat[j]]
        }
      }


      smp_loc <- smp_loc + 1
    }
  }
  # End main MH loop
  acpt_rts <- colMeans(acpt)
  acpt_rts_loc <- colMeans(acpt_loc)
  acpt_rts_scale <- colMeans(acpt_scale)
  acpt_rts_shape <- colMeans(acpt_shape)
  acpt_rts_lB <- colMeans(acpt_lB, na.rm = TRUE)
  acpt_rts_lK <- colMeans(acpt_lK, na.rm = TRUE)
  smp_mat <- cbind(smp_alpha, smp_theta, smp_gvar,smp_gscl, smp_locpar, smp_scalepar, smp_shapepar)
  colnames(smp_mat) <- c("alpha","theta","gvar","gscl",locnames, scalenames, shapenames)
  return(list(smp_par = smp_mat,
              smp_lB = smp_lB,
              smp_lK = smp_lK,
              smp_mu_gp = smp_mu_gp,
              smp_sigma_gp = smp_sigma_gp,
              smp_xi_gp = smp_xi_gp,
              smp_mu_gvar = smp_mu_gvar,
              smp_mu_gscl = smp_mu_gscl,
              smp_sigma_gvar = smp_sigma_gvar,
              smp_sigma_gscl = smp_sigma_gscl,
              smp_xi_gvar = smp_xi_gvar,
              smp_xi_gscl = smp_xi_gscl,
              acpt_rts = acpt_rts,
              acpt_rts_loc = acpt_rts_loc,
              acpt_rts_scale = acpt_rts_scale,
              acpt_rts_shape = acpt_rts_shape,
              acpt_rts_lB = acpt_rts_lB,
              acpt_rts_lK = acpt_rts_lK,
              prop_var = prop_var,
              loc_prop_var = loc_prop_var,
              scale_prop_var = scale_prop_var,
              shape_prop_var = shape_prop_var,
              final_lB_prop_var = lB_prop_var[pos,],
              final_lK_prop_var = lK_prop_var[pos,],
              clusters = clusters))
  }
