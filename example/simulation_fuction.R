example_sim <- function(nmcmc, simset) {
  # Simulate data -----------------------------------------------------------
  s <- matrix(runif(simset$nloc * 2), ncol = 2)
  mar_betas <- list(loc = 0,
                    scale = 0,
                    shape = 0)
  zall <- rstabmix(
    simset$n,
    s,
    nbasis = simset$L,
    alpha = simset$alpha,
    delta = simset$alpha,
    theta = simset$theta,
    gvar = simset$gvar,
    gscl = simset$gscl,
    return_all = T,
    type = "br"
  )
  lZ <- zall[["lZ"]]
  lK <- zall[["lK"]]
  lA <- zall[["lA"]]
  lB <- zall[["lB"]]
  df <- data.frame(x = rep(1, simset$n * simset$nloc))
  mu <- mar_betas[["loc"]]
  sigma <- exp(mar_betas[["scale"]])
  xi <- mar_betas[["shape"]]
  y <- lsm2gev(lZ,simset$alpha, simset$theta, lK, mu, sigma, xi)
  

# Run MCMC ----------------------------------------------------------------
  init <-
    list(
      alpha = simset$alpha,
      theta = simset$theta,
      gvar = simset$gvar,
      gscl = simset$gscl,
      loc = mar_betas[["loc"]],
      scale = mar_betas[["scale"]],
      shape = mar_betas[["shape"]]
    )
  tune_var <- list(
    alpha = 0.01,
    theta = 0.01,
    gvar = 0.1,
    gscl = 0.1,
    loc = 0.1,
    scale = 0.1,
    shape = 0.1,
    lB = 0.1,
    lK = 0.1
  )
  run_time <- system.time(
    mh <- mcmc_lnorm_basis(
      nmcmc,
      s,
      y,
      lB,
      lK,
      init,
      tune_var,
      df,
      loc = ~ 1,
      scale = ~ 1,
      shape = ~ 1,
      thin_int = 1
    )
  )
  mh$run_time <- run_time
  mh$s <- s
  return(mh)
}
