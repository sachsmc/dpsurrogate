
#' Run a full data analysis (not leave one out) using just a bivariate normal hierarchical model
#' @export

run_one_simple_analysis <- function(boo, niter = 50) {

  #boo <- generate_data()

  ldat <- boo$ldat

  subtypes <- levels(ldat$J)

  Xmat <- model.matrix( ~ subtype + subtype:trt:trttype -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)

    Sig <- matrix(.05, nrow = 3, ncol = 3)
    diag(Sig) <- 2
    flexz <- c(rep(0, 80 - 64),boo$rdat$trtZ)

    prior.mu <- matrix(0, nrow = 80, ncol = 2)
    prior.sig <- array(NA, dim = c(2,2,80))

    for(i in 1:80){
      prior.mu[i, ] <- Sig[1:2,3] / Sig[3,3] * (flexz[i])
    }
    prior.sig <- Sig[1:2, 1:2] - Sig[1:2, 3] %*% solve(Sig[3, 3]) %*% Sig[3, 1:2]
    test_data <- list(J = ncol(Xmat),
                      N = nrow(Xmat),
                      y = logYY,
                      s = SS,
                      X = Xmat,
                      prior_mu = prior.mu,
                      prior_sig = prior.sig)


      tmod <- jags.model(system.file("normalmodel.bug", package = "dpsurrogate"),
                         test_data, n.adapt = 100, quiet = TRUE)
      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = niter)


    res.s[1:niter, ] <- t(as.matrix(tsam$beta_s[-c(1:16), , 1]))
    res.y[1:niter, ] <- t(as.matrix(tsam$beta_y[-c(1:16), , 1]))


  list(res.s = res.s, res.y = res.y)

}


#' Run a leave one out data analysis (not leave one out) using just a bivariate normal hierarchical model
#' @export

run_one_loo_simple_analysis <- function(boo, lout, niter = 50) {

  #boo <- generate_data()

  ldat <- boo$ldat

  subtypes <- levels(ldat$J)

  ldat$time[ldat$J == subtypes[lout]] <- NA

  Xmat <- model.matrix( ~ subtype + subtype:trt:trttype -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  res.s <- res.y <- rep(NA, 64)

  Sig <- matrix(.05, nrow = 3, ncol = 3)
  diag(Sig) <- 2
  flexz <- c(rep(0, 80 - 64),boo$rdat$trtZ)

  prior.mu <- matrix(0, nrow = 80, ncol = 2)
  prior.sig <- array(NA, dim = c(2,2,80))

  for(i in 1:80){
    prior.mu[i, ] <- Sig[1:2,3] / Sig[3,3] * (flexz[i])
  }
  prior.sig <- Sig[1:2, 1:2] - Sig[1:2, 3] %*% solve(Sig[3, 3]) %*% Sig[3, 1:2]
  test_data <- list(J = ncol(Xmat),
                    N = nrow(Xmat),
                    y = logYY,
                    s = SS,
                    X = Xmat,
                    prior_mu = prior.mu,
                    prior_sig = prior.sig)


  tmod <- jags.model(system.file("normalmodel.bug", package = "dpsurrogate"),
                     test_data, n.adapt = 100, quiet = TRUE)
  tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = niter)


  res.s <- t(as.matrix(tsam$beta_s[-c(1:16), , 1]))[, lout]
  res.y <- t(as.matrix(tsam$beta_y[-c(1:16), , 1]))[, lout]


  list(res.s = res.s, res.y = res.y)

}

