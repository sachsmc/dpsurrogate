
#' Run a full data analysis (not leave one out) using just a bivariate normal hierarchical model
#' @export

run_one_simple_analysis <- function(boo, niter = 50) {

  #boo <- generate_data(Zeffect = TRUE)

  ldat <- boo$ldat

  subtypes <- levels(ldat$J)

  Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)

    Sig <- matrix(.05, nrow = 2, ncol = 2)
    diag(Sig) <- 4

    test_data <- list(J = ncol(Xmat),
                      N = nrow(Xmat),
                      s = SS,
                      y = logYY,
                      X = Xmat,
                      Z = c(rep(0, 80 - 64), boo$rdat$trtZ),
                      prior_sig = Sig)


      tmod <- jags.model(system.file("normalmodel.bug", package = "dpsurrogate"),
                         test_data, n.adapt = 100, quiet = TRUE)
      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = niter)


    res.s[1:niter, ] <- t(as.matrix(tsam$beta_s[-c(1:16), , 1]))
    res.y[1:niter, ] <- t(as.matrix(tsam$beta_y[-c(1:16), , 1]))


  list(res.s = res.s, res.y = res.y, jags.state = tmod$state())

}


#' Run a leave one out data analysis (not leave one out) using just a bivariate normal hierarchical model
#' @export

run_one_loo_simple_analysis <- function(boo, lout, niter = 50, jags.state = NULL) {

  #boo <- generate_data(effect = "linear")

  ldat <- boo$ldat

  subtypes <- levels(ldat$J)

  ldat$time[ldat$J == subtypes[lout]] <- NA

  Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  res.s <- res.y <- rep(NA, 64)

  Sig <- matrix(.05, nrow = 2, ncol = 2)
  diag(Sig) <- 4

  test_data <- list(J = ncol(Xmat),
                    N = nrow(Xmat),
                    s = SS,
                    y = logYY,
                    X = Xmat,
                    Z = c(rep(0, 80 - 64), boo$rdat$trtZ),
                    prior_sig = Sig)


  if(!is.null(jags.state)) {
    tmod <- jags.model(system.file("normalmodel.bug", package = "dpsurrogate"),
                       inits = jags.state,
                       test_data, n.adapt = 0, quiet = TRUE)
    tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
  } else {
    tmod <- jags.model(system.file("normalmodel.bug", package = "dpsurrogate"),
                     test_data, n.adapt = 100, quiet = TRUE)
  }
  tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = niter)


  res.s <- t(as.matrix(tsam$beta_s[-c(1:16), , 1]))[, lout]
  res.y <- t(as.matrix(tsam$beta_y[-c(1:16), , 1]))[, lout]


  list(res.s = res.s, res.y = res.y)

}

