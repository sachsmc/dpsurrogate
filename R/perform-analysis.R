#' Get initial estimates of treatment effects
#' @export


initialize_mus <- function(ldat, subtypes, lout = NULL) {

  init <- NULL
  for(j in subtypes) {

    # no.pool <- stan_glm(ctDNA_fu ~ trt, data = subset(ldat, J == j),
    #                     family = gaussian, seed = 20201105,
    #                     chains = 1, iter = 2000, refresh = 0)

    no.pool <- glm(ctDNA_fu ~ trt, data = subset(ldat, J == j),
                        family = gaussian)

    if(!is.null(lout) && j == subtypes[lout]) {

      init <- rbind(init, data.frame(seff = coef(no.pool)[-1], yeff = NA, J = j))

    } else {

      # no.pool2 <- stan_glm(log(time) ~ trt, data = subset(ldat, J == j),
      #                      family = gaussian, seed = 20201105,
      #                      chains = 1, iter = 2000, refresh = 0)

      no.pool2 <- glm(log(time) ~ trt, data = subset(ldat, J == j),
                           family = gaussian)

      init <- rbind(init, data.frame(seff = coef(no.pool)[-1], yeff = coef(no.pool2)[-1], J = j))
    }
  }

  init
}

#' Update priors
#' @export
update_priors <- function(dp, boo) {
  ## compute conditional (nu, mu) | Z based on these

  if(ncol(dp$data) == 3) {
    condpriors <- lapply(1:length(dp$clusterLabels), function(i) {

      Sig <- dp$clusterParameters$sig[,, dp$clusterLabels[i]]

      newmu <- dp$clusterParameters$mu[, , dp$clusterLabels[i]][1:2] +
        Sig[1:2,3] / Sig[3,3] * (boo$rdat$trtZ[i] - dp$clusterParameters$mu[, , dp$clusterLabels[i]][3])

      newsig <- Sig[1:2, 1:2] - Sig[1:2, 3] %*% solve(Sig[3, 3]) %*% Sig[3, 1:2]

      list(newmu = newmu, newsig = newsig)

    })


    prior.mu <-  t(sapply(condpriors, function(i) i$newmu))
    prior.sig <- lapply(condpriors, function(i) i$newsig)

  } else if(ncol(dp$data) == 2) {
    ## compute conditional (mu) | Z based on these

    condpriors <- lapply(1:length(dp$clusterLabels), function(i) {

      Sig <- dp$clusterParameters$sig[,, dp$clusterLabels[i]]

      newmu <- dp$clusterParameters$mu[, , dp$clusterLabels[i]][1] +
        Sig[1,2] / Sig[2,2] * (boo$rdat$trtZ[i] - dp$clusterParameters$mu[, , dp$clusterLabels[i]][2])

      newsig <- Sig[1, 1] - Sig[1, 2] %*% solve(Sig[2, 2]) %*% Sig[2, 1]

      list(newmu = newmu, newsig = newsig)

    })


    prior.mu <-  matrix(sapply(condpriors, function(i) i$newmu), ncol = 1)
    prior.sig <- lapply(condpriors, function(i) i$newsig)

  } else {
    stop("Can't update priors")
  }


  list(prior.mu = prior.mu, prior.sig = prior.sig)
}

#' Run a full data analysis (not leave one out)
#' @export

run_one_analysis <- function(boo, niter = 50) {

  #boo <- generate_data()

  ldat <- boo$ldat

  subtypes <- levels(ldat$J)

  init <- initialize_mus(ldat, subtypes)

  initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)

  Xmat <- model.matrix( ~ subtype + subtype:trt:trttype -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  #rstan_options(auto_write = TRUE)

  cur <- init

  dmat <- cbind(as.matrix(cur[, c("seff", "yeff")]), boo$rdat$trtZ)
  dp <- DirichletProcessMvnormal(dmat, alphaPriors = c(2, 4), numInitialClusters = length(initclus$size))
  dp <- Fit(dp, 200, progressBar = FALSE)
  dp <- UpdateAlpha(dp)

  #niter <- 50

  res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)

  for(j in  1:niter){

    upriors <- update_priors(dp, boo)

    prior.mu <- rbind(matrix(0, nrow = 16, ncol = 2), upriors$prior.mu)
    prior.sig0 <- as.array(c(lapply(1:16, function(i) matrix(c(2, .05, .05, 2), nrow = 2)),
                             upriors$prior.sig))
    prior.sig <- array(NA, dim = c(2,2,80))
    for(i in 1:80){
      prior.sig[,,i] <- prior.sig0[[i]]
    }

    test_data <- list(J = ncol(Xmat),
                      N = nrow(Xmat),
                      y = logYY,
                      s = SS,
                      X = Xmat,
                      prior_mu = prior.mu,
                      prior_sig = prior.sig)

    if(j == 1) {
      tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"), test_data, n.adapt = 100, quiet = TRUE)
      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 50)


    } else {

      tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"), test_data, inits = tmod$state(), n.adapt = 0, quiet = TRUE)
      tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 10)


    }


    init.l <- tmod$state()

    res.s[j, ] <- as.matrix(tsam$beta_s[-c(1:16), dim(tsam$beta_s)[2], 1])
    res.y[j, ] <- as.matrix(tsam$beta_y[-c(1:16), dim(tsam$beta_y)[2], 1])

    newdmat <- as.matrix(cbind(seff = res.s[j, ],
                               yeff = res.y[j, ]))
    dmatin <- cbind(newdmat, boo$rdat$trtZ)
    dp <- ChangeObservations(dp, dmatin)
    dpTry <- tryCatch(Fit(dp, 50, progressBar = FALSE),
                      error = function(e) "singular")
    if(identical(dpTry, "singular")) {
      next
    } else {
      dp <- dpTry
    }
    dp <- UpdateAlpha(dp)

    #cat(j, "\n")
  }

  list(res.s = res.s, res.y = res.y)

}

#' Run one analysis leaving one subgroup out
#' @export

run_one_loo <- function(boo, lout, niter = 50) {

  #boo <- generate_data()

  ldat <- boo$ldat
  subtypes <- levels(ldat$J)

  ldat$time[ldat$J == subtypes[lout]] <- NA

  init <- initialize_mus(ldat, subtypes, lout)

  init$yeff[lout] <- mean(init$yeff[-lout])

  initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)


  Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  #rstan_options(auto_write = TRUE)

  cur <- init

  dmat <- cbind(as.matrix(cur[, c("seff", "yeff")]), boo$rdat$trtZ)
  dp <- DirichletProcessMvnormal(dmat, alphaPriors = c(2, 4), numInitialClusters = length(initclus$size))
  dp <- Fit(dp, 200, progressBar = FALSE)
  dp <- UpdateAlpha(dp)

  #niter <- 50

  res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)

  for(j in  1:niter){

    upriors <- update_priors(dp, boo)

    prior.mu <- rbind(matrix(0, nrow = 16, ncol = 2), upriors$prior.mu)
    prior.sig0 <- as.array(c(lapply(1:16, function(i) matrix(c(2, .05, .05, 2), nrow = 2)),
                             upriors$prior.sig))
    prior.sig <- array(NA, dim = c(2,2,80))
    for(i in 1:80){
      prior.sig[,,i] <- prior.sig0[[i]]
    }

    test_data <- list(J = ncol(Xmat),
                      N = nrow(Xmat),
                      y = logYY,
                      s = SS,
                      X = Xmat,
                      prior_mu = prior.mu,
                      prior_sig = prior.sig)

    if(j == 1) {
      tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"), test_data, n.adapt = 100, quiet = TRUE)
      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 50)


    } else {

      tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"), test_data, inits = tmod$state(), n.adapt = 0, quiet = TRUE)
      tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 10)


    }



    init.l <- tmod$state()

    res.s[j, ] <- as.matrix(tsam$beta_s[-c(1:16), dim(tsam$beta_s)[2], 1])
    res.y[j, ] <- as.matrix(tsam$beta_y[-c(1:16), dim(tsam$beta_y)[2], 1])

    ld <- lout + 16
    res.y[j, lout]  <- rnorm(1, mean = prior.mu[ld,2] +
                               sqrt(prior.sig[2,2,ld])/sqrt(prior.sig[1,1,ld]) *
                               prior.sig[1,2,ld] / (sqrt(prior.sig[2,2,ld]) *
                                                      sqrt(prior.sig[1,1,ld])) *
                               (res.s[j,lout] - prior.mu[ld,1]),
                             sd = (1 - prior.sig[1,2,ld] /
                                     (sqrt(prior.sig[2,2,ld]) * sqrt(prior.sig[1,1,ld]))) *
                               sqrt(prior.sig[2,2,ld]) )

    newdmat <- as.matrix(cbind(seff = res.s[j, ],
                               yeff = res.y[j, ]))
    dmatin <- cbind(newdmat, boo$rdat$trtZ)
    dp <- ChangeObservations(dp, dmatin)

    dpTry <- tryCatch(Fit(dp, 40, progressBar = FALSE),
                      error = function(e) "singular")
    if(identical(dpTry, "singular")) {
      next
    } else {
      dp <- dpTry
    }

    dp <- UpdateAlpha(dp)

    #cat(j, "\n")
  }

  list(res.s = res.s[, lout], res.y = res.y[, lout])

}

#' Run one null model analysis
#' @export

run_one_null <- function(boo, niter = 50) {

  #boo <- generate_data()

  ldat <- boo$ldat

  subtypes <- levels(ldat$J)

  init <- initialize_mus(ldat, subtypes)

  initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)


  Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  # rstan_options(auto_write = TRUE)

  cur <- init

  dmat <- cbind(as.matrix(cur[, c("yeff")]), boo$rdat$trtZ)
  dp <- DirichletProcessMvnormal(dmat, alphaPriors = c(2, 4), numInitialClusters = length(initclus$size))
  dp <- Fit(dp, 200, progressBar = FALSE)
  dp <- UpdateAlpha(dp)

  #niter <- 50

  res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)

  for(j in  1:niter){


    upriors <- update_priors(dp, boo)

    prior.mu <- c(matrix(0, nrow = 16, ncol = 1), upriors$prior.mu)
    prior.sig <- c(rep(2, 16), unlist(upriors$prior.sig))

    test_data <- list(J = ncol(Xmat),
                      N = nrow(Xmat),
                      y = logYY,
                      s = SS,
                      X = Xmat,
                      prior_mu = prior.mu,
                      prior_sig = prior.sig)

    if(j == 1) {
      tmod <- jags.model(system.file("regmodel-null.bug", package = "dpsurrogate"), test_data, n.adapt = 100, quiet = TRUE)
      tsam <- jags.samples(tmod, c("beta_y"), n.iter = 50)


    } else {

      tmod <- jags.model(system.file("regmodel-null.bug", package = "dpsurrogate"), test_data, inits = tmod$state(), n.adapt = 0, quiet = TRUE)
      tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
      tsam <- jags.samples(tmod, c("beta_y"), n.iter = 10)


    }


    init.l <- tmod$state()

    #res.s[j, ] <- as.matrix(fit1, pars = "beta_s")[nrow(as.matrix(fit1)), -c(1:16)]
    res.y[j, ] <- as.matrix(tsam$beta_y[-c(1:16), dim(tsam$beta_y)[2], 1])

    newdmat <- as.matrix(cbind(yeff = res.y[j, ]))
    dmatin <- cbind(newdmat, boo$rdat$trtZ)
    dp <- ChangeObservations(dp, dmatin)
    dpTry <- tryCatch(Fit(dp, 50, progressBar = FALSE),
                      error = function(e) "singular")
    if(identical(dpTry, "singular")) {
      next
    } else {
      dp <- dpTry
    }
    dp <- UpdateAlpha(dp)

    #cat(j, "\n")
  }

  list(res.s = res.s, res.y = res.y)

}

#' Run one null model leaving out trial lout
#' @export

run_one_loo_null <- function(boo, lout, niter = 50) {

  #boo <- generate_data()

  ldat <- boo$ldat
  subtypes <- levels(ldat$J)

  ldat$time[ldat$J == subtypes[lout]] <- NA

  init <- initialize_mus(ldat, subtypes, lout)

  init$yeff[lout] <- mean(init$yeff[-lout])

  initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)

  Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
  logYY <- log(ldat$time)
  SS <- ldat$ctDNA_fu

  # rstan_options(auto_write = TRUE)

  cur <- init

  dmat <- cbind(as.matrix(cur[, c("yeff")]), boo$rdat$trtZ)
  dp <- DirichletProcessMvnormal(dmat, alphaPriors = c(2, 4), numInitialClusters = length(initclus$size))
  dp <- Fit(dp, 200, progressBar = FALSE)
  dp <- UpdateAlpha(dp)

  #niter <- 50

  res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)

  for(j in  1:niter){

    upriors <- update_priors(dp, boo)
    prior.mu <- c(matrix(0, nrow = 16, ncol = 1), upriors$prior.mu)
    prior.sig <- c(rep(2, 16), unlist(upriors$prior.sig))

    test_data <- list(J = ncol(Xmat),
                      N = nrow(Xmat),
                      y = logYY,
                #      s = SS,
                      X = Xmat,
                      prior_mu = prior.mu,
                      prior_sig = prior.sig)

    if(j == 1) {
      tmod <- jags.model(system.file("regmodel-null.bug", package = "dpsurrogate"), test_data, n.adapt = 100, quiet = TRUE)
      tsam <- jags.samples(tmod, c("beta_y"), n.iter = 50)


    } else {

      tmod <- jags.model(system.file("regmodel-null.bug", package = "dpsurrogate"), test_data, inits = tmod$state(), n.adapt = 0, quiet = TRUE)
      tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
      tsam <- jags.samples(tmod, c("beta_y"), n.iter = 10)


    }

    res.y[j, ] <- as.matrix(tsam$beta_y[-c(1:16), dim(tsam$beta_y)[2], 1])
    res.y[j, lout] <- rnorm(1, mean = prior.mu[lout + 16],
                            sd = prior.sig[lout + 16])

    newdmat <- as.matrix(cbind(
      yeff = res.y[j, ]))
    dmatin <- cbind(newdmat, boo$rdat$trtZ)
    dp <- ChangeObservations(dp, dmatin)
    dpTry <- tryCatch(Fit(dp, 50, progressBar = FALSE),
                      error = function(e) "singular")
    if(identical(dpTry, "singular")) {
      next
    } else {
      dp <- dpTry
    }
    dp <- UpdateAlpha(dp)

    #cat(j, "\n")
  }

  list(res.s = res.s[, lout], res.y = res.y[, lout])

}


