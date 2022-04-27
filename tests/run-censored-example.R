library(parallel)
library(dpsurrogate)

set.seed(3103)

boo <- generate_data(effect = "nonlinear")
niter <- 50

ldat <- boo$ldat

subtypes <- levels(ldat$J)

init <- initialize_mus(ldat, subtypes)

initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)

Xmat <- model.matrix(  ~ subtype + J:trt -1, data = ldat)
#logYY <- log(ldat$time)
SS <- ldat$ctDNA_fu

## add uniform censoring

CC <- runif(nrow(ldat), 0, 80)

boo$ldat$time.tilde <- pmin(ldat$time, CC)
boo$ldat$delta <- ifelse(ldat$time <= CC, 1, 0)
logYY <- log(boo$ldat$time.tilde)

#rstan_options(auto_write = TRUE)

cur <- init

dmat <- cbind(as.matrix(cur[, c("seff", "yeff")]), boo$rdat$trtZ)


mu0hat <- c(colSums(dmat[, -3] * boo$rdat$ngrp) / sum(boo$rdat$ngrp),
            mean(boo$rdat$trtZ))

v0hat <- c(colSums((dmat[, -3] - mu0hat[-3])^2 * boo$rdat$ngrp) / sum(boo$rdat$ngrp), var(boo$rdat$trtZ))

basepriors <- list(mu0 = mu0hat,
                   Lambda = diag(v0hat) * 4,
                   kappa0 = ncol(dmat),
                   nu = ncol(dmat))

dp <- DirichletProcessMvnormal(dmat,
                               g0Priors = basepriors,
                               alphaPriors = c(1, 2),
                               numInitialClusters = length(initclus$size))
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
                    t.cen = exp(logYY) + boo$ldat$delta,
                    is.censored = 1 - boo$ldat$delta,
                    s = SS,
                    X = Xmat,
                    prior_mu = prior.mu,
                    prior_sig = prior.sig)

  if(j == 1) {
    tmod <- jags.model(system.file("regmodel-cens.bug", package = "dpsurrogate"),
                       test_data, n.adapt = 100, quiet = TRUE)
    tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 50)


  } else {

    tmod <- jags.model(system.file("regmodel-cens.bug", package = "dpsurrogate"), test_data, inits = tmod$state(), n.adapt = 0, quiet = TRUE)
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
  dp <- Fit(dp, 50, progressBar = FALSE)
  dp <- UpdateAlpha(dp)

  cat(j, "\n")
}

allres <- list(res.s = res.s, res.y = res.y, jags.state = init.l,
     clusters = dp)



run_one_loo_cens <- function(lout, niter = 50, jags.state = NULL) {


  ldat <- boo$ldat
  subtypes <- levels(ldat$J)

  ldat$time[ldat$J == subtypes[lout]] <- NA

  init <- initialize_mus(ldat, subtypes, lout)

  init$yeff[lout] <- mean(init$yeff[-lout])

  initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)


  Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
  logYY <- log(ldat$time.tilde)
  SS <- ldat$ctDNA_fu

  #rstan_options(auto_write = TRUE)

  cur <- init

  dmat <- cbind(as.matrix(cur[, c("seff", "yeff")]), boo$rdat$trtZ)

  mu0hat <- c(colSums(dmat[, -3] * boo$rdat$ngrp) / sum(boo$rdat$ngrp),
              mean(boo$rdat$trtZ))

  v0hat <- c(colSums((dmat[, -3] - mu0hat[-3])^2 * boo$rdat$ngrp) / sum(boo$rdat$ngrp), var(boo$rdat$trtZ))

  basepriors <- list(mu0 = mu0hat,
                     Lambda = diag(v0hat) * 4,
                     kappa0 = ncol(dmat),
                     nu = ncol(dmat))

  dp <- DirichletProcessMvnormal(dmat,
                                 g0Priors = basepriors,
                                 alphaPriors = c(1, 2),
                                 numInitialClusters = length(initclus$size))

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
                      t.cen = exp(logYY) + boo$ldat$delta,
                      is.censored = 1 - boo$ldat$delta,
                      s = SS,
                      X = Xmat,
                      prior_mu = prior.mu,
                      prior_sig = prior.sig)


    if(j == 1) {

      if(!is.null(jags.state)) {
        tmod <- jags.model(system.file("regmodel-cens.bug", package = "dpsurrogate"),
                           inits = jags.state,
                           test_data, n.adapt = 0, quiet = TRUE)
        tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
      } else {
        tmod <- jags.model(system.file("regmodel-cens.bug", package = "dpsurrogate"),
                           test_data, n.adapt = 100, quiet = TRUE)
      }

      tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 50)


    } else {

      tmod <- jags.model(system.file("regmodel-cens.bug", package = "dpsurrogate"), test_data,
                         inits = tmod$state(), n.adapt = 0, quiet = TRUE)
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
    dp <- Fit(dp, 10, progressBar = FALSE)
    # dpTry <- tryCatch(Fit(dp, 10, progressBar = FALSE),
    #                   error = function(e) "singular")
    # if(identical(dpTry, "singular")) {
    #   next
    # } else {
    #   dp <- dpTry
    # }

    dp <- UpdateAlpha(dp)

    #cat(j, "\n")
  }

  list(res.s = res.s[, lout], res.y = res.y[, lout])

}



nnodes <- 6
cl <- makeCluster(nnodes)
init <- clusterEvalQ(cl, {
  library(dpsurrogate)
})

clusterExport(cl, c("boo", "run_one_loo_cens"), envir = environment())

resWho <- clusterApply(cl, 1:64, function(i) {

  #resWho <- lapply(cl, 1:64, function(i) {
  tryCatch({
    restmp <- run_one_loo_cens(i)
    data.frame(dpsur = restmp$res.y, leftout = i)
  }, error = function(e){

    list(error = e,
         lout = i
    )

  })

})

stopCluster(cl)

saveRDS(list(resall = allres, resloo = resWho), file = "censored-example.rds")

resloo <- do.call(rbind, resWho)

res1 <- list(trueeffects = boo$rdat, alleffects = allres,
             leaveouts = resloo, clusters = dp)
res1$trueeffects$group <- as.numeric(factor(res1$trueeffects$J))

## get average clustering according to dahl

getAvgCluster <- function(res1) {
  deltmats <- simplify2array(lapply(res1$clusters$labelsChain, function(del) {
    outer(del, del, FUN = function(x1, x2) as.numeric(x1 == x2))
  }))
  pihat <- apply(deltmats, c(1, 2), mean)
  sse <- apply((deltmats - simplify2array(lapply(1:length(res1$clusters$labelsChain), function(i) pihat)))^2,
               3, mean)
  res1$clusters$labelsChain[[which.min(sse)]]
}

res1$trueeffects$cluster <- factor(getAvgCluster(res1))


library(dplyr)
library(ggplot2)

tmp1 <- merge(res1$leaveouts, res1$trueeffects, by.x ="leftout", by.y = "group")

ests <- res1$leaveouts |> group_by(leftout) |> summarize(dpsur = median(dpsur))
ests <- merge(ests, res1$trueeffects, by.x = "leftout", by.y = "group")

ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = dpsur)) +
  geom_density2d(data = tmp1, aes(x = seff, y = dpsur, color = cluster),
                 n = 25, h = c(1, 1)) +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = dpsur, color = cluster)) +
  theme_bw() + ylim(c(-4, 4)) +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none", color = "none") + ggtitle("Nonlinear, with censoring",
                                                   subtitle = "DPM (clusters denoted by colors)")
ggsave("cens-examp.png", width = 5, height = 4.25)

