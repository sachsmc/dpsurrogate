library(ggplot2)
library(parallel)
library(dpsurrogate)

set.seed(401)

boo <- generate_data("twotrt")
lout <- 9


ldat <- boo$ldat
subtypes <- levels(ldat$J)

ldat$time[ldat$J == subtypes[lout]] <- NA

init <- initialize_mus(ldat, subtypes, lout)

init$yeff[lout] <- mean(init$yeff[-lout])

initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)


Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
SS <- ldat$ctDNA_fu

## add uniform censoring

CC <- runif(nrow(ldat), 20, 60)
ldat$CC <- CC

boo$ldat$time.tilde <- pmin(ldat$time, CC)
boo$ldat$delta <- ifelse(ldat$time <= CC, 1, 0)
#boo$ldat$delta[is.na(ldat$time)] <- 0
logYY <- log(boo$ldat$time.tilde)

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


niter <- 100

labelschain <- lapply(1:(niter * 10 + 50), function(x) NA)
labelschain[1:50] <- dp$labelsChain[151:200]
lch <- 51

res.s <- res.y <- matrix(NA, nrow = niter, ncol = 64)
jags.state <- NULL

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
                    OO = which(boo$ldat$delta == 1 | is.na(boo$ldat$delta)),
                    JJ = which(boo$ldat$delta == 0),
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
  labelschain[lch:(lch + 9)] <- dp$labelsChain
  lch <- lch + 10

  dpmu <- dp$clusterParameters$mu[,, dp$clusterLabels[lout]]
  dpsig <- dp$clusterParameters$sig[,, dp$clusterLabels[lout]]

  ppmu <- dpmu[2] + (dpsig[c(2), c(1, 3)] %*% solve(dpsig[c(1,3), c(1,3)])) %*%
    (rbind(res.s[j,lout], boo$rdat$trtZ[lout]) - t(t(dpmu[c(1,3)])))
  ppsig <- sqrt(c(dpsig[2, 2] - dpsig[c(2), c(1, 3)] %*%
                    solve(dpsig[c(1,3), c(1,3)]) %*% dpsig[c(1, 3), c(2)]))

  res.y[j, lout]  <- rnorm(1, mean = ppmu, sd = ppsig)

  cat(j, "\n")
}


deltmats <- simplify2array(lapply(labelschain, function(del) {
    outer(del, del, FUN = function(x1, x2) as.numeric(x1 == x2))
  }))
  pihat <- apply(deltmats, c(1, 2), mean)
  sse <- apply((deltmats - simplify2array(lapply(1:length(labelschain), function(i) pihat)))^2,
               3, mean)
finclust <- labelschain[[which.min(sse)]]
finclust[9]

col <- rep("black", 64)
col[lout] <- "red"
plot(yeff ~ seff, data = boo$rdat, pch = 20, ylim = c(-7, 4), col = col)
points(yeff[lout]~ seff[lout], data = boo$rdat, col = "red")
points(res.y[, 9] ~ res.s[, 9])

est.effs <- data.frame(yeff = apply(res.y, 2, median, na.rm = TRUE),
                       seff = apply(res.s, 2, median, na.rm = TRUE),
                       observed = ifelse(1:64 != 9, "observed", "unobserved"),
                       cluster = factor(finclust))

ggplot() +
#  geom_density_2d(data = data.frame(yeff = res.y[, 9], seff = res.s[, 9]),
#                 aes(x = seff, y = yeff), color = "grey65") +
  geom_point(data = est.effs, aes(x = seff, y = yeff, color = cluster,
                                  shape = observed), size = 2) +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  theme_bw() +
  geom_rug(data = est.effs[9,], aes(x = seff)) +
#  geom_linerange(data = cbind(est.effs[9,]), aes(x = seff, ymin = yeff - .926, ymax = yeff + .926)) +
  scale_color_manual("cluster", values = c("salmon", "slateblue")) +
  scale_shape_manual("", values = c(19, 1))

ggsave("example.png", width = 5.25, height = 3.5)


## leave one outs over the 63 observed


cens_loo <- function(lout) {

ldat <- subset(ldat, J != "--+-:ARSi")
ldat$J <- factor(ldat$J)

subtypes <- levels(ldat$J)

ldat$time[ldat$J == subtypes[lout]] <- NA

init <- initialize_mus(ldat, subtypes, lout)

init$yeff[lout] <- mean(init$yeff[-lout])

initclus <- kmeans(as.matrix(init[, 1:2]), centers = 10)


Xmat <- model.matrix( ~ subtype + J:trt -1, data = ldat)
SS <- ldat$ctDNA_fu

## add uniform censoring

time.tilde <- pmin(ldat$time, ldat$CC)
delta <- ifelse(ldat$time <= ldat$CC, 1, 0)
logYY <- log(time.tilde)

cur <- init

boo$rdat <- subset(boo$rdat, J != "--+-:ARSi")
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


niter <- 100

res.s <- res.y <- matrix(NA, nrow = niter, ncol = 63)
jags.state <- NULL

for(j in  1:niter){

  upriors <- update_priors(dp, boo)

  prior.mu <- rbind(matrix(0, nrow = 16, ncol = 2), upriors$prior.mu)
  prior.sig0 <- as.array(c(lapply(1:16, function(i) matrix(c(2, .05, .05, 2), nrow = 2)),
                           upriors$prior.sig))
  prior.sig <- array(NA, dim = c(2,2,79))
  for(i in 1:79){
    prior.sig[,,i] <- prior.sig0[[i]]
  }



  test_data <- list(J = ncol(Xmat),
                    N = nrow(Xmat),
                    y = logYY,
                    t.cen = exp(logYY) + delta,
                    is.censored = 1 - delta,
                    OO = which(delta == 1 | is.na(delta)),
                    JJ = which(delta == 0),
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

  dpmu <- dp$clusterParameters$mu[,, dp$clusterLabels[lout]]
  dpsig <- dp$clusterParameters$sig[,, dp$clusterLabels[lout]]

  ppmu <- dpmu[2] + (dpsig[c(2), c(1, 3)] %*% solve(dpsig[c(1,3), c(1,3)])) %*%
    (rbind(res.s[j,lout], boo$rdat$trtZ[lout]) - t(t(dpmu[c(1,3)])))
  ppsig <- sqrt(c(dpsig[2, 2] - dpsig[c(2), c(1, 3)] %*%
                    solve(dpsig[c(1,3), c(1,3)]) %*% dpsig[c(1, 3), c(2)]))

  res.y[j, lout]  <- rnorm(1, mean = ppmu, sd = ppsig)

  #cat(j, "\n")
}

  res.y[, lout]

}

cl <- makeCluster(7)
init <- clusterEvalQ(cl, {
  library(dpsurrogate)
})

clusterExport(cl, c("ldat", "boo", "cens_loo"), envir = environment())

resWho <- clusterApply(cl, 1:63, function(i) {
  cens_loo(i)
})



ldat2 <- subset(ldat, J != "--+-:ARSi")
ldat2$J <- factor(ldat2$J)

subtypes <- levels(ldat2$J)


resWho <- readRDS("reswho.rds")

resall <- res.y[, -9]
resall <- data.frame(yeff = c(resall), set = rep(1:63, each = 100))
resloo <- data.frame(yeff.loo = unlist(resWho), set = rep(1:63, each = 10))

dtoo <- merge(resall, resloo, by = "set")

dtoo2 <- merge(dtoo, data.frame(cluster = finclust[-9], set = 1:63), by = "set")

dtoo2$lild <- (abs(dtoo2$yeff - dtoo2$yeff.loo))

median(dtoo2$lild)
by(dtoo2$lild, dtoo2$cluster, median)

dtoo2$null <- sample(c(res.y), nrow(dtoo2), replace = TRUE)

dtoo2$Dnull <- abs(dtoo2$yeff - dtoo2$null)

mean(dtoo2$lild < dtoo2$Dnull)
with(subset(dtoo2, cluster == 1), mean(lild < Dnull))
with(subset(dtoo2, cluster == 2), mean(lild < Dnull))
