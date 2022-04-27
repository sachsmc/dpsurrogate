
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
logYY <- log(ldat$time)
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


niter <- 100

labelschain <- lapply(1:(niter * 10 + 50), function(x) NA)
labelschain[1:50] <- dp$labelsChain[151:200]
lch <- 51

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

    if(!is.null(jags.state)) {
      tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"),
                         inits = jags.state,
                         test_data, n.adapt = 0, quiet = TRUE)
      tmp <- adapt(tmod, n.iter = 0, end.adaptation = TRUE)
    } else {
      tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"),
                         test_data, n.adapt = 100, quiet = TRUE)
    }

    tsam <- jags.samples(tmod, c("beta_s", "beta_y"), n.iter = 50)


  } else {

    tmod <- jags.model(system.file("regmodel.bug", package = "dpsurrogate"), test_data,
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
  labelschain[lch:(lch + 9)] <- dp$labelsChain
  lch <- lch + 10
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

est.effs <- data.frame(yeff = apply(res.y, 2, median),
                       seff = apply(res.s, 2, median),
                       observed = 1:64 != 9,
                       cluster = factor(finclust))

ggplot() +
  geom_density_2d(data = data.frame(yeff = res.y[, 9], seff = res.s[, 9]),
                 aes(x = seff, y = yeff), color = "grey65") +
  geom_point(data = est.effs, aes(x = seff, y = yeff, color = cluster,
                                  shape = observed), size = 2) +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  theme_bw() + scale_shape_manual(values = c(1,  19)) +
  geom_rug(data = est.effs[9,], aes(x = seff)) +
  geom_linerange(data = cbind(est.effs[9,]), aes(x = seff, ymin = yeff - .926, ymax = yeff + .926)) +
  scale_color_manual("cluster", values = c("salmon", "slateblue"))

ggsave("example.png", width = 5.25, height = 3.5)

res2 <- readRDS("tests/sim-res/runout-twotrt-0.0-0.0_16.rds")

getAvgCluster <- function(res1) {
  deltmats <- simplify2array(lapply(res1$clusters$labelsChain, function(del) {
    outer(del, del, FUN = function(x1, x2) as.numeric(x1 == x2))
  }))
  pihat <- apply(deltmats, c(1, 2), mean)
  sse <- apply((deltmats - simplify2array(lapply(1:length(res1$clusters$labelsChain), function(i) pihat)))^2,
               3, mean)
  res1$clusters$labelsChain[[which.min(sse)]]
}

res2$trueeffects$group <- as.numeric(factor(res2$trueeffects$J))
res2$trueeffects$cluster <- factor(getAvgCluster(res2))

tmp2 <- merge(res2$leaveouts, res2$alleffects, by.x = "leftout", by.y = "group")
tmp21 <- merge(res2$leaveouts, res2$trueeffects, by.x ="leftout", by.y = "group")

ests2 <- res2$leaveouts |> group_by(leftout) |> summarize(dpsur = median(dpsur),
                                                          null = median(null),
                                                          simple = median(simple))
ests2 <- merge(ests2, res2$trueeffects, by.x = "leftout", by.y = "group")


cheez2 <-  merge(res2$leaveouts, res2$trueeffects, by.x = "leftout", by.y = "group")

cheez2$D <- abs(cheez2$dpsur - cheez2$yeff)
cheez2$D0 <- abs(cheez2$null - cheez2$yeff)
cheez2$Ds <- abs(cheez2$simple - cheez2$yeff)


by(cheez2$D, cheez2$cluster, summary)

with(cheez2, mean(D < D0))
with(subset(cheez2, cluster == 1), mean(D < D0))
with(subset(cheez2, cluster == 2), mean(D < D0))
