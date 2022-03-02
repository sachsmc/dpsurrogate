library(dpsurrogate)
library(ggplot2)

summarize_sim <- function(path) {

  res1 <- readRDS(file.path("tests", "sim-res", path))

  setting <- paste(strsplit(path, "-", fixed = TRUE)[[1]][c(2,4)], collapse = "-")

  res1$trueeffects$group <- as.numeric(factor(res1$trueeffects$J))
  check.errors <- sapply(res1$leaveouts, function(x) !is.null(x$error))

if( any(check.errors)){

  for(j in which(check.errors)) {
    idat <- res1$leaveouts[[j]]$dset
    niter <- 50
    restmp <- run_one_loo(idat, j, niter = niter, jags.state = NULL)
    restmp2 <- run_one_loo_null(idat, j, niter = niter, jags.state = NULL)
    restmp3 <- run_one_loo_simple_analysis(idat, j, niter = niter, jags.state = NULL)
    res1$leaveouts[[j]] <- data.frame(dpsur = restmp$res.y, null = restmp2$res.y,
               simple = restmp3$res.y, leftout = i)

  }

}
res1$leaveouts <- do.call(rbind, res1$leaveouts)

ssize <- quantile(res1$trueeffects$ngrp, c(0, .25, .5, .75, 1))
# P(D_a > D_0)

tmp <- merge(res1$leaveouts, res1$alleffects, by.x = "leftout", by.y = "group")
tmp1 <- merge(res1$leaveouts, res1$trueeffects, by.x ="leftout", by.y = "group")


tmp$D.dpm <- abs(tmp$dpsur - tmp$dpsur.yeff)
tmp$D.simp <- abs(tmp$simple - tmp$dpsur.yeff)
tmp$D.null <- abs(tmp$null - tmp$dpsur.yeff)

tmp1$D.dpm <- abs(tmp1$dpsur - tmp1$yeff)
tmp1$D.simp <- abs(tmp1$simple - tmp1$yeff)
tmp1$D.null <- abs(tmp1$null - tmp1$yeff)


out1 <- as.list(c(  ssize,
  sapply(tmp[, c("D.dpm", "D.simp", "D.null")], median),
  sapply(tmp1[, c("D.dpm", "D.simp", "D.null")], median),
  p.dpm = mean(tmp$D.dpm < tmp$D.null),
  p.simp = mean(tmp$D.simp < tmp$D.null),
  p.dpm.true = mean(tmp1$D.dpm < tmp1$D.null),
  p.simp.true = mean(tmp1$D.simp < tmp1$D.null)
))

out1$setting <- setting
as.data.frame(out1)

}

allres <- lapply(list.files("tests/sim-res/", "*.rds"), function(d) {
  tryCatch(summarize_sim(d), error = function(e) e)
  })

saveRDS(allres, file = "tests/compiled-sims.rds")


df <- do.call(rbind, allres)

library(data.table)
library(xtable)

df <- data.table(df)

colnames(df) <- c("min", "Q1", "median", "Q3", "max",
                  "median.D.dpm.est", "median.D.simp.est", "median.D.null.est",
                  "median.D.dpm.true", "median.D.simp.true", "median.D.null.true",
                  "p.dpm", "p.simp", "p.dpm.true", "p.simp.true", "setting")


setlo <- c("linear, $c_z = 0$", "linear, $c_z = 0.25$", "nonlinear, $c_z = 0$", "nonlinear, $c_z = 0.25$",
           "null, $c_z = 0$", "null, $c_z = 0.25$")
names(setlo) <- c("lin-0", "lin-1", "nonlin-0", "nonlin-1", "null-0", "null-1")

df$setting <- setlo[df$setting]

print(xtable(df[, .(mean(min), mean(Q1), mean(median), mean(Q3), mean(max)), by = .(setting)],
             digits = 1),
include.rownames = FALSE)

msd <- function(x) {

  sprintf("%.2f (%.2f)", mean(x), sd(x))

}

print(xtable(df[, .(msd(median.D.dpm.est), msd(median.D.simp.est), msd(median.D.dpm.true),
                    msd(median.D.simp.true)), by = .(setting)],
             digits = 2),
      include.rownames = FALSE)



print(xtable(df[, .(msd(p.dpm), msd(p.simp), msd(p.dpm.true), msd(p.simp.true)), by = .(setting)]),
      sanitize.text.function = function(x) x, include.rownames = FALSE)
