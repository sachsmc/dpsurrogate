library(dpsurrogate)
library(ggplot2)


getAvgCluster <- function(res1) {
  deltmats <- simplify2array(lapply(res1$clusters$labelsChain, function(del) {
    outer(del, del, FUN = function(x1, x2) as.numeric(x1 == x2))
  }))
  pihat <- apply(deltmats, c(1, 2), mean)
  sse <- apply((deltmats - simplify2array(lapply(1:length(res1$clusters$labelsChain), function(i) pihat)))^2,
               3, mean)
  res1$clusters$labelsChain[[which.min(sse)]]
}

summarize_sim <- function(path) {

  res1 <- readRDS(file.path("tests", "sim-res", path))

  setting <- paste(strsplit(path, "(-|_)")[[1]][c(2,3,4)], collapse = "-")

  res1$trueeffects$group <- as.numeric(factor(res1$trueeffects$J))
  clusters <- getAvgCluster(res1)

  clap <- unlist(lapply(1:64, function(i) {
    thisc <- res1$trueeffects$trueclusters[i]
    cdex <- which(res1$trueeffects$trueclusters == thisc)
    mean(res1$clusteri[[i]][50, i] == res1$clusteri[[i]][50, cdex])
    }))

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
  p.simp.true = mean(tmp1$D.simp < tmp1$D.null),
  nclusters = length(unique(clusters)),
  cluster.olap = mean(clap)
))

out1$setting <- setting
as.data.frame(out1)

}

filelist <- list.files("tests/sim-res/", "*.rds")
#filelist <- grep("manybiom|simple", filelist, value = TRUE)

allres <- lapply(filelist, function(d) {
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
                  "p.dpm", "p.simp", "p.dpm.true", "p.simp.true","nclusters",
                  "cluster.olap", "setting")


print(xtable(df[, .(mean(min), mean(Q1), mean(median), mean(Q3), mean(max)), by = .(setting)],
             digits = 1),
include.rownames = FALSE, sanitize.text.function = \(x) x)

msd <- function(x) {

  sprintf("%.2f (%.2f)", mean(x), sd(x))

}

print(xtable(df[, .(msd(cluster.olap),
                    msd(median.D.dpm.est), msd(median.D.simp.est), msd(median.D.dpm.true),
                    msd(median.D.simp.true)), by = .(setting)],
             digits = 2),
      include.rownames = FALSE, sanitize.text.function = \(x) x)



print(xtable(df[, .(msd(p.dpm), msd(p.simp), msd(p.dpm.true), msd(p.simp.true)), by = .(setting)]),
      sanitize.text.function = function(x) x, include.rownames = FALSE)
