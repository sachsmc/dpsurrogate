library(dpsurrogate)
library(ggplot2)


idat <- generate_data()

resall <- run_one_analysis(idat)
#resnull <- run_one_null(idat, niter = 10)
ressimp <- run_one_simple_analysis(idat)
plot(idat$rdat$yeff ~ idat$rdat$seff)

#coo1 <- rep(1, 64)
#coo1[lout] <- 2
#points(colMeans(res.y) ~ colMeans(res.s), pch = 20, col = coo1)
 points(colMeans(resall$res.y) ~ colMeans(resall$res.s), pch = 20, col = resall$clusters)
 points(colMeans(ressimp$res.y) ~ colMeans(ressimp$res.s), pch = 20, col = "blue")
#
# mean((colMeans(resall$res.y) - idat$rdat$yeff)^2)
# mean((colMeans(resnull$res.y) - idat$rdat$yeff)^2)

resloo.s <- resloo.y <- matrix(NA, nrow = nrow(resall$res.s), ncol = ncol(resall$res.s))
resloo.s.n <- resloo.y.n <- matrix(NA, nrow = 50, ncol = ncol(resnull$res.s))
for(i in 1:nrow(idat$rdat)) {
  restmp <- run_one_loo(idat, i)
  restmp2 <- run_one_loo_null(idat, i)
  resloo.s[, i] <- restmp$res.s
  resloo.y[, i] <- restmp$res.y
  resloo.y.n[, i] <- restmp2$res.y
  cat(i, "\n")
}


saveRDS(list(resall = resall, resnull = resnull,
             resloo.s = resloo.s, resloo.y = resloo.y,
             resloo.y.n = resloo.y.n), file = "testsim3.rds")

tim <- readRDS("testsim3.rds")

idat$rdat <- cbind(idat$rdat, do.call(rbind, strsplit(idat$rdat$J, ":")))
idat$rdat <- idat$rdat[order(idat$rdat$`2`, idat$rdat$`1`),]

mean(abs(colMeans(resloo.y) - idat$rdat$yeff))
mean(abs(colMeans(resloo.y.n) - idat$rdat$yeff))

mean(abs(colMeans(resloo.y) - colMeans(resall$res.y)))
mean(abs(colMeans(resloo.y.n) - colMeans(resall$res.y)))

xord <- order(sapply(strsplit(levels(idat$ldat$J), ":"), function(x) paste0(x[2], ":", x[1])))

cheez <- merge(data.frame(yeff = colMeans(tim$resall$res.y),
                          seff = colMeans(tim$resall$res.s),
                          J = levels(idat$ldat$J)[xord]),
               data.frame(loo.yeff = colMeans(tim$resloo.y),
                          loo.yeff.n = colMeans(tim$resloo.y.n),
                          J = levels(idat$ldat$J)), by = "J")

mean(abs(cheez$loo.yeff - cheez$yeff))
mean(abs(cheez$loo.yeff.n - cheez$yeff))

library(ggplot2)


ggplot(cheez, aes(x = seff, y = yeff, xend = seff, yend = loo.yeff)) +
  geom_segment() +
  geom_point(shape = 1) +
  geom_point(aes(x = seff, y = loo.yeff)) + theme_bw() +
  xlab("estimated nu") + ylab("estimated mu") +
  guides(alpha = "none") + ggtitle("Leave one outs (filled) vs full (open)")
ggsave("sim-full.png")

ggplot(cheez, aes(x = seff, y = yeff, xend = seff, yend = loo.yeff.n)) +
  geom_segment() +
  geom_point(shape = 1) +
  geom_point(aes(x = seff, y = loo.yeff.n)) + theme_bw() +
  xlab("estimated nu") + ylab("estimated mu") +
  guides(alpha = "none")+ ggtitle("Null model: Leave one outs (filled) vs full (open)")
ggsave("sim-null.png")
