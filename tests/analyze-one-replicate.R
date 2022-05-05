library(dplyr)
library(ggplot2)
library(patchwork)
library(dpsurrogate)

res1 <- readRDS("tests/sim-res/runout-twotrt-0.0-0.0_01.rds")
res2 <- readRDS("tests/sim-res/runout-manybiom-0.0-0.0_00.rds")

name1 <- "Good for two treatments"
name2 <- "Differential by biomarker"

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

tmp <- merge(res1$leaveouts, res1$alleffects, by.x = "leftout", by.y = "group")
tmp1 <- merge(res1$leaveouts, res1$trueeffects, by.x ="leftout", by.y = "group")

res2$trueeffects$group <- as.numeric(factor(res2$trueeffects$J))
res2$trueeffects$cluster <- factor(getAvgCluster(res2))

tmp2 <- merge(res2$leaveouts, res2$alleffects, by.x = "leftout", by.y = "group")
tmp21 <- merge(res2$leaveouts, res2$trueeffects, by.x ="leftout", by.y = "group")

ests <- res1$leaveouts |> group_by(leftout) |> summarize(dpsur = median(dpsur),
                                                         null = median(null),
                                                         simple = median(simple))
ests <- merge(ests, res1$trueeffects, by.x = "leftout", by.y = "group")

ests2 <- res2$leaveouts |> group_by(leftout) |> summarize(dpsur = median(dpsur),
                                                         null = median(null),
                                                         simple = median(simple))
ests2 <- merge(ests2, res2$trueeffects, by.x = "leftout", by.y = "group")


p1 <- ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = dpsur)) +
  geom_density2d(data = tmp1, aes(x = seff, y = dpsur, color = cluster)) +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = dpsur, color = cluster)) +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none", color = "none") + ggtitle(name1, subtitle = "DPM (clusters denoted by colors)")


p2 <- ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = null)) +
  geom_density2d(data = tmp1, aes(x = seff, y = null), color = "grey65") +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = null), color = "grey55") +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("",subtitle = "NULL")



p3 <- ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = simple)) +
  geom_density2d(data = tmp1, aes(x = seff, y = simple), color = "grey65") +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = simple), color = "grey55") +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("", subtitle = "Simple")


p12 <- ggplot(ests2, aes(x = seff, y = yeff, xend = seff, yend = dpsur)) +
  geom_density2d(data = tmp21, aes(x = seff, y = dpsur, color = cluster)) +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = dpsur, color = cluster)) +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none", color = "none") + ggtitle(name2, subtitle = "DPM (clusters denoted as colors)")


p22 <- ggplot(ests2, aes(x = seff, y = yeff, xend = seff, yend = null)) +
  geom_density2d(data = tmp21, aes(x = seff, y = null), color = "grey65") +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = null), color = "grey55") +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("", subtitle = "NULL")



p32 <- ggplot(ests2, aes(x = seff, y = yeff, xend = seff, yend = simple)) +
  geom_density2d(data = tmp21, aes(x = seff, y = simple), color = "grey65") +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = simple), color = "grey55") +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("", subtitle = "Simple")

p1 + p12 +  p3 + p32 + p2 + p22 +plot_layout(ncol = 2)


ggsave("preds-median-twotrt.png", width = 7.5, height = 9.75)


cheez <-  merge(res1$leaveouts, res1$trueeffects, by.x = "leftout", by.y = "group")

cheez$D <- abs(cheez$dpsur - cheez$yeff)
cheez$D0 <- abs(cheez$null - cheez$yeff)
cheez$Ds <- abs(cheez$simple - cheez$yeff)

cheez2 <-  merge(res2$leaveouts, res2$trueeffects, by.x = "leftout", by.y = "group")

cheez2$D <- abs(cheez2$dpsur - cheez2$yeff)
cheez2$D0 <- abs(cheez2$null - cheez2$yeff)
cheez2$Ds <- abs(cheez2$simple - cheez2$yeff)

q1 <- data.frame(err = c(cheez$D, cheez$D0, cheez$Ds),
           type = rep(c("DPM", "null", "simple"), each = nrow(cheez))) |>
  ggplot(aes(x = err, color = type)) + geom_density() + theme_bw() + xlab("Absolute error") +
  xlim(c(0, 10)) + ggtitle(name1)


q2 <- data.frame(err = c(cheez2$D, cheez2$D0, cheez2$Ds),
                 type = rep(c("DPM", "null", "simple"), each = nrow(cheez2))) |>
  ggplot(aes(x = err, color = type)) + geom_density() + theme_bw() + xlab("Absolute error") +
  xlim(c(0,10)) + ggtitle(name2)

q1 + q2 + plot_layout(ncol = 1, guides = "collect")
ggsave("D-density-twotrt.png", width = 6.5, height = 6.5)

cheez2$trt <- sapply(strsplit(cheez2$J, ":"), "[", 2)
cheez2$genome <- sapply(strsplit(cheez2$J, ":"), "[", 1)

ggplot(data.frame(D = c(cheez2$D, cheez2$Ds, cheez2$D0),
                  model = rep(c("DPM", "Simple", "Null"), each = nrow(cheez2)),
                  treatment = c(cheez2$trt, cheez2$trt, cheez2$trt)),
       aes(x = D, y = treatment, color = model)) + geom_boxplot() + theme_bw()

ggsave("onetrt-box.png", width = 5.25, height = 4.5)


tmp3 <- merge(res1$alleffects, res1$trueeffects, by = "group")

ests <- res1$alleffects |> group_by(group) |> summarize(dpsur = median(dpsur.yeff))
ests <- merge(ests, res1$trueeffects, by.x = "group", by.y = "group")

ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = dpsur)) +
  geom_density2d(aes(x = seff, y = dpsur), color = "grey65") +
  geom_segment(color = "grey85") +
  geom_point() +
  geom_point(aes(x = seff, y = dpsur), color = "grey55") +
  theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("Nonlinear", subtitle = "DPM")




ggplot( rbind(cbind(cheez, setting = name1),
              cbind(cheez2, setting = name2)),
             aes(x = D, color = factor(cluster))) + geom_density() +
  theme_bw() + xlab("Absolute error") + facet_wrap(~ setting, ncol = 1) +
  ggtitle("by cluster") + guides(color = "none") + xlim(c(0,10))

ggsave("d-dens-by-cluster-twotrt.png", width = 5.25, height = 4.25)

