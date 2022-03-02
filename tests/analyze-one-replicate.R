library(dplyr)
library(ggplot2)
library(patchwork)
library(dpsurrogate)

res1 <- readRDS("tests/sim-res/res-nonlin-50-1-0-003.rds")
res1$trueeffects$group <- as.numeric(factor(res1$trueeffects$J))
res1$leaveouts <- do.call(rbind, res1$leaveouts)

tmp <- merge(res1$leaveouts, res1$alleffects, by.x = "leftout", by.y = "group")
tmp1 <- merge(res1$leaveouts, res1$trueeffects, by.x ="leftout", by.y = "group")

summary(with(tmp, abs(dpsur - dpsur.yeff)))
summary(with(tmp, abs(null - dpsur.yeff)))
summary(with(tmp, abs(simple - dpsur.yeff)))


summary(with(tmp1, abs(dpsur - yeff)))
summary(with(tmp1, abs(null - yeff)))
summary(with(tmp1, abs(simple - yeff)))

ests <- res1$leaveouts |> group_by(leftout) |> summarize(dpsur = median(dpsur),
                                                         null = median(null),
                                                         simple = median(simple))
ests <- merge(ests, res1$trueeffects, by.x = "leftout", by.y = "group")


p1 <- ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = dpsur)) +
  geom_segment() +
  geom_point(shape = 1) +
  geom_point(aes(x = seff, y = dpsur)) + theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("Leave one outs medians (filled) vs full (open) - DPM") +
  ylim(c(-3.5, 3))


p2 <- ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = null)) +
  geom_segment() +
  geom_point(shape = 1) +
  geom_point(aes(x = seff, y = null)) + theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("Leave one outs median (filled) vs full (open) - NULL")+
  ylim(c(-3.5, 3))



p3 <- ggplot(ests, aes(x = seff, y = yeff, xend = seff, yend = simple)) +
  geom_segment() +
  geom_point(shape = 1) +
  geom_point(aes(x = seff, y = simple)) + theme_bw() +
  xlab(expression(hat(nu))) + ylab(expression(hat(mu))) +
  guides(alpha = "none") + ggtitle("Leave one outs median (filled) vs full (open) - Simple")+
  ylim(c(-3.5, 3))

p1 + p2 + p3 + plot_layout(ncol = 1)


ggsave("preds-median.png", width = 7.5, height = 9.75)


p1 <- ggplot(ests, aes(x = ngrp, y = dpsur - yeff)) + geom_point() +
  ylim(c(-3, 3)) + xlab("Group sample size") + ylab("Error") + ggtitle("DPM") + theme_bw()
p2 <- ggplot(ests, aes(x = ngrp, y = simple - yeff)) + geom_point() +
  ylim(c(-3, 3)) +  xlab("Group sample size") + ylab("Error") + ggtitle("Simple") + theme_bw()

p1 + p2 + plot_layout(ncol = 2)
ggsave("preds-byN.png", width = 7.5, height = 4.5)



cheez2 <-  merge(res1$leaveouts, res1$trueeffects, by.x = "leftout", by.y = "group")

cheez2$D <- abs(cheez2$dpsur - cheez2$yeff)
cheez2$D0 <- abs(cheez2$null - cheez2$yeff)
cheez2$Ds <- abs(cheez2$simple - cheez2$yeff)

data.frame(err = c(cheez2$D, cheez2$D0, cheez2$Ds),
           type = rep(c("DPM", "null", "simple"), each = nrow(cheez2))) |>
  ggplot(aes(x = err, color = type)) + geom_density() + theme_bw() + xlab("Absolute error")
ggsave("D-density.png", width = 6.5, height = 4.5)

p1 <- ggplot(cheez2, aes(x = seff, y = dpsur)) + geom_density2d() +
  geom_point(data = res1$trueeffects, aes(x = seff, y = yeff)) + theme_bw() +
  ggtitle("DPM") + ylim(c(-3.5, 3)) + xlab(expression(hat(nu))) + ylab(expression(hat(mu)))

p2 <- ggplot(cheez2, aes(x = seff, y = null)) + geom_density2d() +
  geom_point(data = res1$trueeffects, aes(x = seff, y = yeff)) + theme_bw() +
  ggtitle("NULL") + ylim(c(-3.5, 3))+ xlab(expression(hat(nu))) + ylab(expression(hat(mu)))


p3 <-ggplot(cheez2, aes(x = seff, y = simple)) + geom_density2d() +
  geom_point(data = res1$trueeffects, aes(x = seff, y = yeff)) + theme_bw() +
  ggtitle("Simple") + ylim(c(-3.5, 3))+ xlab(expression(hat(nu))) + ylab(expression(hat(mu)))


p1+p2+p3 + plot_layout(ncol = 1)
ggsave("pred-density.png", width = 7.5, height = 9.75)


rmarkdown::render("tests/summarize-replicates.Rmd",
                  params = list(path = "sim-res/res-nonlin-take2-0-0-001.rds"),
                  output_file = "noZnoU.html")

rmarkdown::render("tests/summarize-replicates.Rmd",
                  params = list(path = "sim-res/res-nonlin-take2-1-0-001.rds"),
                  output_file = "yesZnoU.html")

rmarkdown::render("tests/summarize-replicates.Rmd",
                  params = list(path = "sim-res/res-nonlin-take2-1-1-001.rds"),
                  output_file = "yesZyesU.html")
