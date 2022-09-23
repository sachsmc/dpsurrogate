args <- commandArgs(trailingOnly = TRUE)
runi <- args[1]
effect <- args[2]
Zeff <- as.numeric(args[3])
Ueff <- as.numeric(args[4])
pid <- as.numeric(args[5])

library(dpsurrogate)

niter <- 50

if(effect == "nonlinearskew") {

  library(sn)
  idat <- generate_data("nonlinear", Zeff, Ueff,
                        Zgen = function(n) rsn(n, xi = 0, omega = 1, alpha = 1),
                        Sgen = function(n) rsn(n, xi = 0, omega = 1.5, alpha = 1))

} else {

  idat <- generate_data(effect, Zeff, Ueff)
}


resall0 <- run_one_analysis(idat, niter = 50)
resnull <- run_one_null(idat, niter = 50)
ressimp <- run_one_simple_analysis(idat, niter = 50)

resall <- data.frame(dpsur.seff = c(resall0$res.s),
             dpsur.yeff = c(resall0$res.y),
             #null.seff = c(resnull$res.s),
             null.yeff = c(resnull$res.y),
             simple.seff = c(ressimp$res.s),
             simple.yeff = c(ressimp$res.y),
             group = rep(1:64, each = niter))

runstr <- sprintf("%s-%.1f-%.1f_%.2d", effect,
		  Zeff, Ueff, as.numeric(runi))
dir.create(paste0("tmp", pid, "-", runstr))
setwd(paste0("tmp", pid, "-", runstr))
saveRDS(idat, paste0("data-", runstr, ".rds"))
saveRDS(resall, paste0("resall-", runstr, ".rds"))
saveRDS(resall0, paste0("resall0-", runstr, ".rds"))

