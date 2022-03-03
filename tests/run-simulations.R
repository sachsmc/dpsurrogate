library(dpsurrogate)


for(i in 1:1) {


res <- run_one_replicate(effect = "nonlinear", Zeffect = FALSE,
                        Ueffect = FALSE, niter = 50, nnodes = 7)



saveRDS(res, file = file.path(sprintf("tests/sim-res/res-nonlin-take2-0-0-%.3d.rds", i)))


res <- run_one_replicate(effect = "nonlinear", Zeffect = TRUE,
                         Ueffect = FALSE, niter = 50, nnodes = 8)



saveRDS(res, file = file.path(sprintf("tests/sim-res/res-nonlin-take2-1-0-%.3d.rds", i)))

for(i in 2:10) {

  cat(i, "\n")
  res <- run_one_replicate(effect = "nonlinear", Zeffect = TRUE,
                         Ueffect = TRUE, niter = 50, nnodes = 4)


saveRDS(res, file = file.path(sprintf("tests/sim-res/res-nonlin-take2-1-1-%.3d.rds", i)))

}




}
