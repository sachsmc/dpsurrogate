library(dpsurrogate)


for(i in 5:10) {


res <- run_one_replicate(effect = "nonlinear", Zeffect = FALSE,
                        Ueffect = FALSE, niter = 50, nnodes = 6)



saveRDS(res, file = file.path(sprintf("tests/sim-res/res-nonlin-0-0-%.3d.rds", i)))

}
