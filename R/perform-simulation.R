#' Run one complete simulation
#'
#' @param scenario
#' @export

run_one_replicate <- function(effect = "nonlinear", Zeffect = FALSE,
                              Ueffect = FALSE, nnodes = 8,
                              niter = 50) {

  idat <- generate_data(effect, Zeffect, Ueffect)

  resall <- run_one_analysis(idat, niter = niter)
  resnull <- run_one_null(idat, niter = niter)
  ressimp <- run_one_simple_analysis(idat, niter = niter)

  resall <- data.frame(dpsur.seff = c(resall$res.s),
             dpsur.yeff = c(resall$res.y),
             #null.seff = c(resnull$res.s),
             null.yeff = c(resnull$res.y),
             simple.seff = c(ressimp$res.s),
             simple.yeff = c(ressimp$res.y),
             group = rep(1:64, each = niter))

  cl <- makeCluster(nnodes)

  init <- clusterEvalQ(cl, {
    library(dpsurrogate)
  })
  clusterExport(cl, c("idat"), envir = environment())
#for(i in 1:64) {
 resWho <- clusterApply(cl, 1:64, function(i) {
    restmp <- run_one_loo(idat, i, niter = niter, jags.state = resall$jags.state)
    restmp2 <- run_one_loo_null(idat, i, niter = niter, jags.state = resnull$jags.state)
    restmp3 <- run_one_loo_simple_analysis(idat, i, niter = niter, jags.state = ressimp$jags.state)
    data.frame(dpsur = restmp$res.y, null = restmp2$res.y,
               simple = restmp3$res.y, leftout = i)
    cat(i, "\n")
#}
 })
  resWho <- do.call(rbind, resWho)

  stopCluster(cl)

  ## organize and save samples from loos to compute prediction error
  ## also save the original data generated
  list(trueeffects = idat$rdat,
       alleffects = resall,
       leaveouts = resWho,
       setting = match.call())

}
