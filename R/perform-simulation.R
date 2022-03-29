#' Run one complete simulation
#'
#' @param scenario
#' @export

run_one_replicate <- function(effect = "nonlinear", Zeffect = 0,
                              Ueffect = 0,
                              niter = 50, nnodes) {

  idat <- generate_data(effect, Zeffect, Ueffect)

  resall0 <- run_one_analysis(idat, niter = niter)
  resnull <- run_one_null(idat, niter = niter)
  ressimp <- run_one_simple_analysis(idat, niter = niter)

  resall <- data.frame(dpsur.seff = c(resall0$res.s),
             dpsur.yeff = c(resall0$res.y),
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

  resWho <- clusterApply(cl, 1:64, function(i) {

 #resWho <- lapply(cl, 1:64, function(i) {
   tryCatch({
    restmp <- run_one_loo(idat, i, niter = niter, jags.state = resall0$jags.state)
    restmp2 <- run_one_loo_null(idat, i, niter = niter, jags.state = resnull$jags.state)
    restmp3 <- run_one_loo_simple_analysis(idat, i, niter = niter, jags.state = ressimp$jags.state)
    data.frame(dpsur = restmp$res.y, null = restmp2$res.y,
               simple = restmp3$res.y, leftout = i)
   }, error = function(e){

     list(error = e,
          lout = i,
          dset = idat
          )

   })
    #cat(i, "\n")
 })

  check.errors <- sapply(resWho$leaveouts, function(x) !is.null(x$error))

  if( any(check.errors)){

    for(j in which(check.errors)) {
      idat <- resWho$leaveouts[[j]]$dset
      restmp <- run_one_loo(idat, j, niter = niter, jags.state = NULL)
      restmp2 <- run_one_loo_null(idat, j, niter = niter, jags.state = NULL)
      restmp3 <- run_one_loo_simple_analysis(idat, j, niter = niter, jags.state = NULL)
      resWho$leaveouts[[j]] <- data.frame(dpsur = restmp$res.y, null = restmp2$res.y,
                                        simple = restmp3$res.y, leftout = j)

    }

  }
  #resWho <- do.call(rbind, resWho)

 stopCluster(cl)

  ## organize and save samples from loos to compute prediction error
  ## also save the original data generated
  list(trueeffects = idat$rdat,
       alleffects = resall,
       leaveouts = resWho,
       setting = match.call(),
       clusters = resall0$clusters)

}
