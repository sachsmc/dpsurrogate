args <- commandArgs(trailingOnly = TRUE)
runi <- args[1]
effect <- args[2]
Zeff <- as.numeric(args[3])
Ueff <- as.numeric(args[4])
pid <- as.numeric(args[5])

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

options(error = traceback)

library(dpsurrogate)

niter <- 50

runstr <- sprintf("%s-%.1f-%.1f_%.2d", effect,
                  Zeff, Ueff, as.numeric(runi))

setwd(paste0("tmp",pid,"-", runstr))
idat <- readRDS(paste0("data-", runstr, ".rds"))

restmp <- run_one_loo(idat, task_id, niter = niter)

restmp2 <- tryCatch(run_one_loo_null(idat, task_id, niter = niter),
                    error = function(e) "error")
if(identical(restmp2, "error")) {
  rerun <- 1
  while(identical(restmp2, "error")) {

    restmp2 <- tryCatch(run_one_loo_null(idat, task_id, niter = niter),
                        error = function(e) "error")

    if(rerun > 10) stop("You done goofed")
    rerun <- rerun + 1
  }
}

restmp3 <- run_one_loo_simple_analysis(idat, task_id, niter = niter)
oresi <- data.frame(dpsur = restmp$res.y, null = restmp2$res.y,
                    simple = restmp3$res.y, leftout = task_id)

saveRDS(oresi, paste0("oresi-", task_id, ".rds"))
saveRDS(restmp$clusters, paste0("cluster-", task_id, ".rds"))

