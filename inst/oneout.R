task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
x <- readRDS("test1.rds")
mean(x[-task_id])
