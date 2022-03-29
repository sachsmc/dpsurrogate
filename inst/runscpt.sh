#!/bin/bash

module load R/4.1.2

srun -n1 Rscript -e "x <- rnorm(10); saveRDS(x, file = 'test1.rds')"
sbatch --array=1-10 --wrap="Rscript oneout.R"
