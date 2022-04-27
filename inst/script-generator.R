basescpt <- '#!/bin/bash

module load gcc/11.2.0
module load R/4.1.2
module load jags/4.3.0

SCEN=%s

for i in {0..%.0f}
do
srun -n1 Rscript setup-generate.R $i $SCEN %.1f %.1f
RES=$(sbatch --parsable --array=1-64 --wrap="Rscript oneout.R $i $SCEN %.1f %.1f")
echo $RES
srun --dependency=afterok:$RES -n1 Rscript cleanup.R $i $SCEN %.1f %.1f
printf -v j "%%02d" $i
rm -r tmp$SCEN-%.1f-%.1f_$j/
rm *.out
done
'

makescpt <- function(scen, num, p1, p2, file) {
  myscpt <- sprintf(basescpt, scen, num, p1, p2, p1, p2, p1, p2, p1, p2)
  cat(myscpt, file = file)
}

settings <- rbind(data.frame(scen = c("nonlinear", "linear", "simple", "null", "inter", "interhide",
                    "onetrt", "twotrt", "manybiom"),
           num = 99, p1 = 0, p2 = 0),
      data.frame(scen = c("nonlinear", "linear", "simple"),
                 num = 99, p1 = 0.3, p2 = 0),
      data.frame(scen = c("nonlinear", "linear", "simple"),
                 num = 99, p1 = 0.3, p2 = 0.3),
      data.frame(scen = c("inter"),
                 num = 99, p1 = 0, p2 = 0.3),
      data.frame(scen = c("interhide"),
                 num = 99, p1 = 0.3, p2 = 0))

settings$fname <- gsub(".", "", make.names(settings$scen, unique = TRUE), fixed = TRUE)

for(i in 1:nrow(settings)){
  makescpt(settings$scen[i], settings$num[i], settings$p1[i],
           settings$p2[i], paste0(settings$fname[i], ".sh"))
}
