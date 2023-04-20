##=============================================================================##
## Load libraries----------------------
##=============================================================================##
library(tidyverse)
library(pomp)
set.seed(123456)

##=============================================================================##
## Split simulations----------------------
##=============================================================================##
reps = 1000
simulation_IDs = 1:reps
output_dir = '../output_experiments_joint'
## Remove previous simulation outputs
if(!dir.exists(output_dir)){
    dir.create(output_dir)
}
calib_params_df = read_csv('../output/calibrated_parameters_forecast.csv')

system(sprintf('rm %s/calibration_sim_*.csv', output_dir,sep = ''))

sobol_df=sobol_design(lower=c(R0=1,disp_log=-6),upper=c(R0=30,disp_log=4),nseq=4000)
sobol_df$dispersion=10^(sobol_df$disp_log)
sobol_df$R0_adj=1 #calib_params_df$dispersion_mag
sobol_df$simID = sort(rep(1:reps,length.out = nrow(sobol_df)))
sobol_df$uID=1:nrow(sobol_df)

#params_df = data.frame(simID = sim_ids, R0 = R0_sweep, , stringsAsFactors = F)
write_csv(sobol_df, '../output_experiments_joint/params_file.csv')

##=============================================================================##
## Submit to CRC array---------------------------
##=============================================================================##
cores_total = 8
jobname = "COVID_DOD_R0_FIT"
submission_template = "#!/bin/csh
#$ -q long
#$ -t 1:JOBSQUEUE:1
#$ -N JOBNAME
#$ -pe smp JOBCORES

Rscript calibrate_model_dispersion_sobol_sweep_joint.R ${SGE_TASK_ID}
"    
submission_str = submission_template %>%
    str_replace_all(pattern="JOBNAME", replacement = jobname) %>%
    str_replace_all(pattern="JOBSQUEUE", replacement = as.character(reps)) %>%
    str_replace_all(pattern="JOBCORES", replacement = as.character(cores_total))

submission_file = sprintf("%s.sh",jobname)
file.connection = file(submission_file)
write(submission_str,file.connection)
close(file.connection)

system(sprintf("qsub %s", submission_file))
