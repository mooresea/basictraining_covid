##===================================##
## Author: Guido España
##===================================##
## Libraries-------------
##===================================##
library(tidyverse)
output_dir = '../output_experiments_joint'
output_dir_sims = '../output_experiments_joint'

##===================================##
## Gather data-------------
##===================================##
parameters_df = read_csv(file.path(output_dir_sims,'params_file.csv'))
parameters_df$Finished = 0
combined_output_file_t = file.path(output_dir, "calibration_sim_training.csv")
combined_output_file_q = file.path(output_dir, "calibration_sim_quarantine.csv")
combined_output_df_t = combined_output_df_q = data.frame()
sims_list = sort(unique(parameters_df$simID))
## For each parameter, determine if it finished or not
for(ii in sims_list){
    ## Check if simulation finished
    out_file_t = file.path(output_dir_sims,sprintf("calibration_sim_%d_training.csv",ii))
    out_file_q = file.path(output_dir_sims,sprintf("calibration_sim_%d_quarantine.csv",ii))
    if(file.exists(out_file_t)){
        output_df_t = read.csv(out_file_t, stringsAsFactors = F)
        output_df_t$simID=ii
        combined_output_df_t = bind_rows(combined_output_df_t, output_df_t)
        output_df_q = read.csv(out_file_q, stringsAsFactors = F)
        output_df_q$simID=ii
        combined_output_df_q = bind_rows(combined_output_df_q, output_df_q)
    }else{
        print(sprintf("File %s doesn't exist", out_file_t))
    }
}
write_csv(combined_output_df_t, combined_output_file_t)
write_csv(combined_output_df_q, combined_output_file_q)

