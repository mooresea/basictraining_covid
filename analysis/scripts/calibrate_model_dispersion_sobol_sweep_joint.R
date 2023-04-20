##==================================================#
## load function to simulate outbreaks--------------
##==================================================#
library(reshape2)
library(scam)
library(mgcv)
library(MASS)
library(tidyverse)
library(rriskDistributions)
library(foreach)
library(doParallel)
library(parallel)
print("Start 1")
source('simOutbreak_data_betabinomR0_vaccination_updated2.R')
nreps = 200
sim_ID = -1
output_dir = '../output'
args = commandArgs(TRUE)

if(length(args) < 1){
    set.seed(123456)
}
if(length(args) >= 1){
    output_dir = '../output_experiments_joint_loAsymp'
    sim_ID = as.numeric(args[1])
    set.seed(runif(1,min = 1, max = as.integer(Sys.time())))    
}

##Check to see if file exists - if so don't need to rerun
ofiles=list.files("../output_experiments_joint_loAsymp")
if(any(grepl(paste0("sim_",as.character(sim_ID),"_"),ofiles))){
	#File already exists - exit
 print("File already exists")
	quit("no") 
}else{



##Recruit data
wr_dat=read_csv(file="../output/estimated_recruits_testing_positive.csv")
wr_dat=wr_dat[which(!is.na(wr_dat$Battalion_fill)&!is.na(wr_dat$Training_infected)&wr_dat$Battalion_fill>0),]
##Infection estimates
inf_est=read_csv(file="../output/daily_estimated_infections_recruits_aggregated_01-11-22.csv")
##Infection estimates in Fort Jackson area (Richland County, SC)
sc_dat=read_csv(file="../output/daily_estimated_infections_Fort_Jackson_01-11-22.csv")
##Vaccination estimates
vacc_dat=read_csv(file="../data/military_recruit_vaccination_timeseries.csv")
vacc_dat_sc=read_csv(file="../data/Richland_SC_vaccination_timeseries.csv")
importations_list = rep(NA,nrow(wr_dat))

##==================================================#
## run sims across different sites----------
##==================================================#
if(sim_ID <= 0){
    R0_sweep = seq(from = 2, to = 30, length.out = 1e3)
    disp_sweep=10^(seq(-5,5,length.out = 0.5e3))
    #calib_params_df = read_csv('../output/calibrated_parameters_forecast.csv')
    # FB_init_sweep = seq(from = calib_params_df$initial_inf_low[calib_params_df$base_name == 'FB'],
    #                        to = calib_params_df$initial_inf_low[calib_params_df$base_name == 'FB'], length.out = 2e5)
    # FLW_init_sweep = seq(from = calib_params_df$initial_inf_low[calib_params_df$base_name == 'FLW'],
    #                             to = calib_params_df$initial_inf_low[calib_params_df$base_name == 'FLW'], length.out = 2e5)

}else{
    R0_df = read.csv('../output_experiments_joint/params_file.csv')
    R0_sweep = R0_df$R0[R0_df$simID == sim_ID]
    disp_sweep = R0_df$dispersion[R0_df$simID == sim_ID]
    adj_sweep=R0_df$R0_adj[R0_df$simID == sim_ID]
    uids=R0_df$uID[R0_df$simID == sim_ID]
}
time_max=73
testpos.time = expand.grid(uID = R0_df$uID[R0_df$simID==sim_ID], rep=1:nreps,stringsAsFactors = F)
testpos.time$Ind=1:length(R0_sweep)
testpos.time$R0_mag = R0_sweep[testpos.time$Ind]
testpos.time$dispersion = disp_sweep[testpos.time$Ind]
testpos.time$R0_adj = adj_sweep[testpos.time$Ind]
wr_weeks=wr_dat$Week
cnt=1
for(i in ncol(testpos.time):(ncol(testpos.time)+length(wr_weeks)-1)) { 
    testpos.time[,i+1]=0
    names(testpos.time)[i+1]=wr_weeks[cnt]
    cnt=cnt+1
}
testpos.time2=testpos.time

#cl <- makeCluster(16,outfile="",setup_strategy="sequential")  
#registerDoParallel(cl) 

for(nn_r in 1:length(R0_sweep)){
    if(nn_r %% 50 == 0){print(nn_r)}
    for(ii_param in 1:length(importations_list)){ 
        print(ii_param)
        source('paramDefaults.R')
        start_date=wr_dat$Start_date[ii_param]
        inf_pp=inf_est[which(inf_est$date<=start_date),]
        inf_pp=inf_pp[rev(order(inf_pp$date)),]
        inf_sc=sc_dat[which(sc_dat$date<=start_date),]
        inf_sc=inf_sc[rev(order(inf_sc$date)),]
        
        ##Estimate beta distribution for daily infection probability
        inf_pp$inf_logmn=NA
        inf_pp$inf_logsd=NA
        for(ii in 1:nrow(inf_pp)){
            nii=get.norm.par(p=c(0.025,.5,.975),q=log(c(inf_pp$infections_lo_pct[ii],inf_pp$infections_pct[ii],inf_pp$infections_hi_pct[ii])),show.output = F,plot=F)
            inf_pp$inf_logmn[ii]=nii[1]
            inf_pp$inf_logsd[ii]=nii[2]
        }
        
        ##Need to create a new contact network for each week b/c recruit numbers differ
        numrecruits=wr_dat$Battalion_fill[ii_param]
        source("network_BT_data.R")
        time = 1:time_max
        R0.prob = R0_sweep[nn_r] / mean(average.contacts) / (symp.mag + (1 - symp.mag) * asymp.adjust)
        arrivebase = c(sample(1:1, numrecruits, replace = T), rep(0, numdrillserg + numStaff))
        #print(length(arrivebase))
        # ## Importations
        ## Importations - daily importation rate over course of training camp
        imp_dates=which(sc_dat$date==start_date):which(sc_dat$date==(start_date+time_max-1))
        imp_rate=sc_dat$pct_inf[imp_dates]
        #importation = c(rep(0, numrecruits), rep(imp_rate, numdrillserg + numStaff))
        importation_mat = matrix(0,nrow=time_max,ncol=( numrecruits+ numdrillserg + numStaff))
        for(im in (numrecruits+1):ncol(importation_mat)){
          importation_mat[,im]=imp_rate
        }
        #Change arrival information based on date
        initial.immune=sum(inf_pp$infections_pct)
        initial.infected=sum(inf_pp$infections_pct)
        initial.staff.immune=initial.staff.infected=sum(sc_dat$pct_inf[sc_dat$date<=start_date])
        
        BOOL_testDaily = FALSE
        BOOL_testOnArrival = 1
        testdates = c(14)
        quarantine.contacts = 0

        ## Assume mask compliance = 0
        #mask.compliance = 0.2
        compliance.avg = mask.protection * mask.compliance
        compliance = rep(compliance.avg, numrecruits + numdrillserg + numStaff)
        
        ## simulate
        disp=disp_sweep[nn_r]
        num.variants=5
        R0_adj_factor=adj_sweep[nn_r]
        #s = replicate(nreps,simOutbreak(inf_pp,inf_sc,vacc_dat,vacc_dat_sc))
        sim_out_l=list()
        sim_out_l=foreach(iter=1:nreps,.packages = 'rmutil') %do% {
	          .GlobalEnv$R0_adj_factor <- R0_adj_factor
	          sim_o<-simOutbreak(inf_pp,inf_sc,vacc_dat,vacc_dat_sc,sampleR=F)
	          return(sim_o)
        }
        
        ##Switch to counting infections during quarantine period including symptomatic testing
        pos.quarantine=unlist(lapply(sim_out_l, function(x) sum((x[[3]]$init.testpos.dat>1&x[[3]]$init.testpos.dat<=14),na.rm = T)))
        pos.training=unlist(lapply(sim_out_l, function(x) sum(x[[3]]$init.testpos.dat>14,na.rm = T)))
        
        testpos.time[which(testpos.time$uID == uids[nn_r]),wr_weeks[ii_param]] = pos.quarantine         
        testpos.time2[which(testpos.time$uID == uids[nn_r]),wr_weeks[ii_param]] = pos.training  
    }
}
#stopCluster(cl)
file_outsims = ifelse(sim_ID != -1, sprintf('calibration_sim_%.0f_quarantine.csv',sim_ID),'calibration_sim_quarantine.csv')
write.csv(testpos.time, file.path(output_dir,file_outsims), row.names = F)
file_outsims2 = ifelse(sim_ID != -1, sprintf('calibration_sim_%.0f_training.csv',sim_ID),'calibration_sim_training.csv')
write.csv(testpos.time2, file.path(output_dir,file_outsims2), row.names = F)
}
