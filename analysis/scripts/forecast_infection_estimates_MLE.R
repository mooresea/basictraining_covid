##==================================================#
## load function to simulate outbreaks--------------
##==================================================#
library(reshape2)
library(tidyverse)
library(rriskDistributions)
library(foreach)
library(doParallel)
library(parallel)
source('simOutbreak_data_betabinomR0_vaccination_updated2.R')
training_site="full"

hiAsymp=T
source('paramDefaults.R')
if(hiAsymp){
  symp.mag=0.31
}else{
  symp.mag=0.57  
}

set.seed(22)

nreps = 200 #1e3

##Recruit data
wr_dat=read_csv(file="../output/estimated_recruits_testing_positive.csv")
wr_datA=wr_dat
wr_dat=wr_dat[which(!is.na(wr_dat$Battalion_fill)&!is.na(wr_dat$Quarantine_positives)&wr_dat$Battalion_fill>0),]
##Infection estimates
inf_est=read_csv(file="../output/daily_estimated_infections_recruits_aggregated_4wk_forecast.csv")
##Infection estimates in Fort Jackson area (Richland County, SC)
sc_dat=read_csv(file="../output/daily_estimated_infections_Fort_Jackson.csv")
##Vaccination estimates
vacc_dat=read_csv(file="../data/military_recruit_vaccination_timeseries.csv")
vacc_dat_sc=read_csv(file="../data/Richland_SC_vaccination_timeseries.csv")

#R0 posterior samples
if(hiAsymp){
  R0_samples=readRDS(file='../output/R0_adjust_posterior_samples_dispersion_hiAsymp.rds')
  
}else{
  R0_samples=readRDS(file='../output/R0_adjust_posterior_samples_dispersion_loAsymp.rds')

}

##Extend vaccination time series assuming that it remains constant to end of infection time series
vdates_add=as.Date((max(vacc_dat$Date)+1):max(inf_est$date),origin="1970-01-01")
vdates_add_sc=as.Date((max(vacc_dat_sc$Date)+1):max(inf_est$date),origin="1970-01-01")
vacc_dat[(nrow(vacc_dat)+1):(nrow(vacc_dat)+length(vdates_add)),]<-NA
vacc_dat$Date[(nrow(vacc_dat)-length(vdates_add)+1):nrow(vacc_dat)]=vdates_add
vacc_dat[(nrow(vacc_dat)-length(vdates_add)+1):nrow(vacc_dat),-1]=vacc_dat[(nrow(vacc_dat)-length(vdates_add)),-1]
vacc_dat_sc[(nrow(vacc_dat_sc)+1):(nrow(vacc_dat_sc)+length(vdates_add_sc)),]<-NA
vacc_dat_sc$Date[(nrow(vacc_dat_sc)-length(vdates_add_sc)+1):nrow(vacc_dat_sc)]=vdates_add_sc
vacc_dat_sc[(nrow(vacc_dat_sc)-length(vdates_add_sc)+1):nrow(vacc_dat_sc),-1]=vacc_dat_sc[(nrow(vacc_dat_sc)-length(vdates_add_sc)),-1]

##==================================================#
## Functions----------
##==================================================#
error.bar <- function(x, y, upper, lower, length=0.1,...){
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}

##
## Extend simulations
##
new_dates=seq(wr_dat$Start_date[nrow(wr_dat)]+7,max(inf_est$date),by=7)
wr_dat[(nrow(wr_dat)+1):(nrow(wr_dat)+length(new_dates)),]<-NA
wr_dat$Start_date[(nrow(wr_dat)-length(new_dates)+1):nrow(wr_dat)]=new_dates
wr_dat$Arrivals[(nrow(wr_dat)-length(new_dates)+1):nrow(wr_dat)]=1200
wr_dat$Battalion_fill[(nrow(wr_dat)-length(new_dates)+1):nrow(wr_dat)]=1200
wr_dat$Total_size=NA
wr_dat$pct_inf_prior_week=NA
wr_dat$pct_inf_prior_2week=NA
##==================================================#
## run sims across different immunity levels----------
##==================================================#
importations_list = rep(NA,nrow(wr_dat)) #c(0,1,10)

testdates_list = list(c(),c(14),c(),c(14),c(14),c(),c(14),
                      c(14),c(14),c(14),c(14),c(14),c(),c(),c(),c(),c(14),c(14),c(),c(14),c(14))
testarrival_list = c(0,0,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,1,0,1,1)
test_type = list(c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),
                 c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'),c('pcr'))
test_scenarios = c('None','Quarantine','Arrival','Both','No Masks','No Masks or Tests',
                   'High Masks','No Vaccine','High Vaccine', 'Staff Vaccine','Early Vaccine',
                   'Staff Tests','Staff Tests Only','Early Vaccine Only','High Vaccine Only','Staff Vaccine Only',
                   'Staff Tests Daily Antigen','Staff Tests Weekly Antigen','Staff Tests Daily Antigen Only',
                   'Early Vaccine and Staff Tests Daily Antigen','Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')

time_max =73

init.testDay = inf.testDay = matrix(NA,nreps * length(importations_list),length(test_scenarios) + 2)
init.testDay[,length(test_scenarios) + 2] = rep(1:length(importations_list), each = nreps)
init.testDay[,length(test_scenarios) + 1] = rep(1:nreps, length(importations_list))

inf.testDay[,length(test_scenarios) + 2] = rep(1:length(importations_list), each = nreps)
inf.testDay[,length(test_scenarios) + 1] = rep(1:nreps, length(importations_list))

infA.testDay=infA2.testDay=infQ.testDay=infT.testDay=sympQ.testDay=sympT.testDay=imp.testDay=pos.testDay=posT.testDay=inf.testDay

inf.time = matrix(NA,nreps*time_max*length(importations_list),length(test_scenarios) + 2)
inf.time[,length(test_scenarios) + 1] = rep(rep(1:nreps, each = time_max), length(importations_list))
inf.time[,length(test_scenarios) + 2] = rep(1:length(importations_list), each = time_max * nreps)

cl <- makeCluster(10)  
registerDoParallel(cl) 

compliance_init=compliance
for(pp_import in 1:nrow(wr_dat)){
    #print(sprintf("Importaion multiplier = %.0f", importations_list[pp_import]))

    start_date=wr_dat$Start_date[pp_import]
    inf_pp=inf_est[which(inf_est$date<=start_date),]
    inf_pp=inf_pp[rev(order(inf_pp$date)),]
    inf_sc=sc_dat[which(sc_dat$date<=start_date),]
    inf_sc=inf_sc[rev(order(inf_sc$date)),]
    wr_dat$pct_inf_prior_week[pp_import]=sum(inf_pp$infections_pct[1:7])
    wr_dat$pct_inf_prior_2week[pp_import]=sum(inf_pp$infections_pct[1:14])
    ##Estimate beta distribution for daily infection probability
    inf_pp$inf_logmn=NA
    inf_pp$inf_logsd=NA
    for(ii in 1:nrow(inf_pp)){
      nii=get.norm.par(p=c(0.025,.5,.975),q=log(c(inf_pp$infections_lo_pct[ii],inf_pp$infections_pct[ii],inf_pp$infections_hi_pct[ii])),show.output = F,plot=F)
      inf_pp$inf_logmn[ii]=nii[1]
      inf_pp$inf_logsd[ii]=nii[2]
    }
    print(sprintf("Week %s",pp_import))
    ##Need to create a new contact network for each week b/c recruit numbers differ
    numrecruits=wr_dat$Battalion_fill[pp_import]
    cocoon.size=ifelse(((start_date>as.Date("2020-12-01")&start_date<as.Date("2020-07-01"))|
                          start_date>as.Date("2021-12-25")),37,60) #Use 37 for 12/1/20-6/30/21 and after 12/25/21

    source("network_BT_data.R")
    wr_dat$Total_size[pp_import]=numrecruits + numdrillserg + numStaff
    R0.prob=R0_samples$R0/ mean(average.contacts) / (symp.mag + (1 - symp.mag) * asymp.adjust)
    ## specify arrival dates on campus
    arrivebase = c(sample(1:1, numrecruits, replace = T), rep(0, numdrillserg + numStaff))

    ## Importations - daily importation rate over course of training camp
    if((start_date+time_max-1)>max(sc_dat$date)){
      ##If training period goes beyond end of SC inf time series - assume constant infection rate for remainder of period
      imp_dates=which(sc_dat$date>=start_date)
      if(length(imp_dates)==0){
        ##If training period starts after time series (forecasting) - assume average of last two months
        imp_rate=mean(tail(sc_dat$pct_inf,60))
      }else{
        imp_rate=sc_dat$pct_inf[imp_dates]
      }

      imp_rate=c(imp_rate,rep(imp_rate[length(imp_rate)],times=(time_max-length(imp_rate))))
    }else{
      imp_dates=which(sc_dat$date==start_date):which(sc_dat$date==(start_date+time_max-1))
      imp_rate=sc_dat$pct_inf[imp_dates]
    }

    #importation = c(rep(0, numrecruits), rep(imp_rate, numdrillserg + numStaff))
    importation_mat = matrix(0,nrow=time_max,ncol=( numrecruits+ numdrillserg + numStaff))
    for(im in (numrecruits+1):ncol(importation_mat)){
      importation_mat[,im]=imp_rate
    }

    ## Default quarantine contacts = 0
    quarantine.contacts = 0

    for(ii_param in 1:length(testdates_list)){
        print(ii_param)
        #source('paramDefaults.R')
        time = 1:time_max

        ##Mask compliance turned off for No Masking scenarios
        if(grepl("No Mask",test_scenarios[ii_param])){
          compliance = rep(0, numrecruits + numdrillserg + numStaff)
        }else{
          compliance = rep(compliance.avg, numrecruits + numdrillserg + numStaff)
        }
        if(grepl("High Mask",test_scenarios[ii_param])){
          BOOL_compliance_date = 0
        }else{
          BOOL_compliance_date = 1
        }

        ## Change values
        BOOL_testOnArrival = testarrival_list[ii_param]

        testdates_type = test_type[[ii_param]]
        testdates= testdates_list[[ii_param]]

        ##Adjust vaccination rate for unvaccinated arrivals starting over time based on anecdotal comments from Maj Guido
        vacc_rate=ifelse(start_date<as.Date("2021-07-01"),0.25,
                         ifelse(start_date<as.Date("2021-12-26"),0.4,0.6))
        ##Set vacc_rate to 0 for no vaccination at beginning of camp scenario
        if(test_scenarios[ii_param]=="No Vaccine"){
          vacc_rate=0
        }
        if(test_scenarios[ii_param] %in% c("High Vaccine","High Vaccine Only")){
          vacc_rate=1
        }

        vacc_dat_sc_ii=vacc_dat_sc
        #If vaccinate staff and trainers set vacc rate to 100%
        if(grepl("Staff Vaccine",test_scenarios[ii_param])){
          vd1=as.Date("2021-05-22")
          vd2=as.Date("2021-05-22")+21
          vd3=as.Date("2021-12-01")
          vacc_dat_sc_ii$pct_dose_1[vacc_dat_sc_ii$Date==vd1]=1.0-sum(vacc_dat_sc_ii$pct_dose_1[vacc_dat_sc_ii$Date<vd1])
          vacc_dat_sc_ii$pct_dose_2[vacc_dat_sc_ii$Date==vd2]=1.0-sum(vacc_dat_sc_ii$pct_dose_2[vacc_dat_sc_ii$Date<vd2])
          vacc_dat_sc_ii$pct_dose_3[vacc_dat_sc_ii$Date==vd3]=1.0-sum(vacc_dat_sc_ii$pct_dose_3[vacc_dat_sc_ii$Date<vd3])
          vacc_dat_sc_ii$pct_dose_1[vacc_dat_sc_ii$Date>vd1]=0
          vacc_dat_sc_ii$pct_dose_2[vacc_dat_sc_ii$Date>vd2]=0
          vacc_dat_sc_ii$pct_dose_3[vacc_dat_sc_ii$Date>vd3]=0
        }
        
        vacc_dat_ii=vacc_dat
        #If vaccinate recruits set vacc rate to 100% 
        if(grepl("Early Vacc",test_scenarios[ii_param])&start_date>as.Date("2021-05-21")){
          vd1=start_date-15-21
          vd2=start_date-15
          vacc_dat_ii$pct_dose_1[vacc_dat_ii$Date==vd1]=1.0-sum(vacc_dat_ii$pct_dose_1[vacc_dat_ii$Date<vd1])
          vacc_dat_ii$pct_dose_2[vacc_dat_ii$Date==vd2]=1.0-sum(vacc_dat_ii$pct_dose_2[vacc_dat_ii$Date<vd2])
          vacc_dat_ii$pct_dose_1[vacc_dat_ii$Date>vd1]=0
          vacc_dat_ii$pct_dose_2[vacc_dat_ii$Date>vd2]=0
          if(start_date>=as.Date("2021-12-01")){
            vd3=as.Date("2021-12-01")
            vacc_dat_ii$pct_dose_3[vacc_dat_ii$Date==vd3]=1.0-sum(vacc_dat_ii$pct_dose_3[vacc_dat_ii$Date<vd3])
            vacc_dat_ii$pct_dose_3[vacc_dat_ii$Date>vd3]=0            
          }

        }
        
        ##Check to see if staff tested
        if(grepl("Staff Test",test_scenarios[ii_param])){
          testStaffFreq=ifelse(grepl("Daily",test_scenarios[ii_param]),1,7) #weekly testing
          testStaffType=ifelse(grepl("Antigen",test_scenarios[ii_param]),"antigen","pcr")
        }else{
          testStaffFreq=0

        }

        ## simulate
        #disp=calib_params_list$dispersion_val
        disp=R0_samples$disp
        #R0_adj_factor depends on scenario (default R0, MLE adjustment, 2x MLE adj)
        R0_adj=R0_samples$R0_adj #calib_params_list$R0_adj_factor

        num.variants=5
        # s = replicate(nreps,simOutbreak(inf_pp,inf_sc,vacc_dat,vacc_dat_sc,sampleR=T))
        sim_out_l=list()
        sim_out_l=foreach(iter=1:nreps,.packages = 'rmutil') %dopar% {
          simOutbreak(inf_pp,inf_sc,vacc_dat_ii,vacc_dat_sc_ii,sampleR=T)
        }
        ss_summary=t(sapply(sim_out_l,function(x) x[[1]],simplify = 'array'))
 

        ##Positive tests on day 1, days 2-14, and days 15+
        pos.arrival=unlist(lapply(sim_out_l, function(x) sum(x[[3]]$init.testpos.dat==1,na.rm = T)))
        pos.quarantine=unlist(lapply(sim_out_l, function(x) sum((x[[3]]$init.testpos.dat>1&x[[3]]$init.testpos.dat<=14),na.rm = T)))
        pos.training=unlist(lapply(sim_out_l, function(x) sum(x[[3]]$init.testpos.dat>14,na.rm = T)))
        num.imports=unlist(lapply(sim_out_l, function(x) sum(x[[2]]$numImported)))

        ##Infections and symptomatic infections in two time periods
        inf.quarantine=unlist(lapply(sim_out_l, function(x) sum(x[[2]]$numInfected[x[[2]]$time<=14])))
        inf.training=unlist(lapply(sim_out_l, function(x) sum(x[[2]]$numInfected[x[[2]]$time>14])))
        inf.total=unlist(lapply(sim_out_l,function(x) sum(x[[3]]$date.infected>0,na.rm=T)))
        symp.quarantine=unlist(lapply(sim_out_l, function(x) sum(x[[2]]$numSymptoms[x[[2]]$time<=14])))
        symp.training=unlist(lapply(sim_out_l, function(x) sum(x[[2]]$numSymptoms [x[[2]]$time>14])))

        ##Summary infections per rep (initial,total,positive tests on day 14)
        ind_testMat = which(init.testDay[,ncol(init.testDay)] == pp_import)
        init.testDay[ind_testMat,ii_param] = ss_summary[,2]
        inf.testDay[ind_testMat,ii_param] = ss_summary[,4]
        pos.testDay[ind_testMat,ii_param] = pos.quarantine #ss_summary[,3]
        posT.testDay[ind_testMat,ii_param] = pos.training
        imp.testDay[ind_testMat,ii_param] = num.imports
        infQ.testDay[ind_testMat,ii_param] = inf.quarantine
        infT.testDay[ind_testMat,ii_param] = inf.training
        infA.testDay[ind_testMat,ii_param] = inf.total
        infA2.testDay[ind_testMat,ii_param] = ss_summary[,4]
        sympQ.testDay[ind_testMat,ii_param] = symp.quarantine
        sympT.testDay[ind_testMat,ii_param] = symp.training
    }
}
stopCluster(cl)
init.testDay = as.data.frame(init.testDay) #number positive on arrival
inf.testDay = as.data.frame(inf.testDay) #Number of infected during simulation
pos.testDay = as.data.frame(pos.testDay) #Number positive at end of quarantine
posT.testDay = as.data.frame(posT.testDay) #Number positive after end of quarantine
imp.testDay = as.data.frame(imp.testDay) #Number imported infections from trainers and support staff
infQ.testDay = as.data.frame(infQ.testDay) #Number of infected during simulation
infT.testDay = as.data.frame(infT.testDay) #Number of infected during simulation
infA.testDay = as.data.frame(infA.testDay) #Number of infected during simulation
infA2.testDay = as.data.frame(infA2.testDay) #Number of infected during simulation
sympQ.testDay = as.data.frame(sympQ.testDay) #Number of infected during simulation
sympT.testDay = as.data.frame(sympT.testDay) #Number of infected during simulation
inf.time = as.data.frame(inf.time)

##==================================================#
## Reshape and save----------
##==================================================#
wr_dat$week_num=1:nrow(wr_dat)
colnames(inf.testDay) = c(test_scenarios, c('rep','week'))
colnames(init.testDay) = colnames(pos.testDay) = colnames(posT.testDay)=colnames(imp.testDay)= colnames(inf.testDay)
colnames(sympT.testDay)=colnames(sympQ.testDay)=colnames(infT.testDay)=colnames(infQ.testDay)= colnames(inf.testDay)
colnames(infA.testDay)=colnames(infA2.testDay)= colnames(inf.testDay)

colnames(inf.time) = c(test_scenarios, c('rep','week'))
inf.time$time = rep(1:time_max, nreps * length(importations_list))

inf.testDay2 = reshape2::melt(inf.testDay, id.vars = c('rep','week'),
                              value.name = 'infections', variable.name = 'testing')

init.testDay2 = reshape2::melt(init.testDay, id.vars = c('rep','week'),
                               value.name = 'arrival_positives', variable.name = 'testing')

pos.testDay2 = reshape2::melt(pos.testDay, id.vars = c('rep','week'),
                              value.name = 'quarantine_positives', variable.name = 'testing')

posT.testDay2 = reshape2::melt(posT.testDay, id.vars = c('rep','week'),
                               value.name = 'training_positives', variable.name = 'testing')

imp.testDay2 = reshape2::melt(imp.testDay, id.vars = c('rep','week'),
                              value.name = 'imported_infections', variable.name = 'testing')

inf.time = reshape2::melt(inf.time, id.vars = c('time','rep','week'),
                          value.name = 'infections', variable.name = 'testing')

infQ.testDay2 = reshape2::melt(infQ.testDay, id.vars = c('rep','week'),
                               value.name = 'quarantine_infections', variable.name = 'testing')

infA.testDay2 = reshape2::melt(infA.testDay, id.vars = c('rep','week'),
                               value.name = 'total_infections', variable.name = 'testing')

infA2.testDay2 = reshape2::melt(infA2.testDay, id.vars = c('rep','week'),
                                value.name = 'total_infections2', variable.name = 'testing')

infT.testDay2 = reshape2::melt(infT.testDay, id.vars = c('rep','week'),
                               value.name = 'training_infections', variable.name = 'testing')

sympQ.testDay2 = reshape2::melt(sympQ.testDay, id.vars = c('rep','week'),
                                value.name = 'symptomatic_quarantine_infections', variable.name = 'testing')

sympT.testDay2 = reshape2::melt(sympT.testDay, id.vars = c('rep','week'),
                                value.name = 'symptomatic_training_infections', variable.name = 'testing')

if(hiAsymp){
  save(ls(),file="../output/forecasted_infection_estimates_hiasymp.rda")
}else{
  save(ls(),file="../output/forecasted_infection_estimates_loasymp.rda")
}

