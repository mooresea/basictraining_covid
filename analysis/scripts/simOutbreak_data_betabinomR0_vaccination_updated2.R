## simOutbreak.R
## 
## Simulates SARS-CoV-2 transmission and interventions in a basic training setting, tracks a single 
## cohort of recruits through a 12-week basic training session
## 
## INPUTS: contact network built by network_BT.R, parameter values for infection processes, 
## parameter values for control processes(testing, isolation/quarantine/contact tracing, masks), 
## number of replicates to simulate
##
## OUTPUTS: output contains summary and time series data on number of infections, number of symptomatic 
## infections, number in isolation and quarantine, number of tests given, who infected who, and individual 
## level dates of infection/symptoms/test
##
## FEATURES:
##
## staggered arrival dates
## quarantine with coccoon from day 4 to day 14
## random contacts allowed throughout 
## test symptomatics
## test return delay
## isolate positive cases
## contact tracing and quarantine test-positive cases
## test quarantine after 3 days and release if test-negative
## individual compliance
## recent pre-arrival infection
## clinical and test-based criteria to leave isolation
## https://www.who.int/news-room/commentaries/detail/criteria-for-releasing-covid-19-patients-from-isolation 
## https://www.journalofinfection.com/article/S0163-4453(20)30119-5/fulltext#seccesectitle0014
## https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
## Contact with drill sergeants/staff

## TO DO LIST:
##
## [x] Test trainers
## [x] Implement TOST and incubation period
## [x] Implement a more thorough sensitivity
##
##
## UPDATED: Made pre-training camp infection process a day-by-day process to account for reinfections and probability of infection based on immune status and variants
## UPDATED 2: Added uncertainty in infection estimates and posterior samples for R0
## UPATED 3: ADDED adjusted probability of symptomatic infection for vaccinated individuals
##
library(rmutil)

simOutbreak = function(inf_est_s,inf_est_sc,vacc_est_s,vacc_est_sc,sampleR=F){
    
    sim_start_date=max(inf_est_s$date)
    staff_ind=c(rep(F,numrecruits),rep(T,numdrillserg+numStaff))
    
    #Make sure infection TS in chronological order (old code used reverse chronological order)
    inf_est_s=inf_est_s[order(inf_est_s$date),]
    inf_est_sc=inf_est_sc[order(inf_est_sc$date),]
    
    ## specify secondary infection distribution
    #If sampleR==T use posterior samples, otherwise use median estimate
    num.samples=1 #1000
    if(sampleR){
        R0.draw=sample(1:nrow(R0_samples),num.samples,replace=T)
        R0 = R0.prob[R0.draw] * genint #%*% t(genint) #R0.prob[R0.draw] %*% t(genint) 
        disp.val=disp[R0.draw]
        R0_adj_factor=R0_adj[R0.draw]
    }else{
        #R0.draws=1:num.samples
        R0 = R0.prob * genint #%*% t(genint) #dpois(1:21, R0.tim) / sum(dpois(1:21, R0.tim))       
        disp.val=disp
        R0_adj_factor=R0_adj_factor
    }


    ## specify delay between infection and symptom onset
    ## Replacing this line with the incub function
    ##prop.inf.symp = symp.mag * dpois(1:28, symp.tim) / sum(dpois(1:28, symp.tim))
    ## INCUBATION PERIOD
    incub = pgamma(1:21,shape=k_inc,scale=0.948) - pgamma(0:20,shape=k_inc,scale=0.948)
    prop.inf.symp = incub / sum(incub)
    
    ##Set initial variant values
    variant.mult=c(var_p$trans_mult[var_p$variant=="alpha"],
                   var_p$trans_mult[var_p$variant=="gamma"],
                   var_p$trans_mult[var_p$variant=="delta"],
                   var_p$trans_mult[var_p$variant=="omicron"],
                   var_p$trans_mult[var_p$variant=="other"])
    variant.crossprotect=matrix(1,nrow=length(variant.mult)+1,ncol=length(variant.mult)+1)
    variant.crossprotect[1,-1]=var_p$cross_protect_prob[var_p$variant=="alpha"] #Probability of cross-protection if infected with alpha
    variant.crossprotect[2,-2]=var_p$cross_protect_prob[var_p$variant=="gamma"] #Probability of cross-protection if infected with gamma
    variant.crossprotect[3,-3]=var_p$cross_protect_prob[var_p$variant=="delta"] #Probability of cross-protection if infected with delta
    variant.crossprotect[4,-4]=var_p$cross_protect_prob[var_p$variant=="omicron"]#Probability of cross-protection if infected with omicron
    variant.crossprotect[5,-5]=var_p$cross_protect_prob[var_p$variant=="other"] #Probability of cross-protection if infected with original or other
    #Protection from infection by vaccination for different variants (protect is relative to peak protection for original virus)
    variant.crossprotect[1,6]=var_p$vacc_protect_prob[var_p$variant=="alpha"] #Probability of protection via vaccination if infected with alpha
    variant.crossprotect[2,6]=var_p$vacc_protect_prob[var_p$variant=="gamma"] #Probability of protection via vaccination if infected with gamma
    variant.crossprotect[3,6]=var_p$vacc_protect_prob[var_p$variant=="delta"] #Probability of protection via vaccination if infected with delta
    variant.crossprotect[4,6]=var_p$vacc_protect_prob[var_p$variant=="omicron"] #Probability of protection via vaccination if infected with omicron
    variant.crossprotect[5,6]=var_p$vacc_protect_prob[var_p$variant=="other"]
    
    ## previously vaccinated
    ##Assume protection 14 days after dose administered
    vacc_est_s$Date_protect=vacc_est_s$Date+14
    vacc_est_sc$Date_protect=vacc_est_sc$Date+14   
    vacc_est_s=vacc_est_s[which(vacc_est_s$Date_protect<=sim_start_date),]
   
    vacc_est_s=vacc_est_s[rev(order(vacc_est_s$Date)),]
    vacc_est_sc=vacc_est_sc[rev(order(vacc_est_sc$Date)),]  
    
    initial.vaccinated.1d=sum(vacc_est_s$pct_dose_1)
    initial.vaccinated.2d=sum(vacc_est_s$pct_dose_2)
    initial.staff.vaccinated.1d=sum(vacc_est_sc$pct_dose_1[vacc_est_sc$Date_protect<sim_start_date])
    initial.staff.vaccinated.2d=sum(vacc_est_sc$pct_dose_2[vacc_est_sc$Date_protect<sim_start_date])
    #Correct for dates when 2-dose completion % higher than 1-dose
    initial.staff.vaccinated.1d=ifelse(initial.staff.vaccinated.2d>initial.staff.vaccinated.1d,initial.staff.vaccinated.2d,initial.staff.vaccinated.1d)
    
    vaccinated.m= matrix(0, numrecruits+numdrillserg+numStaff,ncol=3)
    ##vaccinated/immune among recruits
    vaccinated.m[sample(numrecruits, rbinom(1, numrecruits, initial.vaccinated.1d), replace = F),1] = 1
    if(any(colSums(vaccinated.m)>0)){
        vaccinated.m[sample(which(vaccinated.m[1:numrecruits,1]==1), rbinom(1, length(which(vaccinated.m[1:numrecruits,1]==1)), initial.vaccinated.2d/initial.vaccinated.1d), replace = F),2] = 1
    }
    ##vaccinated/immune among staff
    vaccinated.m[sample((numrecruits+1):nrow(vaccinated.m), size=rbinom(1, numdrillserg+numStaff, initial.staff.vaccinated.1d), replace = F),1] = 1
    staff2=which(vaccinated.m[,1]==1&staff_ind)
    if(length(staff2)>0){
        vaccinated.m[sample(staff2, size=rbinom(1, length(staff2), initial.staff.vaccinated.2d/initial.staff.vaccinated.1d), replace = F),2] = 1        
    }

    ##Add boosters if available
    if(any(grepl("pct_dose_3",names(vacc_est_s)))){
        initial.vaccinated.3d=sum(vacc_est_s$pct_dose_3)
        initial.staff.vaccinated.3d=sum(vacc_est_sc$pct_dose_3[vacc_est_sc$Date_protect<sim_start_date])        
        vaccinated.m[sample(numrecruits, rbinom(1, numrecruits, initial.vaccinated.3d), replace = F),3] = 1
        vaccinated.m[sample((numrecruits+1):nrow(vaccinated.m), size=rbinom(1, numdrillserg+numStaff, initial.staff.vaccinated.3d), replace = F),3] = 1
    }
    
    ## simulate vaccination before start of training camp
    date.vaccinated=matrix(NA,nrow=nrow(vaccinated.m),ncol=3)
    recruits.dose1=which(vaccinated.m[,1] > 0&!staff_ind)
    recruits.dose2=which(vaccinated.m[,2] > 0&!staff_ind)
    recruits.dose3=which(vaccinated.m[,3] > 0&!staff_ind)
    staff.dose1=which(vaccinated.m[,1] > 0&staff_ind)
    staff.dose2=which(vaccinated.m[,2] > 0&staff_ind)
    staff.dose3=which(vaccinated.m[,3] > 0&staff_ind)
    if(length(recruits.dose1)>0){
        date.vaccinated[recruits.dose1,1] = arrivebase[recruits.dose1] - 
            sample(1:nrow(vacc_est_s), length(recruits.dose1), prob=vacc_est_s$pct_dose_1,replace = T)        
    }
    if(length(recruits.dose2)>0){
        date.vaccinated[recruits.dose2,2] = arrivebase[recruits.dose2] - 
            sample(1:nrow(vacc_est_s), length(recruits.dose2), prob=vacc_est_s$pct_dose_2,replace = T)
        date.vaccinated[recruits.dose2,1] = date.vaccinated[recruits.dose2,2] - 14
    }
    if(length(recruits.dose3)>0){
        date.vaccinated[recruits.dose3,3] = arrivebase[recruits.dose3] - 
            sample(1:nrow(vacc_est_s), length(recruits.dose3), prob=vacc_est_s$pct_dose_3,replace = T)   
        date.vaccinated[recruits.dose3,2] = date.vaccinated[recruits.dose3,3] - 180
        date.vaccinated[recruits.dose3,1] = date.vaccinated[recruits.dose3,2] - 14
    }
    ##Vaccination dates for staff
    if(length(staff.dose1)>0){
        date.vaccinated[staff.dose1,1] = arrivebase[staff.dose1] - 
            sample(1:nrow(vacc_est_sc[which(vacc_est_sc$Date_protect<sim_start_date),]), length(staff.dose1), prob=vacc_est_sc$pct_dose_1[which(vacc_est_sc$Date_protect<sim_start_date)],replace = T)
    }
    if(length(staff.dose2)>0){
        date.vaccinated[staff.dose2,2] = arrivebase[staff.dose2] - 
            sample(1:nrow(vacc_est_sc[which(vacc_est_sc$Date_protect<sim_start_date),]), length(staff.dose2), prob=vacc_est_sc$pct_dose_2[which(vacc_est_sc$Date_protect<sim_start_date)],replace = T)
        date.vaccinated[staff.dose2,1] = date.vaccinated[staff.dose2,2] - 14
    }
    if(length(staff.dose3)>0){
        date.vaccinated[staff.dose3,3] = arrivebase[staff.dose3] - 
            sample(1:nrow(vacc_est_sc[which(vacc_est_sc$Date_protect<sim_start_date),]), length(staff.dose3), prob=vacc_est_sc$pct_dose_3[which(vacc_est_sc$Date_protect<sim_start_date)],replace = T)  
        date.vaccinated[staff.dose3,2] = date.vaccinated[staff.dose3,3] - 180
        date.vaccinated[staff.dose3,1] = date.vaccinated[staff.dose3,2] - 14
    }
    
    ## previously infected and (partially) immune
    infected=rep(0,numrecruits+numdrillserg+numStaff)
    variant.date.infected=matrix(NA,nrow=numrecruits+numdrillserg+numStaff,ncol=num.variants)
    date.infected = rep(NA, numrecruits+numdrillserg+numStaff)
    secondary.infected = rep(NA, numrecruits+numdrillserg+numStaff) 
    immune.source=rep(NA, numrecruits+numdrillserg+numStaff) #Vector that tracks what variant (or vaccine is most recent source of immunity)
    variant.inf=rep(NA, numrecruits+numdrillserg+numStaff) #Vector tracking which variant was most recent infection
    
    ##Before initial infections, immunity by vaccine only
    immune.source[which(rowSums(vaccinated.m)>0)]=nrow(variant.crossprotect)

    ##
    ##Infected/immune among recruits
    #immune[sample(numrecruits, rbinom(1, numrecruits, initial.immune), replace = F)] = 1
    inf.prob=exp(rnorm(nrow(inf_est_s),inf_est_s$inf_logmn,inf_est_s$inf_logsd)) #Incorporate uncertainty in infection estimate
    daily.num.inf.r=rbinom(nrow(inf_est_s),size=numrecruits,prob=inf.prob)
    #Determine who is infected for each day of time series
    for(ir in 1:length(daily.num.inf.r)){
        #determine which recruits infected on day ir
        if(daily.num.inf.r[ir]>0){
            ir.date=inf_est_s$date[ir]
            ir.sim.day=as.numeric(ir.date-sim_start_date)
            ir.var.wk=which.min(abs(ir.date-var_dat$week))
            ir.var.infs=sample(1:num.variants,size=daily.num.inf.r[ir],prob=c(var_dat$alpha_pct[ir.var.wk],
                                                                          var_dat$gamma_pct[ir.var.wk],
                                                                          var_dat$delta_pct[ir.var.wk],
                                                                          var_dat$omicron_pct[ir.var.wk],
                                                                          var_dat$other_pct[ir.var.wk]),replace=T)
            ##Set individual susceptibility levels (based on immune status) at start of each day
            #Susceptibility based on infection
            days_since_pk_titer=ir.sim.day-date.infected-imm_decay_start
            days_since_pk_titer=ifelse(days_since_pk_titer>0,days_since_pk_titer,0)
            suscept.day.inf=ifelse(is.na(date.infected),1,
                                   1-exp(-imm_decay_rate*days_since_pk_titer))
            #Susceptibility based on vaccination
            ir.date.vaccinated=apply(date.vaccinated,2,function(x) ifelse(x<=ir.sim.day,x,NA)) #Only consider vaccines after IR date
            days_since_vacc=ir.sim.day-ir.date.vaccinated #-vacc_decay_start
            vacc_susc_d1=1-(vacc_1_dose_eff/vacc_2_dose_eff)*ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,1])))
            vacc_susc_d2=1-ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,2])))
            vacc_susc_d3=1-(vacc_3_dose_eff/vacc_2_dose_eff)*ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,3])))
            vacc_susc_d1[is.na(vacc_susc_d1)]=1
            vacc_susc_d2[is.na(vacc_susc_d2)]=1
            vacc_susc_d3[is.na(vacc_susc_d3)]=1
            suscept.day.v=pmin(vacc_susc_d1,vacc_susc_d2,vacc_susc_d3)
            suscept.ir=pmin(suscept.day.inf,suscept.day.v)
            
            ##Determine who infected
            for(jr in 1:daily.num.inf.r[ir]){
                jr.var=ir.var.infs[jr]
                suscept.jr=1-(1-suscept.ir[!staff_ind])*ifelse(is.na(immune.source[!staff_ind]),0,variant.crossprotect[jr.var,immune.source[!staff_ind]])
                jr.inf=sample(1:numrecruits,size=1,prob=suscept.jr)
                infected[jr.inf]=1
                if(!is.na(date.infected[jr.inf])){
                    secondary.infected[jr.inf]=date.infected[jr.inf]
                }
                date.infected[jr.inf]=ir.sim.day
                immune.source[jr.inf]=jr.var
                variant.inf[jr.inf]=jr.var
                variant.date.infected[jr.inf,jr.var]=ir.sim.day
            }
        }
    }
    
    ##Infected/immune among staff
    #immune[sample((numrecruits+1):length(immune), size=rbinom(1, numdrillserg+numStaff, initial.staff.immune), replace = F)] = 1
    daily.num.inf.s=rbinom(nrow(inf_est_sc),size=numdrillserg+numStaff,prob=inf_est_sc$pct_inf)
    for(iis in 1:length(daily.num.inf.s)){
        #determine which recruits infected on day ir
        if(daily.num.inf.s[iis]>0){
            is.date=inf_est_s$date[iis]
            is.sim.day=as.numeric(is.date-sim_start_date)
            is.var.wk=which.min(abs(is.date-var_dat$week))
            is.var.infs=sample(1:num.variants,size=daily.num.inf.s[iis],prob=c(var_dat$alpha_pct[is.var.wk],
                                                                              var_dat$gamma_pct[is.var.wk],
                                                                              var_dat$delta_pct[is.var.wk],
                                                                              var_dat$omicron_pct[is.var.wk],
                                                                              var_dat$other_pct[is.var.wk]),replace=T)
            ##Set individual susceptibility levels (based on immune status) at start of each day
            #Susceptibility based on infection
            days_since_pk_titer=is.sim.day-date.infected-imm_decay_start
            days_since_pk_titer=ifelse(days_since_pk_titer>0,days_since_pk_titer,0)
            suscept.day.inf=ifelse(is.na(date.infected),1,
                                   1-exp(-imm_decay_rate*days_since_pk_titer))
            #Susceptibility based on vaccination
            is.date.vaccinated=apply(date.vaccinated,2,function(x) ifelse(x<=is.sim.day,x,NA)) #Only consider vaccines after IR date
            days_since_vacc=is.sim.day-is.date.vaccinated #-vacc_decay_start
            vacc_susc_d1=1-(vacc_1_dose_eff/vacc_2_dose_eff)*ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,1])))
            vacc_susc_d2=1-ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,1])))
            vacc_susc_d3=1-(vacc_3_dose_eff/vacc_2_dose_eff)*ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,1])))
            vacc_susc_d1[is.na(vacc_susc_d1)]=1
            vacc_susc_d2[is.na(vacc_susc_d2)]=1
            vacc_susc_d3[is.na(vacc_susc_d3)]=1
            suscept.day.v=pmin(vacc_susc_d1,vacc_susc_d2,vacc_susc_d3)
            suscept.is=pmin(suscept.day.inf,suscept.day.v)
            
            ##Determine who infected
            for(js in 1:daily.num.inf.s[iis]){
                js.var=is.var.infs[js]
                suscept.js=1-(1-suscept.is[staff_ind])*ifelse(is.na(immune.source[staff_ind]),0,variant.crossprotect[js.var,immune.source[staff_ind]])
                js.inf=sample((numrecruits+1):(numrecruits+numdrillserg+numStaff),size=1,prob=suscept.js)
                infected[js.inf]=1
                if(!is.na(date.infected[js.inf])){
                    secondary.infected[js.inf]=date.infected[js.inf]
                }
                date.infected[js.inf]=is.sim.day
                immune.source[js.inf]=js.var
                variant.inf[js.inf]=js.var
                variant.date.infected[js.inf,js.var]=is.sim.day
            }
        }
    }  
    init.infect.r=sum(infected[!staff_ind],na.rm=T)
    init.infect.s=sum(infected[staff_ind],na.rm=T)
    
    # Update immune status to account for individuals who have been vaccinated since they were infected
    recent.vacc.dates=suppressWarnings(apply(date.vaccinated,1,function(x) min(x,na.rm=T))) #Suppress warnings about infinite values when all NA
    recent.vacc.dates=ifelse(is.infinite(recent.vacc.dates),-9999,recent.vacc.dates)
    inf.dates=ifelse(is.na(date.infected),-8888,date.infected)
    immune.source[which(recent.vacc.dates>inf.dates)]=nrow(variant.crossprotect)

    ## simulate date of symptom onset and duration
    incub.period.inf = rep(NA, numrecruits+numdrillserg+numStaff)
    date.symptoms = rep(NA, numrecruits+numdrillserg+numStaff)
    symptom.duration = rep(NA, numrecruits+numdrillserg+numStaff)
    for(ii in which(infected > 0)){
        ii_symp_mag=ifelse(vaccinated.m[ii,1]==1,vacc.symp.adj,symp.mag)
        if(rbinom(1, 1, ii_symp_mag) == 1){
            date.symptoms[ii] =
                date.infected[ii] +
                sample(1:length(prop.inf.symp), 1, prob = prop.inf.symp)
            symptom.duration[ii] = rpois(1,symp.mean)
            incub.period.inf[ii] = date.symptoms[ii]
        }else{
            incub.period.inf[ii] = date.infected[ii] +
                sample(1:length(prop.inf.symp),1,replace = T,prob = prop.inf.symp)
        }
    }
    
    ## individual status storage
    tested.date = rep(NA, numrecruits+numdrillserg+numStaff)
    isolate.date = rep(NA, numrecruits+numdrillserg+numStaff)
    quarantine.date = rep(NA, numrecruits+numdrillserg+numStaff)    
    testneg.date = rep(NA, numrecruits+numdrillserg+numStaff)
    init.testpos.date=rep(NA, numrecruits+numdrillserg+numStaff)
    suscept.day=rep(NA, numrecruits+numdrillserg+numStaff)
    init.pos = 0
    retest.pos = 0

    ## aggregate numbers over time
    numIsolated = rep(0, length(time))
    numQuarantine = rep(0, length(time))
    numTested = rep(0, length(time))
    numInfected = rep(0, length(time))
    numSymptoms = rep(0, length(time))
    numImported = rep(0, length(time))
    numTestPositive_daily = rep(0, length(time))
    numTestPositive_sympTesting = rep(0, length(time))
    
    ## keep track of who infected whom
    edges = cbind(rep(0, sum(infected)), which(infected > 0))
    
    ##
    ## loop over each day of basic training
    ##
    inf_cntr=c()
    for(tt in time){
        compliance.mult = 16/24
        if(BOOL_compliance_time == 1){
            compliance.mult = ((compliance.final.avg - compliance.avg)/(length(time) - 1) * (tt-1) + compliance.avg) / compliance.avg
        }
        if(BOOL_compliance_date == 1){
            #Compliance based on date - from DELPHI surveys - Assume compliance was constant before first date in Feb 2021
            compliance.mult = mask_dat$value[which.min(abs(mask_dat$date-(sim_start_date+tt-1)))]*(16/24) #Assume compliance reduced by time in barracks, etc
        }
        ##Update vaccination status of staff - go by date of protection (so reflects vaccination from several weeks ago)
        tt_1st_dose_pct=vacc_est_sc$pct_dose_1[vacc_est_sc$Date_protect==(sim_start_date+tt)]
        tt_2nd_dose_pct=vacc_est_sc$pct_dose_2[vacc_est_sc$Date_protect==(sim_start_date+tt)]
        tt_1st_dose_pct=ifelse(length(tt_1st_dose_pct)==0,0,tt_1st_dose_pct)
        tt_2nd_dose_pct=ifelse(length(tt_2nd_dose_pct)==0,0,tt_2nd_dose_pct)
        if(any(grep("pct_dose_3",names(vacc_est_sc)))){
            tt_3rd_dose_pct=vacc_est_sc$pct_dose_3[vacc_est_sc$Date_protect==(sim_start_date+tt)]            
        }else{
            tt_3rd_dose_pct=0
        }
        tt_3rd_dose_pct=ifelse(length(tt_3rd_dose_pct)==0,0,tt_3rd_dose_pct)
        
        tt_num_dose_1=rbinom(1,size=numStaff+numdrillserg,tt_1st_dose_pct)
        tt_num_dose_2=rbinom(1,size=numStaff+numdrillserg,tt_2nd_dose_pct)
        tt_num_dose_3=rbinom(1,size=numStaff+numdrillserg,tt_3rd_dose_pct)
        if(tt_num_dose_1>length(which(vaccinated.m[,1]==0&staff_ind))){
            tt_num_dose_1=length(which(vaccinated.m[,1]==0&staff_ind))
        }
        if(tt_num_dose_2>length(which(vaccinated.m[,2]==0&staff_ind&vaccinated.m[,1]==1))){
            tt_num_dose_2=length(which(vaccinated.m[,2]==0&staff_ind&vaccinated.m[,1]==1))
        }
        if(tt_num_dose_3>length(which(vaccinated.m[,3]==0&staff_ind&vaccinated.m[,2]==1))){
            tt_num_dose_3=length(which(vaccinated.m[,3]==0&staff_ind&vaccinated.m[,2]==1))
        }
        tt_staff_dose_1=sample(which(vaccinated.m[,1]==0&staff_ind),size=tt_num_dose_1,replace=F)
        tt_staff_dose_2=sample(which(vaccinated.m[,2]==0&staff_ind&vaccinated.m[,1]==1),size=tt_num_dose_2,replace=F)
        tt_staff_dose_3=sample(which(vaccinated.m[,3]==0&staff_ind&vaccinated.m[,2]==1),size=tt_num_dose_3,replace=F)
        vaccinated.m[tt_staff_dose_1,1]=1
        date.vaccinated[tt_staff_dose_1,1]=tt
        vaccinated.m[tt_staff_dose_2,2]=1
        date.vaccinated[tt_staff_dose_2,2]=tt
        vaccinated.m[tt_staff_dose_3,3]=1
        date.vaccinated[tt_staff_dose_3,3]=tt
        immune.source[c(tt_staff_dose_1,tt_staff_dose_2,tt_staff_dose_3)]=ifelse(immune.source[c(tt_staff_dose_1,tt_staff_dose_2,tt_staff_dose_3)]==4,immune.source[c(tt_staff_dose_1,tt_staff_dose_2,tt_staff_dose_3)],nrow(variant.crossprotect))
        
        ##Vaccinate recruits near start of training - build in 14 day lag for first dose to provide protection
        if(include_vacc&sim_start_date>=vacc_start_date&tt==(vacc_start+14)){
            num_doses=round(length(which(vaccinated.m[,2]==0&!staff_ind))*vacc_rate,0) #Vaccinate only those who have not already received a 2nd dose
            
            tt_recruit_vacc=sample(which(vaccinated.m[,2]==0&!staff_ind),size=num_doses,replace=F)
            tt_recruit_dose1=which(vaccinated.m[tt_recruit_vacc,1]==0)
            tt_recruit_dose2=which(vaccinated.m[tt_recruit_vacc,1]>0)
            #Booster doses
            if(include_booster){
                num_doses_3=round(length(which(vaccinated.m[,3]==0&!staff_ind&vaccinated.m[,2]==1))*booster_vacc_rate,0) #Vaccinate only those who have not already received a 3rd dose
                tt_recruit_dose3=sample(which(vaccinated.m[,3]==0&!staff_ind&vaccinated.m[,2]==1),size=num_doses_3,replace=F)      
                vaccinated.m[tt_recruit_dose3,3]=1
                date.vaccinated[tt_recruit_dose3,3]=tt
                #If recruit already infected with omicron keep that as source of immunity, otherwise make vaccine source of immunity
                immune.source[c(tt_recruit_dose3)]=ifelse(immune.source[c(tt_recruit_dose3)]==4,immune.source[c(tt_recruit_dose3)],nrow(variant.crossprotect))
                
            }

            
            vaccinated.m[tt_recruit_dose1,1]=1
            date.vaccinated[tt_recruit_dose1,1]=tt
            vaccinated.m[tt_recruit_dose2,2]=1
            date.vaccinated[tt_recruit_dose2,2]=tt
            #If recruit already infected with omicron keep that as source of immunity, otherwise make vaccine source of immunity
            immune.source[c(tt_recruit_dose1,tt_recruit_dose2)]=ifelse(immune.source[c(tt_recruit_dose1,tt_recruit_dose2)]==4,immune.source[c(tt_recruit_dose1,tt_recruit_dose2)],nrow(variant.crossprotect))
        }
        ##Time for protection to 2nd dose (for those who received 1st dose after arrival)
        if(include_vacc&sim_start_date>=vacc_start_date&tt==(vacc_start+14+21)){
            vaccinated.m[tt_recruit_dose1,2]=1
            date.vaccinated[tt_recruit_dose1,2]=tt
            #If recruit already infected with omicron keep that as source of immunity, otherwise make vaccine source of immunity
            immune.source[c(tt_recruit_dose1)]=ifelse(immune.source[c(tt_recruit_dose1)]==4,immune.source[c(tt_recruit_dose1)],nrow(variant.crossprotect))
        }
        ##
        ##Set individual susceptibility levels (based on immune status) at start of each day
        ##
        #Susceptibility based on infection
        days_since_pk_titer=tt-date.infected-imm_decay_start
        days_since_pk_titer=ifelse(days_since_pk_titer>0,days_since_pk_titer,0)
        suscept.day.inf=ifelse(is.na(date.infected),1,
                             1-exp(-imm_decay_rate*days_since_pk_titer))
        #Susceptibility based on vaccination
        days_since_vacc=tt-date.vaccinated #-vacc_decay_start
        vacc_susc_d1=1-(vacc_1_dose_eff/vacc_2_dose_eff)*ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,1])))
        vacc_susc_d2=1-ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,2])))
        vacc_susc_d3=1-(vacc_3_dose_eff/vacc_2_dose_eff)*ve_peak/(1+exp(-(ve_p1+ve_dr*days_since_vacc[,3])))
        vacc_susc_d1[is.na(vacc_susc_d1)]=1
        vacc_susc_d2[is.na(vacc_susc_d2)]=1
        vacc_susc_d3[is.na(vacc_susc_d3)]=1
        suscept.day.v=pmin(vacc_susc_d1,vacc_susc_d2,vacc_susc_d3)
        
        suscept.day=pmin(suscept.day.inf,suscept.day.v)
        ##Get variant percentages for particular day
        tt_var_week=which.min(abs(sim_start_date+tt-var_dat$week))
        
        ##Reset individual infection risk by variant for each day
        
        ## test on arrival
        ## What to do with those who test positive for antibodies? ---> they don't get isolated on day 0
        if(tt %in% unique(arrivebase) & BOOL_testOnArrival == 1){
            needtesting = which(arrivebase == tt)
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))
            if(length(trueneg) > 0){ # false pos
                testpos[which(needtesting %in% trueneg)] =
                    rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
            }
            ## Upon arrival, always test with PCR
            if(length(truepos) > 0){ # true pos
                testpos[which(needtesting %in% truepos)] =
                    rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] + 1))
            }
            ## Only trace and isolate those seronegative 
            tested.date[needtesting] = tt
            needtracing = needtesting[which(testpos == 1)] # all pos
            init.testpos.date[needtracing] = tt #Individuals testing positive at initial test
            #ind_testpos = which(testpos == 1)
            #needtracing = needtesting[ind_testpos[which(!(ind_testpos %in% which(pos.immune == 1)))]] # all pos
            isolate.date[needtracing] = tt + testReturn
            
            ## TROUBLESHOOTING CODE
            init.infs=needtracing #which(!is.na(isolate.date))
            
            testneg.date[setdiff(needtesting,needtracing)] = tt
            numTested[tt] = numTested[tt] + length(needtesting)
            init.pos = init.pos + length(needtracing)
        }

        ## Only if want to test everyday, for bookkeeping
        if(BOOL_testDaily == TRUE){
            needtesting = 1:numrecruits
            truepos = needtesting[!is.na(date.infected[needtesting])]
            trueneg = numrecruits - length(truepos)
            testpos = 0
            if(trueneg > 0){
                testpos = sum(rbinom(trueneg, 1, 1 - pcr.spec.commercial))
            }
            if(length(truepos) > 0){
                testpos = testpos + sum(
                    rbinom(length(truepos), 1,
                           pcr.sens.commercial(tt - date.infected[truepos] + 1)
                           ))
            }
            numTestPositive_daily[tt] = testpos
        }
        
        ## test on specified days
        if(tt %in% testdates){
            needtesting = 1:numrecruits
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))

            ind_testdates = which(testdates == tt)
            tmp_test_type = 'pcr'
            if(ind_testdates <= length(testdates_type)){
                tmp_test_type = testdates_type[ind_testdates]
            }
            ## can be antigen or PCR, if not antigen, assume pcr
            if(tmp_test_type == "antigen"){
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.screen)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.screen(tt - date.infected[truepos] +
                                                                   1))
                }
            }else{
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                   1))
                }
            }
            tested.date[needtesting] = tt
            needtracing = needtesting[which(testpos == 1)]
            isolate.date[needtracing] = tt + testReturn
            testneg.date[setdiff(needtesting,needtracing)] = tt
            init.testpos.date[needtracing] = ifelse(is.na(init.testpos.date[needtracing]),tt,init.testpos.date[needtracing]) #Individuals testing positive who weren't already positive
            numTested[tt] = numTested[tt] + length(needtesting)
            retest.pos = retest.pos + length(needtracing)
            
            for(ii in needtracing){
                if(ii <= numrecruits){ # recruits                    
                    if(tt <= 14){
                        contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                    } else {
                        contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                    } 
                } else if(ii > numrecruits + numdrillserg){ # staff
                    contacts.all = contacts.staff.recruit[[ii - numrecruits - numdrillserg]]
                }else { # drill sergeants
                    contacts.all = contacts.drillserg.recruit[[ii - numrecruits]]
                }

                contacts.all = contacts.all[contacts.all != ii]
                contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))

                ## CONTACT TRACING!! quarantine.contacts and set date to isolate if positive
                contacts.needtest = sample(contacts.all,
                       min(
                           rpois(1, quarantine.contacts),
                           length(contacts.all)
                       ), replace = F)
                
                trueneg = contacts.needtest[is.na(date.infected[contacts.needtest])]
                truepos = setdiff(contacts.needtest, trueneg)
                testpos = rep(NA, length(contacts.needtest))

                ## can be antigen or PCR, if not antigen, assume pcr
                if(tmp_test_type == "antigen"){
                    if(length(trueneg) > 0){
                        testpos[which(contacts.needtest %in% trueneg)] =
                            rbinom(length(trueneg), 1, 1 - pcr.spec.screen)
                    }
                    if(length(truepos) > 0){
                        testpos[which(contacts.needtest %in% truepos)] =
                            rbinom(length(truepos), 1, pcr.sens.screen(tt - date.infected[truepos] +
                                                                       1))
                    }
                }else{
                    if(length(trueneg) > 0){
                        testpos[which(contacts.needtest %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                    }
                    if(length(truepos) > 0){
                        testpos[which(contacts.needtest %in% truepos)] =
                            rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                           1))
                    }
                }
                tested.date[contacts.needtest] = tt                
                quarantine.date[contacts.needtest[testpos == 1]] = tt + testReturn # Can also be isolate.date
                testneg.date[contacts.needtest[testpos != 1]] = tt
                numTested[tt] = numTested[tt] + length(contacts.needtest)
                init.testpos.date[contacts.needtest[testpos == 1]] = ifelse(is.na(init.testpos.date[contacts.needtest[testpos == 1]]),tt,init.testpos.date[contacts.needtest[testpos == 1]])
                retest.pos = retest.pos + length(which(testpos == 1))            
            }
        }


        ## test staff on specified frequencies
        if((testStaffFreq > 0) && ((tt - 1) %% testStaffFreq == 0)){
            needtesting = (numrecruits + 1):length(infected) # Is this for staff or drill sergants too?
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))
           
            ## can be antigen or PCR, if not antigen, assume pcr
            if(testStaffType == "antigen"){
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.screen)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.screen(tt - date.infected[truepos] +
                                                                   1))
                }
            }else{
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                       1))
                }
            }
            tested.date[needtesting] = tt
            
            needtracing = needtesting[which(testpos == 1)]
            isolate.date[needtracing] = tt + testReturn
            testneg.date[setdiff(needtesting,needtracing)] = tt
            init.testpos.date[needtracing] = ifelse(is.na(init.testpos.date[needtracing]),tt,init.testpos.date[needtracing])
            numTested[tt] = numTested[tt] + length(needtesting)
            retest.pos = retest.pos + length(needtracing)                        
        }
        
        ##Need to adjust for differential susceptibility based on variant - use pct weighted cross-protection
        suscept.day.alpha=1-(1-suscept.day)*ifelse(is.na(immune.source),0,variant.crossprotect[1,immune.source])
        suscept.day.gamma=1-(1-suscept.day)*ifelse(is.na(immune.source),0,variant.crossprotect[2,immune.source])
        suscept.day.delta=1-(1-suscept.day)*ifelse(is.na(immune.source),0,variant.crossprotect[3,immune.source])
        suscept.day.omicron=1-(1-suscept.day)*ifelse(is.na(immune.source),0,variant.crossprotect[4,immune.source])
        suscept.day.other=1-(1-suscept.day)*ifelse(is.na(immune.source),0,variant.crossprotect[5,immune.source])
        suscept.day.var=suscept.day.alpha*var_dat$alpha_pct[tt_var_week]+suscept.day.gamma*var_dat$gamma_pct[tt_var_week]+
            suscept.day.delta*var_dat$delta_pct[tt_var_week]+suscept.day.omicron*var_dat$omicron_pct[tt_var_week]+suscept.day.other*var_dat$other_pct[tt_var_week]

        ## importation from staff - infection rate in Richland county* individual susceptibility   
        #Number to infect
        infected.offcampus = which(rbinom(numrecruits+numdrillserg+numStaff, 1,suscept.day.var*(1-compliance*compliance.mult)*as.numeric(importation_mat[tt,])) == 1)
        numImported[tt] = length(infected.offcampus)
        if(length(infected.offcampus) > 0){
            #Simplify code to allow and record reinfections 
            edges = rbind(
                    edges,
                    cbind(0,infected.offcampus))
            imp.2nd = infected.offcampus[which(infected[infected.offcampus] == 1)]
            secondary.infected[imp.2nd] = date.infected[imp.2nd]
            infected[infected.offcampus] = 1
            suscept.day[infected.offcampus] = 0 #Need to make sure not infected again by one of their contacts
            ##Set variant
            variant.inf[infected.offcampus]=sample(1:length(variant.mult),size=length(infected.offcampus),
                                                    prob=c(var_dat$alpha_pct[tt_var_week],
                                                           var_dat$gamma_pct[tt_var_week],
                                                           var_dat$delta_pct[tt_var_week],
                                                           var_dat$omicron_pct[tt_var_week],
                                                           var_dat$other_pct[tt_var_week]),replace = T)
            immune.source[infected.offcampus]= variant.inf[infected.offcampus]
            date.infected[infected.offcampus] = tt
            for(ii in infected.offcampus){
                ii_symp_mag=ifelse(vaccinated.m[ii,1]==1,vacc.symp.adj,symp.mag)
                date.symptoms[ii] =
                 ifelse(
                    rbinom(1, 1, ii_symp_mag) == 1,
                    date.infected[ii] +
                    sample(
                        1:length(prop.inf.symp),
                        1,
                        prob = prop.inf.symp,
                        replace = T),
                    NA)
            }
            incub.period.inf[infected.offcampus] = date.symptoms[infected.offcampus]
            
            
            incub.period.inf[infected.offcampus[is.na(date.symptoms[infected.offcampus])]] = date.infected[infected.offcampus[is.na(date.symptoms[infected.offcampus])]] +
                sample(1:length(prop.inf.symp),
                       length(infected.offcampus[is.na(date.symptoms[infected.offcampus])]),
                       replace = T, prob = prop.inf.symp)
        }
        
        ## loop through those who are capable of infecting others
        max_infectious_period = 21
        infectious = which(infected > 0 &
                           date.infected < tt &
                           date.infected >(tt - max_infectious_period) &
                           arrivebase <= tt)
    
        for(ii in infectious){
            # infect.today =
            #     R0[tt - date.infected[ii] + 1] *
            #     ifelse(is.na(date.symptoms[ii]), asymp.adjust, 1) * (1 - compliance[ii])
            #ii.R0.draw=sample(1:num.samples,size=1)
            if(tt<=14){
                infect.today =
                    variant.mult[variant.inf[ii]]*R0[((tt - date.infected[ii]) + 1)] *
                    ifelse(is.na(date.symptoms[ii]), asymp.adjust, 1) * (1 - compliance[ii]*compliance.mult)                
            }else{
                infect.today =
                    R0_adj_factor*variant.mult[variant.inf[ii]]*R0[((tt - date.infected[ii]) + 1)] *
                    ifelse(is.na(date.symptoms[ii]), asymp.adjust, 1) * (1 - compliance[ii]*compliance.mult)
            }

            ## look up this person's contacts today
            if(!is.na(isolate.date[ii]) | !is.na(quarantine.date[ii])){
                if(ii <= numrecruits){
                    which.in.isolation = unique(c(which(!is.na(isolate.date)), which(!is.na(quarantine.date))))
                    contacts.all = which.in.isolation
                }else{
                    contacts.all = numeric()
                }
            } else {
                if(ii <= numrecruits){ # recruits
                    if(tt <= 14){
                        contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                    } else {
                        contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                    }
                } else if(ii > numrecruits + numdrillserg){ # staff
                    contacts.all = unique(c(contacts.staff.recruit[[ii - numrecruits - numdrillserg]]))
                }else { # drill sergeants
                    contacts.all = unique(c(contacts.drillserg.recruit[[ii - numrecruits]]))
                }
            }
            contacts.all = contacts.all[contacts.all != ii]

            ##Susceptibility of contacts adjusted for variant of ii
            variant.cp.ii=ifelse(is.na(variant.crossprotect[variant.inf[ii],immune.source[contacts.all]]),0,variant.crossprotect[variant.inf[ii],immune.source[contacts.all]])
            susc.contacts.ii=1-(1-suscept.day[contacts.all])*variant.cp.ii
            ##Set susceptibility for individuals infected in last 30 days to 0
            susc.contacts.ii=ifelse((tt-date.infected[contacts.all])>imm_total_init|is.na(date.infected[contacts.all]),susc.contacts.ii,0)
                                        
           if(length(contacts.all)==0|infect.today*mean(1 - compliance[contacts.all]*compliance.mult)*mean(susc.contacts.ii)==0){
               infect.num=0
           }else{
               contacts.susc=contacts.all[which(susc.contacts.ii>0)]
               infect.num= rbetabinom(1,size=length(contacts.susc),
                                      m=infect.today*mean(1 - compliance[contacts.susc]*compliance.mult)*mean(susc.contacts.ii[which(susc.contacts.ii>0)]),s=disp.val)             
                #infect.num= rbinom(1,length(contacts.all),
                #         infect.today *(1 - compliance[contacts.all]*compliance.mult))
           }
            inf_cntr=c(inf_cntr,infect.num)
            if(infect.num > 0){
                infect.who=sample(1:length(contacts.susc),infect.num,prob=(1 - compliance[contacts.susc]*compliance.mult)*susc.contacts.ii[which(susc.contacts.ii>0)],replace=F)
                infect.who = contacts.susc[infect.who]
                
                ## update their status if they're not already infected
                infect.2nd = infect.who[which(infected[infect.who] == 1)]
                secondary.infected[infect.2nd] = date.infected[infect.2nd]
                infected[infect.who] = 1
                suscept.day[infect.who] = 0 #Need to make sure not reinfected immediately by one of their contacts
                ##Set variant
                variant.inf[infect.who]=variant.inf[ii]
                immune.source[infect.who]=variant.inf[ii]
                
                date.infected[infect.who] = tt
                for(ii in infect.who){
                    ii_symp_mag=ifelse(vaccinated.m[ii,1]==1,vacc.symp.adj,symp.mag)
                    date.symptoms[ii] =
                        ifelse(
                            rbinom(1, 1, ii_symp_mag) == 1,
                            date.infected[ii] +
                                sample(
                                    1:length(prop.inf.symp),
                                    1,
                                    prob = prop.inf.symp,
                                    replace = T),
                            NA)
                }
                incub.period.inf[infect.who] = date.symptoms[infect.who]                    
                incub.period.inf[infect.who[is.na(date.symptoms[infect.who])]] = date.infected[infect.who[is.na(date.symptoms[infect.who])]] +
                        sample(1:length(prop.inf.symp),
                               length(infect.who[is.na(date.symptoms[infect.who])]),
                               replace = T,
                               prob = prop.inf.symp)
                    symptom.duration[infect.who] = rpois(1,symp.mean)
                edges = rbind(edges, cbind(ii, infect.who))
            }
        } # end infectious loop
        
        ## test symptomatics, and perform isolation and quarantine accordingly
        needtesting = which(date.symptoms == tt)
        needtesting = needtesting[which(arrivebase[needtesting] < tt)]
        needtesting = needtesting[which(is.na(isolate.date[needtesting]))]
        needtesting = needtesting[which(is.na(quarantine.date[needtesting]))]
        needtesting = needtesting[which(tested.date[needtesting] != tt|is.na(tested.date[needtesting]))]
        numTested[tt] = numTested[tt] + length(needtesting)
        if(length(needtesting) > 0){
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))
            if(length(trueneg) > 0){
                testpos[which(needtesting %in% trueneg)] =
                    rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
            }
            if(length(truepos) > 0){
                testpos[which(needtesting %in% truepos)] =
                    rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                   1))
            }
            tested.date[needtesting] = tt
            needtracing = needtesting[which(testpos == 1)]
            ##Get number testing positive through symptomatic testing
            numTestPositive_sympTesting[tt]=sum(testpos == 1)
            isolate.date[needtracing] = tt # isolate on test date
            init.testpos.date[needtracing] = ifelse(is.na(init.testpos.date[needtracing]),tt,init.testpos.date[needtracing])
            testneg.date[setdiff(needtesting,needtracing)] = tt
            for(ii in needtracing){
                if(ii <= numrecruits){ # recruits
                    if(tt <= 14){
                        contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                    } else {
                        contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                    } 
                } else if(ii > numrecruits + numdrillserg){ # staff
                    contacts.all = contacts.staff.recruit[[ii - numrecruits - numdrillserg]]
                }else { # drill sergeants
                    contacts.all = contacts.drillserg.recruit[[ii - numrecruits]]
                }
                contacts.all = contacts.all[contacts.all != ii]
                contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))
                quarantine.date[sample(contacts.all, min(
                                                         rpois(1, quarantine.contacts),
                                                         length(contacts.all)
                                                     ), replace = F)] =
                    tt + 1
            }
        }
        
        ## test those who were recently quarantined and release if negative        
        if(tt > testDelayQuarantine){
            needtesting = which(quarantine.date ==(tt - testDelayQuarantine))
            needtesting = needtesting[which(tested.date[needtesting] != tt)]
            numTested[tt] = numTested[tt] + length(needtesting)
            if(length(needtesting) > 0){
                trueneg = needtesting[is.na(date.infected[needtesting] + 1)]
                truepos = setdiff(needtesting, trueneg)
                testpos = rep(NA, length(needtesting))
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos),
                               1,
                               pcr.sens.commercial(tt - date.infected[truepos] + 1))
                }
                tested.date[needtesting] = tt
                release = which(testpos == 0)
                if(length(release) > 0){
                    quarantine.date[needtesting[release]] = NA
                }
                needtracing = needtesting[which(testpos == 1)]
                ##Add positives in quarantine
                numTestPositive_sympTesting[tt]=numTestPositive_sympTesting[tt]+sum(testpos==1)
                init.testpos.date[needtracing] = ifelse(is.na(init.testpos.date[needtracing]),tt,init.testpos.date[needtracing])
                for(ii in needtracing){
                    if(ii <= numrecruits){ # recruits
                        if( tt <= 14){
                            contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                        } else {
                            contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                        } 
                    } else if(ii > numrecruits + numdrillserg){ # staff
                        contacts.all = contacts.staff.recruit[[ii - numrecruits - numdrillserg]]
                    }else { # drill sergeants
                        contacts.all = contacts.drillserg.recruit[[ii - numrecruits]]
                    }
                    contacts.all = contacts.all[contacts.all != ii]
                    contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))
                    quarantine.date[sample(contacts.all, min(
                                                             rpois(1, quarantine.contacts),
                                                             length(contacts.all)
                                                         ), replace = F)] =
                        tt + 1
                }
            }
        }
        
        ## record numbers for the day
        numIsolated[tt] = sum(!is.na(isolate.date))
        numQuarantine[tt] = sum(!is.na(quarantine.date))
        numInfected[tt] = sum(date.infected == tt, na.rm = T)
        numSymptoms[tt] = sum(date.symptoms == tt, na.rm = T)
        
        ## release from isolation
        if(BOOL_clinicalrelease){ # release based on clinical criteria 
            release = c(which(max(date.symptoms+isolate.length, # 10 days from symptom onset 
                                  date.symptoms+symptom.duration+isolate.nosymp) == tt), # delay from last day of fever
                        which(tested.date[which(is.na(date.symptoms))]+isolate.length == tt), # 10 days from test date for asymptomatics
                        which(testneg.date == tt - testReturn)) # negative tests returned 
        } else { # release based on testing
            needtesting = which(isolate.date <=(tt - isolate.length))
            numTested[tt] = numTested[tt] + length(needtesting)
            if(length(needtesting) > 0){
                trueneg = needtesting[is.na(tt - date.infected[needtesting] + 1)]
                truepos = setdiff(needtesting, trueneg)
                testpos = rep(NA, length(needtesting))
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos),
                               1,
                               pcr.sens.commercial(tt - date.infected[truepos] + 1))
                }
                tested.date[needtesting] = tt
                release = which(testpos == 0)
            } else {
                release = numeric()
            }
        }
        
        if(length(release) > 0){
            isolate.date[release] = NA
        }
        
        ## release from quarantine if time up
        release = which(quarantine.date ==(tt - quarantine.max))
        if(length(release) > 0){
            quarantine.date[release] = NA
        }
        
    } # end time loop
 
    ## outputs
    summary_out = c(init.infect.r, 
                    init.pos,
                    retest.pos,
                    sum(numInfected), 
                    sum(numSymptoms),
                    max(numIsolated + numQuarantine),
                    sum(numTested),
                    length(which(!is.na(secondary.infected))))
    timeseries_out = data.frame(cbind(time,numInfected,numSymptoms,numIsolated,numQuarantine,numTested, numTestPositive_sympTesting, numImported, numTestPositive_daily))
    indivudal_out = data.frame(cbind(date.symptoms,isolate.date,quarantine.date,init.testpos.date,date.infected,secondary.infected,variant.inf))
    return(list(summary_out,timeseries_out,indivudal_out,date.vaccinated,edges))
} # end function definition

