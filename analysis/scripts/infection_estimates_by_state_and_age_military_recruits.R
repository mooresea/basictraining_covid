##
## Generate age-specific infection time series from COVID-Estim estimates
##
library(tidyverse)
library(ggpubr)
library(data.table)
setwd("d:/covid/basictraining_covid")
st_pops=read_csv(file="data/state_pops_2020.csv")

#Age-specific symptomatic (clinical) fractions
age_symp=read_csv(file="data/davies_symptomatic_by_age.csv")
age_symp$age_class=substr(age_symp$age_group,2,3)

## Age-specific infection trends at national-level
age_dat=read_csv(file="data/national_age-specific_incidence_RR_TS.csv")
age_dat = age_dat %>% select(date,age_class,age_group.y,age_pct,incid_OR)
age_dat = rename(age_dat,age_group=age_group.y)
age_dat$date=as.Date(age_dat$date,format="%m/%d/%Y")

###Infection estimates
###########
# Covid-estim
###########
inf_est=read_csv(file="data/covid-estim_infection_estimates.csv")
all(unique(inf_est$state) %in% unique(st_pops$state))
which(!unique(inf_est$state) %in% unique(st_pops$state))
inf_est=rename(inf_est, state_name=state)
inf_est$date=as.Date(inf_est$date,format="%m/%d/%Y")
inf_est=inf_est %>% 
  group_by(state_name) %>%
  mutate(state=unlist(lapply(state_name,function(x) {st_pops$code[st_pops$state==x]})),
         pop=unlist(lapply(state_name,function(x) {st_pops$pop2020[st_pops$state==x]}))) %>%
  ungroup()

##Get number of infections for each day and state
inf_est=inf_est %>% select(state,pop,date,cases.fitted,cases.fitted.hi,cases.fitted.lo,cum.incidence,infections,infections.hi,infections.lo,
                           pop.infectiousness,pop.infectiousness.hi,pop.infectiousness.lo)


##Total infections
inf_sum = inf_est %>% group_by(date) %>% summarize(cases=sum(cases.fitted,na.rm = T),
                                                   cases.lo=sum(cases.fitted.lo,na.rm=T),
                                                   cases.hi=sum(cases.fitted.hi,na.rm=T),
                                                   infections=sum(infections,na.rm=T),
                                                   infections.lo=sum(infections.lo,na.rm=T),
                                                   infections.hi=sum(infections.hi,na.rm=T))
##
## Back to age-specific infection estimates
##
##Replicate the last day in the age dataset up until the end of the infections dataset
#Assumes that age-specific ORs remain constant
max(inf_est$date)
final_age=age_dat[which(age_dat$date==max(age_dat$date)),]
add_dates=seq((max(age_dat$date)+1),max(inf_est$date),by=1)
for(dd in 1:length(add_dates)){
  dd_dat=final_age
  dd_dat$date=add_dates[dd]
  age_dat=rbind(age_dat,dd_dat)
}
##
## Create dataset with age classes for each date and state
##
inf_est_age=inf_est %>% left_join(age_dat,by="date")
inf_est_age$age_infections=inf_est_age$infections*inf_est_age$age_pct*inf_est_age$incid_OR
inf_est_age$age_infections_lo=inf_est_age$infections.lo*inf_est_age$age_pct*inf_est_age$incid_OR
inf_est_age$age_infections_hi=inf_est_age$infections.hi*inf_est_age$age_pct*inf_est_age$incid_OR
#inf_est_age$age_infectiousness=inf_est_age$pop.infectiousness*inf_est_age$age_pct*inf_est_age$incid_OR
inf_est_age$age_cases=inf_est_age$cases.fitted*inf_est_age$age_pct*inf_est_age$incid_OR




##Check that sum of age-specific infections  match total infections
inf_age_sum=inf_est_age %>% group_by(state,date) %>% summarize(infections_sum=sum(age_infections))
inf_age_sum=inf_age_sum %>% left_join(inf_est,by=c("state","date"))
inf_age_sum$diff=inf_age_sum$infections_sum-inf_age_sum$infections
hist(inf_age_sum$diff)

##Get pop by age class
inf_est_age$pop_age=inf_est_age$pop*inf_est_age$age_pct
inf_est_age=inf_est_age %>% select(date,state,age_class,age_group,pop_age,age_cases,age_infections,age_infections_lo,age_infections_hi)
inf_est_age=rename(inf_est_age,
                   cases=age_cases,
                   infections=age_infections,
                   infections_lo=age_infections_lo,
                   infections_hi=age_infections_hi)

##
##Military recruit age
##
inf_mr=inf_est_age[inf_est_age$age_class=="20",]
inf_mr$case_pct=inf_mr$cases/inf_mr$pop_age
inf_mr$infection_pct=inf_mr$infections/inf_mr$pop_age
inf_mr$infection_pct_lo=inf_mr$infections_lo/inf_mr$pop_age
inf_mr$infection_pct_hi=inf_mr$infections_hi/inf_mr$pop_age

ggplot(inf_est_age[inf_est_age$state=="FL",],aes(x=date,y=infections,color=age_group))+geom_point()
ggplot(inf_mr[inf_mr$state=="FL",],aes(x=date,y=infections,color=age_group))+geom_point()
plot(inf_mr$date[inf_mr$state=="FL"],inf_mr$infections[inf_mr$state=="FL"])
lines(inf_mr$date[inf_mr$state=="FL"],inf_mr$infections_hi[inf_mr$state=="FL"])

### Daily infections by state - use to calculate how many will test positive upon arrival
write_csv(inf_mr,file="data/daily_estimated_infections_recruits_by_state.csv")

### Aggregate infections by state to total infections 

##Recruits by state
rec_dat=read_csv(file="data/recruits_by_state_2018.csv")

inf_mr_total=inf_mr %>% inner_join(rec_dat,by=c("state"="code"))
inf_mr_total$infections_frac=inf_mr_total$infections*inf_mr_total$recruit_pct/100
inf_mr_total$infections_lo_frac=inf_mr_total$infections_lo*inf_mr_total$recruit_pct/100
inf_mr_total$infections_hi_frac=inf_mr_total$infections_hi*inf_mr_total$recruit_pct/100
inf_mr_total$pop_age_frac=inf_mr_total$pop_age*inf_mr_total$recruit_pct/100
inf_mr_total$cases_frac=inf_mr_total$cases*inf_mr_total$recruit_pct/100

inf_mr_dat=inf_mr_total %>% group_by(date) %>% 
  summarize(pop_age=sum(pop_age),
            pop_age_frac=sum(pop_age_frac),
            cases=sum(cases_frac),
            infections=sum(infections_frac),
            infections_lo=sum(infections_lo_frac),
            infections_hi=sum(infections_hi_frac))
inf_mr_dat$infections_pct=inf_mr_dat$infections/inf_mr_dat$pop_age_frac
inf_mr_dat$infections_lo_pct=inf_mr_dat$infections_lo/inf_mr_dat$pop_age_frac
inf_mr_dat$infections_hi_pct=inf_mr_dat$infections_hi/inf_mr_dat$pop_age_frac

### Daily infections by state - use to calculate how many will test positive upon arrival
write_csv(inf_mr_dat,file="data/daily_estimated_infections_recruits_aggregated.csv")

##
##Incorporate COVI-19 Forecasting Hub forecasts
##https://github.com/reichlab/covid19-forecast-hub/tree/master/data-processed/COVIDhub-4_week_ensemble
##
# case.lag=which.max(inf_mr_dat$cases)-which.max(inf_mr_dat$infections) #Lag from infection to case report based on omicron peak
# 
# inf_mr_dat$rep_frac=inf_mr_dat$cases/inf_mr_dat$infections

for_dat=read_csv(file="data/2022-02-28-COVIDhub-4_week_ensemble.csv")
for_dat = for_dat %>% filter(location=="US",quantile==0.5)
for_dat=for_dat[grepl("case",for_dat$target),]
for_dat$target_start_date=for_dat$target_end_date-6
for_dat$daily_cases=for_dat$value/7

##Add forecasted cases to COVID-estim estimates
case.lag.us=which.max(inf_sum$cases)-which.max(inf_sum$infections)
inf_sum$rep_frac=c(rep(inf_sum$cases[case.lag.us]/inf_sum$infections[1],case.lag.us-1),
                   inf_sum$cases[case.lag.us:nrow(inf_sum)]/inf_sum$infections[1:(nrow(inf_sum)-case.lag.us+1)])
inf_sum$rep_frac_lo=c(rep(inf_sum$cases[case.lag.us]/inf_sum$infections.lo[1],case.lag.us-1),
                   inf_sum$cases[case.lag.us:nrow(inf_sum)]/inf_sum$infections.lo[1:(nrow(inf_sum)-case.lag.us+1)])
inf_sum$rep_frac_hi=c(rep(inf_sum$cases[case.lag.us]/inf_sum$infections.hi[1],case.lag.us-1),
                   inf_sum$cases[case.lag.us:nrow(inf_sum)]/inf_sum$infections.hi[1:(nrow(inf_sum)-case.lag.us+1)])

inf_sum_fcast=data.frame(date=seq(inf_sum$date[nrow(inf_sum)]+1,max(for_dat$target_end_date),by=1),
                         cases=0,cases.lo=0,cases.hi=0,infections=NA,infections.lo=NA,infections.hi=NA,
                         rep_frac=mean(inf_sum$rep_frac[(nrow(inf_sum)-6):nrow(inf_sum)]),
                         rep_frac_lo=mean(inf_sum$rep_frac_lo[(nrow(inf_sum)-6):nrow(inf_sum)]),
                         rep_frac_hi=mean(inf_sum$rep_frac_hi[(nrow(inf_sum)-6):nrow(inf_sum)]))
                         
for(i in 1:nrow(for_dat)){
  iinds=which(inf_sum_fcast$date >= for_dat$target_start_date[i] & inf_sum_fcast$date<=for_dat$target_end_date[i])
  inf_sum_fcast$cases[iinds]=for_dat$daily_cases[i]
}
inf_sum_fcast$infections[1:(nrow(inf_sum_fcast)-case.lag.us+1)]=inf_sum_fcast$cases[case.lag.us:nrow(inf_sum_fcast)]/inf_sum_fcast$rep_frac[1]
inf_sum_fcast$infections.lo[1:(nrow(inf_sum_fcast)-case.lag.us+1)]=inf_sum_fcast$cases[case.lag.us:nrow(inf_sum_fcast)]/inf_sum_fcast$rep_frac_lo[1]
inf_sum_fcast$infections.hi[1:(nrow(inf_sum_fcast)-case.lag.us+1)]=inf_sum_fcast$cases[case.lag.us:nrow(inf_sum_fcast)]/inf_sum_fcast$rep_frac_hi[1]

##Remove end of time series with NA infections
inf_sum_fcast=inf_sum_fcast[!is.na(inf_sum_fcast$infections),]
inf_sum$forecast=F
inf_sum_fcast$forecast=T
inf_sum=rbind(inf_sum,inf_sum_fcast)

##Now combine recruit_infections with total_us_infections
inf_mr_dat_fcast=inf_mr_dat %>% right_join(inf_sum,by="date")
inf_mr_dat_fcast=inf_mr_dat_fcast[which(!is.na(inf_mr_dat_fcast$infections.x)|inf_mr_dat_fcast$forecast==T),]
inf_mr_dat_fcast$rec_inf_pct=inf_mr_dat_fcast$infections.x/inf_mr_dat_fcast$infections.y
recent_inf_pct=mean(tail(inf_mr_dat_fcast$rec_inf_pct[which(inf_mr_dat_fcast$forecast==F)],7))
inf_mr_dat_fcast$infections.x[which(inf_mr_dat_fcast$forecast==T)]=inf_mr_dat_fcast$infections.y[which(inf_mr_dat_fcast$forecast==T)]*recent_inf_pct
inf_mr_dat_fcast$infections_lo[which(inf_mr_dat_fcast$forecast==T)]=inf_mr_dat_fcast$infections.lo[which(inf_mr_dat_fcast$forecast==T)]*recent_inf_pct
inf_mr_dat_fcast$infections_hi[which(inf_mr_dat_fcast$forecast==T)]=inf_mr_dat_fcast$infections.hi[which(inf_mr_dat_fcast$forecast==T)]*recent_inf_pct

pop_frac=inf_mr_dat_fcast$pop_age_frac[nrow(inf_mr_dat)]
inf_mr_dat_fcast$infections_pct=inf_mr_dat_fcast$infections.x/pop_frac
inf_mr_dat_fcast$infections_lo_pct=inf_mr_dat_fcast$infections_lo/pop_frac
inf_mr_dat_fcast$infections_hi_pct=inf_mr_dat_fcast$infections_hi/pop_frac


### Daily infections by state - use to calculate how many will test positive upon arrival
write_csv(inf_mr_dat_fcast,file="data/daily_estimated_infections_recruits_aggregated_4wk_forecast.csv")
