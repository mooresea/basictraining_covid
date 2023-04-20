###
### State-level time series of variants
###
library(tidyverse)
library(jsonlite)
library(data.table)
setwd("d:/covid/basictraining_covid/data")

###Variants data from covariants.org (compiled from GSAID)
###data pulled from: https://github.com/hodcroftlab/covariants/tree/master/cluster_tables
var_dat=fromJSON(txt="USAClusters_data.json")

var_state_l=list()
for(i in 1:length(var_dat$countries)){
  i_df=as.data.frame(var_dat$countries[[i]])
  i_df$state=names(var_dat$countries)[i]
  var_state_l[[i]]=i_df
}
var_st=rbindlist(var_state_l,fill=T)

##Variant proportions
var_st$alpha_pct=var_st$X20I..Alpha..V1./var_st$total_sequences
var_st$beta_pct=var_st$X20H..Beta..V2./var_st$total_sequences
var_st$beta_pct=ifelse(is.na(var_st$beta_pct),0,var_st$beta_pct)
var_st$gamma_pct=var_st$X20J..Gamma..V3./var_st$total_sequences
var_st$delta_pct=(var_st$X21A..Delta.+var_st$X21I..Delta.+var_st$X21J..Delta.)/var_st$total_sequences
var_st$omicron_pct=(var_st$X21K..Omicron.+var_st$X21L..Omicron.)/var_st$total_sequences
var_st$other_pct=1-var_st$alpha_pct-var_st$beta_pct-var_st$delta_pct-var_st$gamma_pct-var_st$omicron_pct
var_st$week=as.Date(var_st$week)

##Clean up NAs
var_st$state[which(is.na(var_st$other_pct))]
var_st$state[which(is.na(var_st$alpha_pct))]
var_st$state[which(is.na(var_st$beta_pct))]
var_st$state[which(is.na(var_st$gamma_pct))]
var_st$state[which(is.na(var_st$delta_pct))]

var_usa=var_st %>% group_by(week) %>% summarize(alpha=mean(alpha_pct,na.rm=T),
                                                      beta=mean(beta_pct,na.rm=T),
                                                      gamma=mean(gamma_pct,na.rm=T),
                                                      delta=mean(delta_pct,na.rm=T),
                                                      omicron=mean(omicron_pct,na.rm=T),
                                                      other=mean(other_pct,na.rm=T))

plot(var_usa$week,var_usa$other,typ="l")
lines(var_usa$week,var_usa$alpha,col="red")
lines(var_usa$week,var_usa$beta,col="blue")
lines(var_usa$week,var_usa$delta,col="green")
lines(var_usa$week,var_usa$gamma,col="purple")
lines(var_usa$week,var_usa$omicron,col="black",lty=2)

var_st_dat=var_st %>% select(week,state,total_sequences,alpha_pct,beta_pct,gamma_pct,delta_pct,omicron_pct,other_pct)

# state_plot="Florida"
# plot(var_st_dat$week[var_st_dat$state==state_plot],var_st_dat$other_pct[var_st_dat$state==state_plot],typ="l")
# lines(var_st_dat$week[var_st_dat$state==state_plot],var_st_dat$alpha_pct[var_st_dat$state==state_plot],col="red")
# lines(var_st_dat$week[var_st_dat$state==state_plot],var_st_dat$beta_pct[var_st_dat$state==state_plot],col="blue")
# lines(var_st_dat$week[var_st_dat$state==state_plot],var_st_dat$delta_pct[var_st_dat$state==state_plot],col="green")
# lines(var_st_dat$week[var_st_dat$state==state_plot],var_st_dat$gamma_pct[var_st_dat$state==state_plot],col="purple")

write_csv(var_st_dat,file="covid_variants_state_time_series.csv")

###
### Military recruit-specific variant fractions
###

##Recruits by state
rec_dat=read_csv(file="recruits_by_state_2018.csv")
all(rec_dat$State %in% unique(var_st_dat$state))
unique(rec_dat$State)[which(!unique(rec_dat$State) %in% unique(var_st_dat$state))]
rec_dat$State[rec_dat$State=="District of Columbia"]="Washington DC"


var_st_dat_mr=var_st_dat %>% left_join(rec_dat,by=c("state"="State"))

var_st_dat_mr$alpha_recruit_pct=var_st_dat_mr$alpha_pct*var_st_dat_mr$recruit_pct/100
var_st_dat_mr$gamma_recruit_pct=var_st_dat_mr$gamma_pct*var_st_dat_mr$recruit_pct/100
var_st_dat_mr$delta_recruit_pct=var_st_dat_mr$delta_pct*var_st_dat_mr$recruit_pct/100
var_st_dat_mr$omicron_recruit_pct=var_st_dat_mr$omicron_pct*var_st_dat_mr$recruit_pct/100
var_st_dat_mr$other_recruit_pct=(var_st_dat_mr$beta_pct+var_st_dat_mr$other_pct)*var_st_dat_mr$recruit_pct/100

var_mr_dat=var_st_dat_mr %>% group_by(week) %>% 
  summarize(alpha_pct=sum(alpha_recruit_pct,na.rm=T),
            gamma_pct=sum(gamma_recruit_pct,na.rm=T),
            delta_pct=sum(delta_recruit_pct,na.rm=T),
            omicron_pct=sum(omicron_recruit_pct,na.rm=T),
            other_pct=sum(other_recruit_pct,na.rm=T))
var_mr_dat$total_pct=var_mr_dat$alpha_pct+var_mr_dat$gamma_pct+var_mr_dat$delta_pct+var_mr_dat$other_pct+var_mr_dat$omicron_pct
var_mr_dat$alpha_pct=var_mr_dat$alpha_pct/var_mr_dat$total_pct
var_mr_dat$gamma_pct=var_mr_dat$gamma_pct/var_mr_dat$total_pct
var_mr_dat$delta_pct=var_mr_dat$delta_pct/var_mr_dat$total_pct
var_mr_dat$omicron_pct=var_mr_dat$omicron_pct/var_mr_dat$total_pct
var_mr_dat$other_pct=var_mr_dat$other_pct/var_mr_dat$total_pct

##Remove final week with limited data
#var_mr_dat$week

plot(var_mr_dat$week,var_mr_dat$total_pct,ylim=c(0,1))
lines(var_mr_dat$week,var_mr_dat$alpha_pct,col="red")
lines(var_mr_dat$week,var_mr_dat$gamma_pct,col="purple")
lines(var_mr_dat$week,var_mr_dat$delta_pct,col="green")
lines(var_mr_dat$week,var_mr_dat$omicron_pct,col="orange")
lines(var_mr_dat$week,var_mr_dat$other_pct,lty=2)

write_csv(var_mr_dat,file="covid_variants_recruits.csv")