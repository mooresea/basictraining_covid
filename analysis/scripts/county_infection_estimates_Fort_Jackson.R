##
## Generate age-specific infection time series from COVID-Estim estimates
##
library(tidyverse)
library(ggpubr)
library(data.table)
setwd("d:/covid/basictraining_covid")
st_pops=read_csv(file="data/state_pops_2020.csv")


###Infection estimates
###########
# Covid-estim
###########
inf_est_c=read_csv(file="data/covid-estim_county_infection_estimates.csv") 
inf_est_c=inf_est_c %>% filter(fips=="45079")
inf_est_c$date=as.Date(inf_est_c$date,format="%m/%d/%Y")

##Get number of infections for each day
inf_est_c=inf_est_c %>% select(date,cases.fitted,cum.incidence,infections,infections.lo,infections.hi)
inf_est_c$pop=415759
inf_est_c$pct_inf=inf_est_c$infections/inf_est_c$pop
inf_est_c$pct_inf_lo=inf_est_c$infections.lo/inf_est_c$pop
inf_est_c$pct_inf_hi=inf_est_c$infections.hi/inf_est_c$pop

write_csv(inf_est_c,file="output/daily_estimated_infections_Fort_Jackson.csv")