##
## State-level vaccination rates
##
library(tidyverse)
library(lubridate)

setwd("d:/covid/basictraining_covid")

##Recruits by state
rec_dat=read_csv(file="data/recruits_by_state_2018.csv")

##
st_pops=read_csv(file="data/state_pops_2020.csv")
us_pop=read_csv(file="data/us_census_pop_2020.csv",skip_empty_rows = T)
# us_pop$age_grp=ifelse(us_pop$age<5,"<5 yrs",
#                       ifelse(us_pop$age>=75,"75+ yrs",
#                              ifelse(us_pop$age>=65,"65-74 yrs",
#                                     ifelse(us_pop$age>=50,"50-64 yrs",
#                                            ifelse(us_pop$age>=40,"40-49 yrs",
#                                                   ifelse(us_pop$age>=25,"25-39 yrs",
#                                                          ifelse(us_pop$age>=18,"18-24 yrs",
#                                                                 ifelse(us_pop$age>=12,"12-17 yrs","5-11 yrs"))))))))
us_pop$age_grp=ifelse(us_pop$age<5,"grp_0_4",
                      ifelse(us_pop$age>=75,"grp_75",
                             ifelse(us_pop$age>=65,"grp_65_74",
                                    ifelse(us_pop$age>=50,"grp_50_64",
                                           ifelse(us_pop$age>=40,"grp_40_49",
                                                  ifelse(us_pop$age>=25,"grp_25_39",
                                                         ifelse(us_pop$age>=18,"grp_18_24",
                                                                ifelse(us_pop$age>=12,"grp_12_17","grp_5_11"))))))))
pop_age=us_pop %>% group_by(age_grp) %>% summarise(pop=sum(Middle)*1000)
pop_age=pop_age[!is.na(pop_age$pop),]
pop_age$age_pct=pop_age$pop/sum(pop_age$pop)

##National age-specific vaccination coverage trends
vacc_age=read_csv(file="data/COVID-19_Vaccination_and_Case_Trends_by_Age_Group__United_States.csv")
vacc_age$Date=as.Date(vacc_age$`Date Administered`,format="%m/%d/%Y")
old_age_grps=sort(unique(vacc_age$AgeGroupVacc))
new_age_grps=sort(unique(pop_age$age_grp))[-1] #No 0-4 in vaccination data
vacc_age$age_grp=new_age_grps[unlist(lapply(vacc_age$AgeGroupVacc, function(x) which(old_age_grps==x)))]

vacc_age=vacc_age %>% left_join(pop_age,by=c("age_grp"="age_grp"))
vacc_age$num_dose1=vacc_age$Administered_Dose1_pct_agegroup*vacc_age$pop
vacc_age$num_dose2=vacc_age$Series_Complete_Pop_pct_agegroup*vacc_age$pop

#vacc_age_18plus=vacc_age %>% filter(!age_grp %in% c("grp_0_4","grp_12_17","grp_5_11"))
##replace vacc_age with age_18plus -> issue is state specific data on 18+ only available starting in mid-Feb
vacc_avg=vacc_age %>% group_by(Date) %>% summarize(pop=sum(pop),
                                                   num_dose1=sum(num_dose1),
                                                   num_dose2=sum(num_dose2))
vacc_avg$pct_dose1=vacc_avg$num_dose1/vacc_avg$pop
vacc_avg$pct_dose2=vacc_avg$num_dose2/vacc_avg$pop

#Military age-adjusted
# vacc_18=vacc_age %>% filter(age_grp=="grp_18_24") %>% select(Date,Administered_Dose1_pct_agegroup,
#                                                              Series_Complete_Pop_pct_agegroup,pop,age_pct,num_dose1,num_dose2)
vacc_18=vacc_age %>% filter(age_grp=="grp_18_24") %>% select(Date,Administered_Dose1_pct_agegroup,
                                                                   Series_Complete_Pop_pct_agegroup,age_pct)
names(vacc_18)[2:3]=c("pct_1st_dose_18a","pct_2nd_dose_18a")
# vacc_18=vacc_18 %>% left_join(vacc_avg,by="Date")
# vacc_18$rel_dose1=vacc_18$pct_1st_dose_18a/vacc_18$pct_dose1
# vacc_18$rel_dose2=vacc_18$pct_2nd_dose_18a/vacc_18$pct_dose2
  
## State-level vaccination rates over time
vacc_dat=read_csv(file="data/COVID-19_Vaccination_Trends_in_the_United_States_National_and_Jurisdictional.csv") %>% 
            filter(date_type=="Admin")
vacc_dat$Date=as.Date(vacc_dat$Date,format="%m/%d/%Y")
#us_vacc=vacc_dat %>% filter(Location=="US") %>% select(Date,Admin_Per_100K,Admin_Per_100k_18Plus) ##Older vaccine file details
us_vacc=vacc_dat %>% filter(Location=="US") %>% select(Date,Admin_Dose_1_Cumulative,
                                                                   Series_Complete_Cumulative,
                                                                   Booster_Cumulative) 
us_vacc$pop_total=sum(st_pops$pop2020)
us_vacc$total_dose_1=us_vacc$Admin_Dose_1_Cumulative/us_vacc$pop_total
us_vacc$total_dose_2=us_vacc$Series_Complete_Cumulative/us_vacc$pop_total
us_vacc$total_dose_3=us_vacc$Booster_Cumulative/us_vacc$pop_total
us_vacc=us_vacc %>% left_join(vacc_18,by="Date")
us_vacc$rel_dose1_age18=us_vacc$pct_1st_dose_18a/us_vacc$total_dose_1
us_vacc$rel_dose2_age18=us_vacc$pct_2nd_dose_18a/us_vacc$total_dose_2

#us_vacc$rel_dose1_age18=us_vacc$pct_1st_dose_18a/us_vacc$total_dose_1
st_vacc=vacc_dat %>% filter(Location %in% rec_dat$code) %>% select(Date,Location,Admin_Dose_1_Cumulative,
                                                                   Series_Complete_Cumulative,
                                                                   Booster_Cumulative) 

st_vacc=st_vacc %>% left_join(st_pops,by=c("Location"="code"))
st_vacc$pct_dose_1=st_vacc$Admin_Dose_1_Cumulative/st_vacc$pop2020
st_vacc$pct_dose_2=st_vacc$Series_Complete_Cumulative/st_vacc$pop2020
st_vacc$pct_dose_3=st_vacc$Booster_Cumulative/st_vacc$pop2020
st_vacc=st_vacc %>% left_join(us_vacc,by="Date")

##
##Method #1 - Adjust 18-24 year old vaccination coverage up/down based on # of % doses administered at state-level vs. national average
##
st_vacc1=st_vacc
st_vacc1$pct_dose_1_18age=st_vacc1$pct_1st_dose_18a*(st_vacc1$pct_dose_1/st_vacc1$total_dose_1)
st_vacc1$pct_dose_2_18age=st_vacc1$pct_2nd_dose_18a*(st_vacc1$pct_dose_2/st_vacc1$total_dose_2)
#st_vacc1$pct_dose_3_18age=st_vacc1$pct_3rd_dose_18a*(st_vacc1$pct_dose_3/st_vacc1$total_dose_3)
st_vacc1=st_vacc1 %>% left_join(rec_dat,by=c("Location"="code"))
##
## Method #2 - Adjust state level vaccination coverage up/down based on %of 18-24 year olds vaccinated relative to overall population (16+)
##
st_vacc2=st_vacc
st_vacc2$pct_dose_1_18age=st_vacc2$pct_dose_1*st_vacc2$rel_dose1_age18
st_vacc2$pct_dose_2_18age=st_vacc2$pct_dose_2*st_vacc2$rel_dose2_age18
st_vacc2$pct_dose_3_18age=st_vacc2$pct_dose_3*st_vacc2$rel_dose2_age18 #Assume booster uptake at same relative pct as 2nd dose
st_vacc2=st_vacc2 %>% left_join(rec_dat,by=c("Location"="code"))

##Get recruit-weighted national averages
st_vacc1$pct_dose_1_frac=st_vacc1$pct_dose_1_18age*st_vacc1$recruit_pct/100
st_vacc1$pct_dose_2_frac=st_vacc1$pct_dose_2_18age*st_vacc1$recruit_pct/100

st_vacc2$pct_dose_1_frac=st_vacc2$pct_dose_1_18age*st_vacc2$recruit_pct/100
st_vacc2$pct_dose_2_frac=st_vacc2$pct_dose_2_18age*st_vacc2$recruit_pct/100
st_vacc2$pct_dose_3_frac=st_vacc2$pct_dose_3_18age*st_vacc2$recruit_pct/100

recruit_vacc1=st_vacc1 %>% group_by(Date) %>%
                summarize(pct_dose_1_cumulative=sum(pct_dose_1_frac),
                          pct_dose_2_cumulative=sum(pct_dose_2_frac))

recruit_vacc2=st_vacc2 %>% group_by(Date) %>%
  summarize(pct_dose_1_cumulative=sum(pct_dose_1_frac),
            pct_dose_2_cumulative=sum(pct_dose_2_frac),
            pct_dose_3_cumulative=sum(pct_dose_3_frac))

recruit_vacc2$pct_dose_1=c(0,diff(recruit_vacc2$pct_dose_1_cumulative,lag = 1))
recruit_vacc2$pct_dose_2=c(0,diff(recruit_vacc2$pct_dose_2_cumulative,lag = 1))
recruit_vacc2$pct_dose_3=c(0,diff(recruit_vacc2$pct_dose_3_cumulative,lag = 1))
recruit_vacc2=recruit_vacc2[!is.na(recruit_vacc2$pct_dose_1),]
recruit_vacc2$pct_dose_1=ifelse(recruit_vacc2$pct_dose_1<0,0,recruit_vacc2$pct_dose_1)
recruit_vacc2$pct_dose_2=ifelse(recruit_vacc2$pct_dose_2<0,0,recruit_vacc2$pct_dose_2)
recruit_vacc2$pct_dose_3=ifelse(recruit_vacc2$pct_dose_3<0,0,recruit_vacc2$pct_dose_3)
write_csv(recruit_vacc2,file="data/military_recruit_vaccination_timeseries.csv")

###
### Get vaccination rates for Richland County, SC (Fort Jackson)
###
vacc_dat_c=read_csv(file="data/COVID-19_Vaccinations_in_the_United_States_County.csv")
vacc_dat_c$Date=as.Date(vacc_dat_c$Date,format="%m/%d/%Y")
vacc_dat_c=vacc_dat_c %>% filter(FIPS=="45079")

##Use 18+ measures
vacc_dat_c=vacc_dat_c %>% select(Date,Administered_Dose1_Pop_Pct,Series_Complete_18PlusPop_Pct,Booster_Doses_18Plus_Vax_Pct)
names(vacc_dat_c)[2:4]=c("pct_dose_1_cumulative","pct_dose_2_cumulative","pct_dose_3_cumulative")
vacc_dat_c=vacc_dat_c[order(vacc_dat_c$Date),]
vacc_dat_c$pct_dose_1_cumulative=vacc_dat_c$pct_dose_1_cumulative/100
vacc_dat_c$pct_dose_2_cumulative=vacc_dat_c$pct_dose_2_cumulative/100
vacc_dat_c$pct_dose_3_cumulative=vacc_dat_c$pct_dose_3_cumulative/100

##merge with us_vacc to get booster rates before 12-15-21
vacc_dat_c=vacc_dat_c %>% left_join(us_vacc,by="Date")
##get relative booster rate for Richland County on first date data is available
sc_booster_ratio=vacc_dat_c$pct_dose_3_cumulative[min(which(vacc_dat_c$pct_dose_3_cumulative>0))]/vacc_dat_c$total_dose_3[min(which(vacc_dat_c$pct_dose_3_cumulative>0))]
vacc_dat_c$pct_dose_3_cumulative[1:(min(which(vacc_dat_c$pct_dose_3_cumulative>0))-1)]=sc_booster_ratio*vacc_dat_c$total_dose_3[1:(min(which(vacc_dat_c$pct_dose_3_cumulative>0))-1)]

vacc_dat_c$pct_dose_1=c(0,diff(vacc_dat_c$pct_dose_1_cumulative,lag=1))
vacc_dat_c$pct_dose_2=c(0,diff(vacc_dat_c$pct_dose_2_cumulative,lag=1))
vacc_dat_c$pct_dose_3=c(0,diff(vacc_dat_c$pct_dose_3_cumulative,lag=1))

write_csv(vacc_dat_c,file="data/Richland_SC_vaccination_timeseries.csv")
