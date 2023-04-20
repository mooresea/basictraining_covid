#### Set default seed
##set.seed(123456)
library(tidyverse)

load_network_params = function(site_name, params_file){
    params_df_calib = read.csv(params_file)
    if(site_name == 'full'){
        network_file = '../output/network.RData'
        R0.mag = mean(params_df_calib$R0_mag)
        initial.infected = mean(params_df_calib$initial_inf)
        genint_scale = mean(params_df_calib$genint_scale)
        genint_shape = mean(params_df_calib$genint_shape)
        dispersion_val= params_df_calib$dispersion_mag[params_df_calib$base_name == "base"]
        #R0_adj_factor = params_df_calib$R0_adj_factor[params_df_calib$base_name == "base"]
    }else{        
        network_file = sprintf('../output/network_%s.RData', training_site)
        R0.mag = params_df_calib$R0_mag[params_df_calib$base_name == site_name]
        initial.infected = params_df_calib$initial_inf[params_df_calib$base_name == site_name]
        genint_scale = params_df_calib$genint_scale[params_df_calib$base_name == site_name]
        genint_shape = params_df_calib$genint_shape[params_df_calib$base_name == site_name]
        dispersion_val= params_df_calib$dispersion_mag[params_df_calib$base_name == site_name]
        #R0_adj_factor = params_df_calib$R0_adj_factor[params_df_calib$base_name == site_name]
    }
    R0_adj_factor=1
    return(list(R0.mag = R0.mag, initial.infected = initial.infected, network_file = network_file,
                genint_scale = genint_scale, genint_shape = genint_shape,dispersion_val=dispersion_val,R0_adj_factor=R0_adj_factor))
}
#### load libraries
library('scam')

## Check if training_site is specified, else, specified it as 'full'
if(!exists('training_site')){
    training_site = 'full'
}
if(!(training_site %in% c('FB', 'FLW','base'))){
    training_site = 'full'
}


#### Load data
calib_params_list = load_network_params(training_site, '../output/calibrated_parameters_forecast.csv')

load(calib_params_list$network_file)
## parameter values for infection processes
R0.mag = calib_params_list$R0.mag
initial.infected = calib_params_list$initial.infected
R0.tim = 5.2 ## 7.8
symp.mag = 0.31 #0.57 #0.31- Changed this based on data from Shilpa (changed back to higher value for now)
symp.tim = 5.8
symp.mean = 10
asymp.adjust = 0.35 #0.8 #Relative infectiousness of asymptomatic individuals - Buitrago-Garcia et al.
initial.immune = 0.026
importation = c(rep(0, numrecruits), rep(0.01/84, numdrillserg + numStaff)) # 1% chance of infection over 84 day period

#R0.prob = R0.mag / mean(total.contacts) / (symp.mag + (1 - symp.mag) * asymp.adjust)
R0.prob = R0.mag / mean(average.contacts) / (symp.mag + (1 - symp.mag) * asymp.adjust)

#Adjustment factor that can be used to modulate R0 between initial cocoon/quarantine and training periods
R0_adj_factor=1
## parameter values for control processes
## testing
BOOL_testDaily = FALSE
BOOL_testOnArrival = 1
testDelayQuarantine = 3
testReturn = 1
testdates = c(14)
testdates_type = c('pcr')

## Staff routine testing
testStaffFreq = 0
testStaffType = 'antigen'

## testing upon exit
## Not necessary, just set BOOL_clinicalrelease = 0

## isolation/quarantine and contact tracing
BOOL_clinicalrelease = 1
isolate.max = 10
isolate.min = 10
isolate.length = 10
isolate.nosymp = 0
isolate.effectiveness = 0.99
quarantine.max = 10
quarantine.effectiveness = 0.99
quarantine.contacts = 0

## masks
mask.protection = 0.3
mask.compliance = 1.0 #0.2 - Use compliance.mult to adjust compliance downwards based on date
compliance.avg = mask.protection * mask.compliance
compliance.spread = 1e1
compliance = rep(compliance.avg, numrecruits + numdrillserg + numStaff)

mask.final.compliance = 0.1
compliance.final.avg = mask.protection * mask.final.compliance
BOOL_compliance_time = 0
BOOL_compliance_date = 1
##Mask time series
mask_dat=read_csv(file="../data/CovidCast_mask_wearing.csv")
mask_dat=mask_dat %>% filter(regionId=="SC")
mask_dat$value=mask_dat$value/100

## All the individuals have the same compliance, 
## compliance = rbeta(
##   numrecruits+numdrillserg+numStaff,
##   compliance.avg * compliance.spread,
##   (1 - compliance.avg) * compliance.spread
## )

## specify arrival dates on campus
#arrivebase = c(sample(1:1, numrecruits, replace = T), rep(0, numdrillserg + numStaff))

pcr.spec.screen = 0.998
pcr.spec.commercial = 0.998

pcr.sens.screen.target = 0.859
pcr.sens.commercial.target = 0.859

##======================================#
## Incubation period and sens of tests
##======================================#
## Time from symptom onset to transmission
incub.shape =  5.807
incub.scale = 0.948
tost.mu = 0.46
tost.ke = 2.8
k_inc = incub.shape
gamma = 1 / (k_inc * incub.scale)
mu = tost.mu
k_E = tost.ke
k_I = 1
alpha = 1
k_P = k_inc - k_E
C = k_inc * gamma * mu / (alpha * k_P * mu + k_inc * gamma)

t_tost = seq(-23,23,0.05)
fm = alpha * C * (1 - pgamma(-t_tost, shape=k_P, scale=1/(k_inc*gamma)))
fp = C * (1 - pgamma(t_tost, shape=k_I, scale=1/(k_I*mu)))
f = ifelse(t_tost>0,fp,fm)
f = approxfun(t_tost, f)

t_tost = -22:22
f_tost = sapply(t_tost,function(t)integrate(f,t,t+1)[[1]])
f_tost = f_tost / sum(f_tost)
tost = cbind(t_tost,f_tost)

load('sens_postOnset_extended.RData')
sens_draw = sens_postOnset.median
## sensitivity =
##   rbind(
##     cbind(
##       tost[2:23,1],
##       tost[1:22,2]/tost[22,2]*max(sens_draw)),
##     cbind(
##       1:33,
##       sens_draw[-1]))

if(is.null(calib_params_list$genint_scale)){
    calib_params_list$genint_scale = 5.67
    calib_params_list$genint_shape = 2.89
    
}
genint_shape = calib_params_list$genint_shape
genint_scale = calib_params_list$genint_scale
genint = pweibull(1:21,scale=genint_scale,shape=genint_shape) - pweibull(0:20,scale=genint_scale,shape=genint_shape)
genint = genint / sum(genint)

sensitivity = rbind(
    cbind(
        1:6,
        genint[1:6]/genint[6] * max(sens_draw)),
    cbind(
        (1:(length(sens_draw) - 1))  + 6,
        sens_draw[-1]))



sensitivity_array = sensitivity[,2]
sensitivity_time = sensitivity[,1]

sens.daily = sens_postOnset.median * (pweibull(1:length(sens_postOnset.median),shape=0.9,scale=5.657) -
                                      pweibull(0:(length(sens_postOnset.median)-1),shape=0.9,scale=5.657))
sens.overall.default = sum(sens.daily)

sens.fun = function(sens.overall){
  pmin(1, (sens.overall / sens.overall.default) * sensitivity[,2])
}

pcr.sens.commercial = function(dd){
    ## Day here means day relative to infection
    dd[dd > length(sensitivity_array)] = length(sensitivity_array)
    dd[dd < 1] = 1
    return(sens.fun(pcr.sens.commercial.target)[dd])
}


pcr.sens.screen = function(dd){
    ## Day here means day relative to infection
    dd[dd > length(sensitivity_array)] = length(sensitivity_array)
    dd[dd < 1] = 1
    return(sens.fun(pcr.sens.screen.target)[dd])
}

#########Waning of natural immunity
#For original 
#Based on study of Townshend et al. 21 https://doi.org/10.1016/S2666-5247(21)00219-6
# Fit an exponential decay model of protection for different coronaviruses
sars_cov2_a=-4.795377
sars_cov2_b=-12.14602
sars_cov2_g=0.1298892
imm_decay_rate=1/(1+exp(-(sars_cov2_a+sars_cov2_b*sars_cov2_g)))
imm_decay_start=90 #day after infection at which immune decay starts
imm_total_init=30 #Assume that initial immune response to infection prevents re-infection from any variant in 1st 30 days

###Vaccination parameters

##Vaccination of recruits on arrival
include_vacc=T
include_booster=F
vacc_start_date=as.Date("2021-05-22")
vacc_rate=0.25 ##Fraction of unvaccinated recruits that agree to vaccination
vacc_start=3 #day on which these recruits are vaccinated - just a guess at this point
booster_vacc_rate=0 #vacc_rate

#Protection against infection 
vacc_1_dose_eff=0.368 #95%CI: .332-0.402 #Chemaitelly et al. 21 (NEJM)
vacc_2_dose_eff=0.775 #95%CI: 0.764-0.786
vacc_3_dose_eff=vacc_2_dose_eff #Assume booster returns protection to 1-mo post-2nd dose
#vacc_decay_start=16 #Peak was measured at 1-mo post-2nd dose, we have already built in 14 days to full protection
rel_symp_frac=0.695 # This is the fraction of 20-29 year olds who have clinical case relative to pop average (for US) based on Davies et al
##Probability that a vaccinated individual is symptomatic given that they are infected
## Calculation is prob of clinical case if infected among vaccinated individuals (population average) * relative symptomatic fraction for 20-29
vacc.symp.adj=rel_symp_frac*(1-0.945)/(1-vacc_2_dose_eff) 
  
##Waning of vaccine protection - based on data from Chemaitelly et al. 
ve.dates=(1:7)*30 #months of measurement
ve.est=c(vacc_2_dose_eff,.732,.696,.517,.225,.173,.223)
ve.df=data.frame(date=ve.dates,ve=ve.est)
#estimate logistic decay of vaccine protection
glm.ve=suppressWarnings(glm(ve~date,data=ve.df,family="binomial"))
nls.ve<-nls(ve~phi1/(1+exp(-(phi2+phi3*date))),
            start=list(phi1=ve.est[1],phi2=coef(glm.ve)[1],phi3=coef(glm.ve)[2]),data=ve.df)
ve_peak=coef(nls.ve)["phi1"]
ve_p1=coef(nls.ve)["phi2"]
ve_dr=coef(nls.ve)["phi3"]

###
### Variant specific parameters
###
#Variant time series
var_dat=read_csv(file="../data/covid_variants_recruits.csv")
var_p=read_csv(file="../data/variant_parameters.csv")
#Transmission multiplier relative to original 
#From Campbell et al. 21 - https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2021.26.24.2100509
#var_param_dat=read_csv(file="../data/variant_param_values.csv")
# trans_alpha=1.29
# trans_alpha_l=1.24
# trans_alpha_h=1.33
# trans_gamma=1.38
# trans_gamma_l=1.29
# trans_gamma_h=1.48
# trans_delta=1.97
# trans_delta_l=1.76
# trans_delta_h=2.17
# trans_omicron=NA