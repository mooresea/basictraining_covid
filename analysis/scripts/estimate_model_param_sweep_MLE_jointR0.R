##==================================================#
## load function to simulate outbreaks--------------
##==================================================#
library(reshape2)
library(scam)
library(mgcv)
library(MASS)
library(tidyverse)
hiAsymp=F
if(hiAsymp){
  output_dir = '../output_experiments_joint/single_R0_hiAsymp'
}else{
  output_dir = '../output_experiments_joint/single_R0_loAsymp'
}
##LL function
dLogBetaBinomial<-function(x,n,alpha,beta){
  
  ##Error checking
  if(any(x < 0)){
    stop("The number of successes (x) must be non-negative.")
  }
  
  if(any(x > n)){
    stop("The number of trials (n) must be greater than the number of successes (x)")
  }
  
  if(any(beta <= 0) || any(alpha <= 0)){
    stop("Both alpha and beta must be postive.")
  }
  
  return(sum(lchoose(n,x) + lbeta(x + alpha,n - x + beta) - lbeta(alpha,beta)))
}
betaBinom_LL<-function(referenceCases,referencePop,simCases){
  
  #apply(simCases,1, dLogBetaBinomial,x=referenceCases,n=referencePop)
  ll_sim=rep(NA,nrow(simCases))
  for(j in 1:nrow(simCases)){
    ll_sim[j]=dLogBetaBinomial(referenceCases,
                            referencePop,
                            as.numeric(simCases[j,]) + 1,
                            referencePop - as.numeric(simCases[j,]) + 1)   
  }

  return(mean(ll_sim))
}

binom_ll<-function(referenceCases,referencePop,simCases){
  
  #apply(simCases,1, dLogBetaBinomial,x=referenceCases,n=referencePop)
  ll_sim=rep(NA,nrow(simCases))
  for(j in 1:nrow(simCases)){
    ll_sim[j]=sum(dbinom(referenceCases,
                               referencePop,
                               as.numeric(simCases[j,]+0.01)/referencePop,log=T))   
  }
  
  return(mean(ll_sim))
}

pois_ll<-function(referenceCases,referencePop,simCases){
  
  #apply(simCases,1, dLogBetaBinomial,x=referenceCases,n=referencePop)
  ll_sim=rep(NA,nrow(simCases))
  for(j in 1:nrow(simCases)){
    ll_sim[j]=sum(dpois(referenceCases,
                         as.numeric(simCases[j,]+0.01),log=T))   
  }
  
  return(mean(ll_sim))
}


## Gather testpos.time from outputs
R0.sweep = read_csv(file.path(output_dir,'params_file.csv'))
testpos.time.q = read_csv(file.path(output_dir,'calibration_sim_quarantine.csv'))
testpos.time.t = read_csv(file.path(output_dir,'calibration_sim_training.csv'))
##Recruit data
wr_dat=read_csv(file="../output/estimated_recruits_testing_positive.csv")
wr_dat=wr_dat[which(!is.na(wr_dat$Battalion_fill)&!is.na(wr_dat$Training_infected)&wr_dat$Battalion_fill>0),]
##==================================================#
## Get R0 intervals----------
##==================================================#
source('paramDefaults.R')
if(hiAsymp){
  symp.mag=0.31
}else{
  symp.mag=0.57  
}
params_df = read.csv('../output/calibrated_parameters.csv',stringsAsFactors = F)

first.col=which(names(testpos.time.q)==gsub(" ",".",wr_dat$Week[1]))
last.col=which(names(testpos.time.q)==gsub(" ",".",wr_dat$Week[nrow(wr_dat)]))

R0.sweep$R0_prob = R0.sweep$R0 * symp.mag + R0.sweep$R0 * (1 - symp.mag)*asymp.adjust
R0.sweep$LL_val=NA


library(foreach)
library(parallel)
library(doParallel)

cl <- makeCluster(10)  
registerDoParallel(cl) 

R0.sweep$LL_val=foreach(ii=1:nrow(R0.sweep),.combine=c) %dopar% {
  #if(ii %% 500==0) print(ii)
  tmp_df_q=testpos.time.q[which(testpos.time.q$uID==R0.sweep$uID[ii]),first.col:last.col]
  tmp_df_t=testpos.time.t[which(testpos.time.t$uID==R0.sweep$uID[ii]),first.col:last.col]
  tot_df=tmp_df_q+tmp_df_t
  betaBinom_LL(wr_dat$Training_infected+wr_dat$Quarantine_positives,
                                     wr_dat$Battalion_fill,
                                     tot_df)
}
#stopCluster(cl)

plot(R0.sweep$R0,R0.sweep$LL_val)
plot(R0.sweep$disp_log,R0.sweep$LL_val)

##GAM
library(mgcv)
library(ggplot2)

R0.sweep.full=R0.sweep
R0.sweep=R0.sweep[R0.sweep$disp_log>0,]
AIC.kk=c()
#cl <- makeCluster(12)  
#registerDoParallel(cl) 
AIC.kk=foreach(kk=seq(100,1000,by=100),.packages = 'mgcv',.export = 'R0.sweep') %dopar% {
  AIC(gam(LL_val ~ s(R0, disp_log, bs= "gp", k = kk), data = R0.sweep)) ## Matern spline
}
AIC.kk=unlist(AIC.kk)
kval=seq(100,1000,by=100)[which.min(AIC.kk)]
gam.rdB  <- gam(LL_val ~ s(R0, disp_log, bs= "gp", k = kval), data = R0.sweep) ## Matern spline
ppB <- predict(gam.rdB, type = "response")
plot(LL_val ~ R0, data=R0.sweep)
lines(ppB ~ R0, data=R0.sweep, col="red")
vis.gam(gam.rdB, theta=30, main = "Matern")
plot(R0.sweep$LL_val,ppB,xlim=c(-1000,-500),ylim=c(-1000,-500))


n <- 10000
R0.pred=expand.grid(R0=seq(2,20,by=0.005),disp_log=seq(0,4,by=0.05))
post.R <- predict(gam.rdB, type = "response", newdata=R0.pred)
post_samples=sample(1:nrow(R0.pred),size=n,prob=exp(post.R),replace=T)
#Xp <- predict(gam.rdB,type="lpmatrix") ## map coefs to fitted curves
Xp <- predict(gam.rdB,type="lpmatrix",newdata=R0.pred) ## map coefs to fitted curves

beta <- coef(gam.rdB)
Vb   <- vcov(gam.rdB) ## posterior mean and cov of coefs
library("MASS") ## for mvrnorm
set.seed(10)
mrand <- mvrnorm(n, beta, Vb)

R_opt = disp_opt= adj_opt=rep(NA, n)
ilink <- family(gam.rdB)$linkinv
pred_list=foreach(i=seq_len(n),combine=cbind) %dopar% { 
  pred   <- ilink(Xp %*% mrand[i, ])
  iR_opt <- R0.pred$R0[which.max(pred)]
  idisp_opt <- R0.pred$disp_log[which.max(pred)]
  #adj_opt[i] <- R0.pred$R0_adj[which.max(pred)]
  return(c(iR_opt,idisp_opt))
}
stopCluster(cl)
library(data.table)
R0_samples=rbindlist(lapply(pred_list,function(x) as.data.frame(matrix(x,ncol=2))))
names(R0_samples)=c("R0","disp")
R0_samples$disp=10^R0_samples$disp
R0_samples$R0_adj=1

if(hiAsymp){
  saveRDS(R0_samples,'../output/R0_adjust_posterior_samples_dispersion_hiAsymp.rds')
}else{
  saveRDS(R0_samples,'../output/R0_adjust_posterior_samples_dispersion_loAsymp.rds') 
}


params_df$R0_mag[params_df$base_name == 'base']  = quantile(R0_samples$R0,c(0.5))
params_df$R0_max[params_df$base_name == 'base']  = quantile(R0_samples$R0,c(.975))
params_df$R0_min[params_df$base_name == 'base']  = quantile(R0_samples$R0,c(0.025))
params_df$dispersion_mag[params_df$base_name == 'base']  = quantile(R0_samples$disp,c(0.5))
params_df$dispersion_max[params_df$base_name == 'base']  = quantile(R0_samples$disp,c(0.975))
params_df$dispersion_min[params_df$base_name == 'base']  = quantile(R0_samples$disp,c(0.025))

params_df$R0_prob_med = params_df$R0_mag * symp.mag + params_df$R0_mag * (1 - symp.mag)*asymp.adjust
params_df$R0_prob_min = params_df$R0_min * symp.mag + params_df$R0_min * (1 - symp.mag)*asymp.adjust
params_df$R0_prob_max = params_df$R0_max * symp.mag + params_df$R0_max * (1 - symp.mag)*asymp.adjust

if(hiAsymp){
  write.csv(params_df, '../output/calibrated_parameters_forecast_hiAsymp.csv', row.names = F)
}else{
  write.csv(params_df, '../output/calibrated_parameters_forecast_loAsymp.csv', row.names = F)
}




