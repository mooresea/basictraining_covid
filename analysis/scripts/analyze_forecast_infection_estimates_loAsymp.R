library(reshape2)
library(tidyverse)
library(rriskDistributions)
library(foreach)
library(doParallel)
library(parallel)
source('simOutbreak_data_betabinomR0_vaccination_updated2.R')

load(file="../output/forecasted_infection_estimates_loasymp.rda")

combo_dat=merge(init.testDay2,inf.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,pos.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,posT.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,imp.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,infQ.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,infT.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,infA.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,infA2.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,sympQ.testDay2,by=c('rep','week','testing'))
combo_dat=merge(combo_dat,sympT.testDay2,by=c('rep','week','testing'))


##Add uncertainty values for data

combo_dat=merge(combo_dat,wr_dat,by.x="week",by.y="week_num")

##summary data
combo_sum = combo_dat %>% filter(testing=="Both") %>% group_by(Start_date) %>%
  summarize(Arrival_num=mean(Arrivals),
            Arrival_pos=mean(Positives),
            Battalion_num=mean(Battalion_fill),
            Quarantine_obs=mean(Quarantine_positives),
            Training_obs=mean(Training_infected),
            training_pos_median=median(training_positives),
            training_pos_lo=quantile(training_positives,0.025),
            training_pos_hi=quantile(training_positives,0.975),
            quarantine_pos_median=median(quarantine_positives),
            quarantine_pos_lo=quantile(quarantine_positives,0.025),
            quarantine_pos_hi=quantile(quarantine_positives,0.975),
            arrival_pos_median=median(arrival_positives),
            arrival_pos_lo=quantile(arrival_positives,0.025),
            arrival_pos_hi=quantile(arrival_positives,0.975))

##summary data
combo_sum_test = combo_dat %>% group_by(Start_date,testing) %>%
  summarize(Arrival_num=mean(Arrivals),
            Arrival_pos=mean(Positives),
            Battalion_num=mean(Battalion_fill),
            Total_size=mean(Total_size),
            Quarantine_obs=mean(Quarantine_positives),
            Training_obs=mean(Training_infected),
            training_pos_median=median(training_positives),
            training_pos_lo=quantile(training_positives,0.025),
            training_pos_hi=quantile(training_positives,0.975),
            quarantine_pos_median=median(quarantine_positives),
            quarantine_pos_lo=quantile(quarantine_positives,0.025),
            quarantine_pos_hi=quantile(quarantine_positives,0.975),
            arrival_pos_median=median(arrival_positives),
            arrival_pos_lo=quantile(arrival_positives,0.025),
            arrival_pos_hi=quantile(arrival_positives,0.975),
            quarantine_inf_median=median(quarantine_infections),
            training_inf_median=median(training_infections),
            quarantine_inf_lo=quantile(quarantine_infections,0.025),
            quarantine_inf_hi=quantile(quarantine_infections,0.975),
            training_inf_lo=quantile(training_infections,0.025),
            training_inf_hi=quantile(training_infections,0.975),
            total_inf_median=median(total_infections),
            total_inf_lo=quantile(total_infections,0.025),
            total_inf_hi=quantile(total_infections,0.975),
            total_sympI_median=median(symptomatic_quarantine_infections+symptomatic_training_infections),
            total_sympI_lo=quantile(symptomatic_quarantine_infections+symptomatic_training_infections,0.025),
            total_sympI_hi=quantile(symptomatic_quarantine_infections+symptomatic_training_infections,0.975),
            total_pos_median=median(quarantine_positives+training_positives),
            total_pos_lo=quantile(quarantine_positives+training_positives,0.025),
            total_pos_hi=quantile(quarantine_positives+training_positives,0.975),
            imported_inf_median=median(imported_infections),
            imported_inf_lo=quantile(imported_infections,0.025),
            imported_inf_hi=quantile(imported_infections,0.975))
combo_sum_test$Staff_size=combo_sum_test$Total_size-combo_sum_test$Battalion_num

##Arrival Stats
#View(combo_sum_test[combo_sum_test$testing=="Both",])
combo_sum_test$arrival_pct_m=combo_sum_test$arrival_pos_median/combo_sum_test$Battalion_num
combo_sum_test$arrival_pct_l=combo_sum_test$arrival_pos_lo/combo_sum_test$Battalion_num
combo_sum_test$arrival_pct_h=combo_sum_test$arrival_pos_hi/combo_sum_test$Battalion_num

##Uncertainty for observations
combo_sum_test$q_pct=combo_sum_test$Quarantine_obs/combo_sum_test$Battalion_num
combo_sum_test$q_lo=qbinom(0.025,combo_sum_test$Battalion_num,combo_sum_test$Quarantine_obs/combo_sum_test$Battalion_num)/combo_sum_test$Battalion_num
combo_sum_test$q_pct_hi=qbinom(0.975,combo_sum_test$Battalion_num,combo_sum_test$Quarantine_obs/combo_sum_test$Battalion_num)/combo_sum_test$Battalion_num
combo_sum_test$t_pct=combo_sum_test$Training_obs/combo_sum_test$Battalion_num
combo_sum_test$t_lo=qbinom(0.025,combo_sum_test$Battalion_num,combo_sum_test$Training_obs/combo_sum_test$Battalion_num)/combo_sum_test$Battalion_num
combo_sum_test$t_pct_hi=qbinom(0.975,combo_sum_test$Battalion_num,combo_sum_test$Training_obs/combo_sum_test$Battalion_num)/combo_sum_test$Battalion_num
combo_sum_test$tot_obs_pct=(combo_sum_test$Quarantine_obs+combo_sum_test$Training_obs)/combo_sum_test$Battalion_num
combo_sum_test$tot_obs_lo=qbinom(0.025,combo_sum_test$Battalion_num,(combo_sum_test$Quarantine_obs+combo_sum_test$Training_obs)/combo_sum_test$Battalion_num)/combo_sum_test$Battalion_num
combo_sum_test$tot_obs_pct_hi=qbinom(0.975,combo_sum_test$Battalion_num,(combo_sum_test$Quarantine_obs+combo_sum_test$Training_obs)/combo_sum_test$Battalion_num)/combo_sum_test$Battalion_num

###Phase colors
col_arrive="Blue"
col_q="Gold"
col_t="Green"
col_both="Purple"
col_imp="grey"

##Data correlations
wr_datD=wr_dat
##Calculate infection pressure during training camp for Richland Cty
cweeks=which(!is.na(wr_datD$Training_infected))
wr_datD$offbase_inf=NA
for(ii in cweeks){
  ii_start=wr_datD$Start_date[ii]
  wr_datD$offbase_inf[ii]=sum(inf_sc$pct_inf[which(inf_sc$date %in% seq(ii_start,ii_start+69,by=1))])
  
}

cor.test(wr_datD$Positives/wr_datD$Arrivals,wr_datD$Quarantine_positives/wr_datD$Battalion_fill,use="pairwise.complete.obs",alternative="greater")
cor.test(wr_datD$Positives/wr_datD$Arrivals,(wr_datD$Quarantine_positives+wr_datD$Training_infected)/wr_datD$Battalion_fill,use="pairwise.complete.obs",alternative="greater")
cor.test(wr_datD$Quarantine_positives/wr_datD$Battalion_fill,(wr_datD$Training_infected)/wr_datD$Battalion_fill,use="pairwise.complete.obs",alternative="greater")

cor.test(wr_datD$Quarantine_positives/wr_datD$Battalion_fill,wr_datD$offbase_inf,use="pairwise.complete.obs",alternative="greater")
cor.test(wr_datD$Training_infected/wr_datD$Battalion_fill,wr_datD$offbase_inf,use="pairwise.complete.obs",alternative="greater")
cor.test((wr_datD$Quarantine_positives+wr_datD$Training_infected)/wr_datD$Battalion_fill,wr_datD$offbase_inf,use="pairwise.complete.obs",alternative="greater")
cor.test(wr_datD$Positives/wr_datD$Arrivals,wr_datD$offbase_inf,use="pairwise.complete.obs",alternative="greater")

theme_set(theme_bw())
vcc=c("Cocoon" = "Gold", "Training" = "Purple")
p_arrival_corr<-ggplot(wr_datD,aes(x=Positives/Arrivals,y=Quarantine_positives/Battalion_fill))+geom_point(color=col_q,cex=2)+
  geom_point(aes(x=Positives/Arrivals,y=(Quarantine_positives+Training_infected)/Battalion_fill),color=col_both,cex=2)+
  xlab("Portion of trainees\npositive on arrival")+ylab("Portion positive\nduring training camp")+xlim(0.0,0.13)

p_arrival_corr_l<-ggplot(wr_datD,aes(x=Positives/Arrivals,y=Quarantine_positives/Battalion_fill,color="Cocoon"))+geom_point(cex=2)+
  geom_point(aes(x=Positives/Arrivals,y=(Quarantine_positives+Training_infected)/Battalion_fill,color="Training"),cex=2)+
  xlab("Portion of trainees\npositive on arrival")+ylab("Portion positive\nduring training camp")+xlim(0.0,0.13)+
  scale_color_manual(name="Training period",values=vcc,labels=c("Cocoon","Entire training camp"))+theme(legend.position = "top") #c(0.2, 0.88))

p_q_corr<-ggplot(wr_datD,aes(x=Quarantine_positives/Battalion_fill,y=Training_infected/Battalion_fill))+geom_point(color="darkgreen",cex=2)+
  theme_bw()+xlab("Portion of trainees\npositive during cocoon period")+ylab("Portion positive\nduring post-cocoon period")

p_imp_corr<-ggplot(wr_datD,aes(x=offbase_inf,y=(Quarantine_positives+Training_infected)/Battalion_fill))+geom_point(color="darkgrey",cex=2)+
  theme_bw()+xlab("Local incidence rate\n(during entire training camp)")+ylab("Portion positive\nduring training camp")

##Legend only
library(grid)
library(gridExtra)
library(cowplot)
legend <- get_legend(p_arrival_corr_l)                    
# Create new plot window
grid.newpage()                              
# Draw Only legend 
grid.draw(legend)


plot_corr<-plot_grid(p_arrival_corr,p_q_corr,p_imp_corr,ncol=3,labels="AUTO")
plot_corr2<-plot_grid(legend,plot_corr,ncol=1,labels=NULL)
ggsave("../figures/loAsymp/observation_correlations.pdf",plot_corr2,width=8,height=5,pointsize=14,units="in")


## Confidence intervals for observations
wr_datA$arrival_pct=wr_datA$Positives/wr_datA$Arrivals
wr_datA$arrival_pct_lo=qbinom(0.025,wr_datA$Arrivals,wr_datA$Positives/wr_datA$Arrivals)/wr_datA$Arrivals
wr_datA$arrival_pct_hi=qbinom(0.975,wr_datA$Arrivals,wr_datA$Positives/wr_datA$Arrivals)/wr_datA$Arrivals

plotT<-ggplot(combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=training_pos_median/Battalion_num))+geom_line(color=col_t,lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Recruits positive\npost-cocoon")+
  geom_ribbon(aes(x=Start_date,ymin=training_pos_lo/Battalion_num,ymax=training_pos_hi/Battalion_num),fill=col_t,alpha=0.25)+
  #geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=(Training_obs)/Battalion_num),col="black",cex=2)+
  geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=t_pct),col="black",cex=2)+
  geom_segment(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,xend=Start_date,y=t_lo,yend=t_pct_hi),col="black")+  
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)

plotT2<-ggplot(combo_sum_test[combo_sum_test$testing=="Both"&combo_sum_test$Start_date<"2021-05-10",],aes(x=Start_date,y=training_pos_median/Battalion_num))+geom_line(color=col_t,lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\n positive during training (post-cocoon)")+xlim(wr_dat$Start_date[1],wr_dat$Start_date[22])+
  geom_ribbon(aes(x=Start_date,ymin=training_pos_lo/Battalion_num,ymax=training_pos_hi/Battalion_num),fill=col_t,alpha=0.25)+
  geom_point(data=combo_sum_test[combo_sum_test$testing=="Both"&combo_sum_test$Start_date<"2021-05-10",],aes(x=Start_date,y=(Training_obs)/Battalion_num),col="black",cex=2)

plotA<-ggplot(combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=arrival_pos_median/Battalion_num))+geom_line(color=col_arrive,lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\n positive on arrival")+
  geom_ribbon(aes(x=Start_date,ymin=arrival_pos_lo/Battalion_num,ymax=arrival_pos_hi/Battalion_num),fill=col_arrive,alpha=0.3)+
  geom_point(data=wr_datA,aes(x=Start_date,y=arrival_pct),col="black",cex=2)+
  geom_segment(data=wr_datA,aes(x=Start_date,xend=Start_date,y=arrival_pct_lo,yend=arrival_pct_hi),col="black")+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)

plotQ<-ggplot(combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=quarantine_pos_median/Battalion_num))+geom_line(color=col_q,lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Recruits positive\nduring cocoon period")+
  geom_ribbon(aes(x=Start_date,ymin=quarantine_pos_lo/Battalion_num,ymax=quarantine_pos_hi/Battalion_num),fill=col_q,alpha=0.25)+
  #geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=(Quarantine_obs)/Battalion_num),col="black",cex=2)+
  geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=q_pct),col="black",cex=2)+
  geom_segment(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,xend=Start_date,y=q_lo,yend=q_pct_hi),col="black")+  
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)

plotQT<-ggplot(combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=(quarantine_pos_median+training_pos_median)/Battalion_num))+geom_line(color=col_both,lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\n positive during training camp")+
  geom_ribbon(aes(x=Start_date,ymin=(quarantine_pos_lo+training_pos_lo)/Battalion_num,ymax=(quarantine_pos_hi+training_pos_hi)/Battalion_num),fill=col_both,alpha=0.25)+
  #geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=(Quarantine_obs+Training_obs)/Battalion_num),col="black",cex=2)+
  geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=tot_obs_pct),col="black",cex=2)+
  geom_segment(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,xend=Start_date,y=tot_obs_lo,yend=tot_obs_pct_hi),col="black")+  
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)


pdf(file="../figures/loAsymp/positives_on_arrival_forecast_ribbon3.pdf",width=8,height=6,pointsize=14)
plotA
dev.off()

pdf(file="../figures/loAsymp/positives_in_quarantine_MLE_forecast_ribbon3.pdf",width=8,height=6,pointsize=14)
plotQ
dev.off()
pdf(file="../figures/loAsymp/positives_during_training_MLE_forecast_ribbon3.pdf",width=8,height=6,pointsize=14)
plotT
dev.off()

pdf(file="../figures/loAsymp/positives_during_training_noforecast.pdf",width=8,height=6,pointsize=14)
plotT2
dev.off()

pdf(file="../figures/loAsymp/positives_during_All_training_MLE_forecast_ribbon3.pdf",width=8,height=6,pointsize=14)
plotQT
dev.off()

##Combined figure
#library(gridExtra)
#plotC<-grid.arrange(plotA,plotQ,plotT,ncol=1)
library(cowplot)
plotC<-plot_grid(plotA,plotQ,plotT,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/combo_infections_plot_ribbon3.pdf",plotC,width=8,height=14,pointsize=14,units="in")

##Also look at imports
plotImp<-ggplot(combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=imported_inf_median/Staff_size))+geom_line(color="black",lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Trainers and staff\ninfected off base")+
  geom_ribbon(aes(x=Start_date,ymin=imported_inf_lo/Staff_size,ymax=imported_inf_hi/Staff_size),fill="grey",alpha=0.5)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)

plotC4<-plot_grid(plotA,plotQ,plotT,plotImp,ncol=2,labels="AUTO")

ggsave("../figures/loAsymp/combo_infections_imports_plot_ribbon3.pdf",plotC4,width=14,height=10,pointsize=14,units="in")

plotC2<-plot_grid(plotQ,plotT,plotImp,ncol=1,labels=c("C","D","E"))
plotC1<-plot_grid(plotA,plotQT,ncol=1,labels="AUTO")
plotC1_alt<-plot_grid(plotA,plotQT,plotImp,ncol=1,labels="AUTO")
plotCC<-plot_grid(plotC1,plotC2,ncol=2)
ggsave("../figures/loAsymp/combo_positives_imports_plot_2panel.pdf",plotCC,width=8.5,height=8.5/1.4,pointsize=14,units="in")
ggsave("../figures/loAsymp/combo_arrivals_camp.pdf",plotC1,width=8,height=10,pointsize=14,units="in")
ggsave("../figures/loAsymp/combo_arrivals_camp_imports.pdf",plotC1_alt,width=8,height=15,pointsize=14,units="in")


combo_a=combo_sum_test %>% filter(testing=="Both") %>%
  select(Start_date=Start_date,Arrival_pos=Arrival_pos,Arrival_num=Arrival_num,
         Battalion_num=Battalion_num,Staff_size=Staff_size,
         arrival_pos_median=arrival_pos_median,arrival_pos_lo=arrival_pos_lo,arrival_pos_hi=arrival_pos_hi)
combo_a$value="Trainees"
combo_a$inf_median=combo_a$arrival_pos_median/combo_a$Battalion_num
combo_a$inf_lo=combo_a$arrival_pos_lo/combo_a$Battalion_num
combo_a$inf_hi=combo_a$arrival_pos_hi/combo_a$Battalion_num
combo_i=combo_sum_test %>% filter(testing=="Both") %>%
  select(Start_date,Arrival_pos,Arrival_num,
         Battalion_num,Staff_size,
         imported_inf_median,imported_inf_lo,imported_inf_hi)
combo_i$value="Staff"
combo_i$inf_median=combo_i$imported_inf_median/combo_i$Staff_size
combo_i$inf_lo=combo_i$imported_inf_lo/combo_i$Staff_size
combo_i$inf_hi=combo_i$imported_inf_hi/combo_i$Staff_size
combo_a=combo_a %>% select(Start_date,Arrival_pos,Arrival_num,value,
                           inf_median,inf_lo,inf_hi)
combo_i=combo_i %>% select(Start_date,Arrival_pos,Arrival_num,value,
                           inf_median,inf_lo,inf_hi)
combo_ai=bind_rows(combo_a,combo_i)

plotA2<-ggplot(combo_ai,aes(x=Start_date,y=inf_median,color=as.factor(value)))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\n positive on arrival")+
  geom_ribbon(aes(x=Start_date,ymin=inf_lo,ymax=inf_hi,fill=as.factor(value)),alpha=0.2,color=NA)+
  #geom_ribbon(aes(x=Start_date,ymin=arrival_pos_lo/Battalion_num,ymax=arrival_pos_hi/Battalion_num),fill="blue",alpha=0.3)+
  geom_point(data=combo_sum_test[combo_sum_test$testing=="Both",],aes(x=Start_date,y=(Arrival_pos)/Arrival_num),col="black",cex=2)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  scale_color_discrete(name="Arrival\nTesting",breaks=c("Quarantine","Both"),labels=c("No","Yes"))+
  scale_fill_discrete(name="Arrival\nTesting",breaks=c("Quarantine","Both"),labels=c("No","Yes"))

plotCAlt<-plot_grid(plotA2,plotQ,plotT,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/combo_infections_imports_plot_ribbon3_alt.pdf",plotCAlt,width=8,height=14,pointsize=14,units="in")

###
### Intervention Impacts
###
testing_colors=c("red2", "cyan3", "olivedrab3","orchid3")
##Quarantine
pdf(file="../figures/loAsymp/arrival_testing_impact_quarantine.pdf",width=8,height=6,pointsize=14)
plot_at<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Quarantine","Both"),],aes(x=Start_date,y=quarantine_pos_median/Battalion_num,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\ntesting positive\nduring cocoon")+ylim(0,.6)+
  geom_ribbon(aes(x=Start_date,ymin=quarantine_pos_lo/Battalion_num,ymax=quarantine_pos_hi/Battalion_num,fill=testing),alpha=0.2,color=NA)+
  geom_point(data=combo_sum,aes(x=Start_date,y=(Quarantine_obs)/Battalion_num),col="black",cex=2)+
  geom_vline(xintercept = as.Date("2022-01-24"),lty=2)+scale_color_manual(name="Arrival\nTesting",breaks=c("Quarantine","Both"),labels=c("No","Yes"),values=testing_colors[1:2])+
  scale_fill_manual(name="Arrival\nTesting",breaks=c("Quarantine","Both"),labels=c("No","Yes"),values=testing_colors[1:2])
plot_at
dev.off()

# ggplot(combo_sum_test[combo_sum_test$testing %in% c("Quarantine","Both","Arrival","None"),],aes(x=Start_date,y=quarantine_inf_median/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
#   xlab("Arrival date")+ylab("Portion of recruits infected during quarantine")+
#   geom_ribbon(aes(x=Start_date,ymin=quarantine_inf_lo/Total_size,ymax=quarantine_inf_hi/Total_size,fill=testing),alpha=0.25)+
#   geom_vline(xintercept = as.Date("2022-03-02"),lty=2)

plot_att<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Quarantine","Both"),],aes(x=Start_date,y=training_pos_median/Battalion_num,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\ntesting positive\npost cocoon")+ylim(0,.6)+
  geom_ribbon(aes(x=Start_date,ymin=training_pos_lo/Battalion_num,ymax=training_pos_hi/Battalion_num,fill=testing),alpha=0.2,color=NA)+
  geom_point(data=combo_sum,aes(x=Start_date,y=(Training_obs)/Battalion_num),col="black",cex=2)+
  geom_vline(xintercept = as.Date("2022-01-24"),lty=2)+scale_color_manual(name="Arrival\nTesting",breaks=c("Quarantine","Both"),labels=c("No","Yes"),values=testing_colors[1:2])+
  scale_fill_manual(name="Arrival\nTesting",breaks=c("Quarantine","Both"),labels=c("No","Yes"),values=testing_colors[1:2])


##Training
pdf(file="../figures/loAsymp/quarantine_testing_impact_training.pdf",width=8,height=6,pointsize=14)
plot_qt<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Arrival"),],aes(x=Start_date,y=training_pos_median/Battalion_num,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of recruits\ntesting positive\npost cocoon")+ylim(0,0.6)+
  geom_ribbon(aes(x=Start_date,ymin=training_pos_lo/Battalion_num,ymax=training_pos_hi/Battalion_num,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_point(data=combo_sum,aes(x=Start_date,y=(Training_obs)/Battalion_num),col="black",cex=2)+
  scale_fill_manual(name="Cocoon\nTesting",breaks=c("Arrival","Both"),labels=c("No","Yes"),values=testing_colors[c(1,3)])+
  scale_color_manual(name="Cocoon\nTesting",breaks=c("Arrival","Both"),labels=c("No","Yes"),values=testing_colors[c(1,3)])
plot_qt
dev.off()

##Combined surveillance testing figure
plot_sti<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Quarantine","Both","Arrival","None"),],aes(x=Start_date,y=(total_inf_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion infected\nduring training camp")+ylim(0,1)+
  geom_ribbon(aes(x=Start_date,ymin=(total_inf_lo)/Total_size,ymax=(total_inf_hi)/Total_size,fill=testing),alpha=0.1,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  scale_fill_manual(name="Surveillance\nTesting",breaks=c("Arrival","Quarantine","Both","None"),c("Arrival","Cocoon","Both","None"),values=testing_colors[c(2,3,4,1)])+
  scale_color_manual(name="Surveillance\nTesting",breaks=c("Arrival","Quarantine","Both","None"),labels=c("Arrival","Cocoon","Both","None"),values=testing_colors[c(2,3,4,1)])

plot_stp<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Quarantine","Both","Arrival","None"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion with symptomatic\n infection during training camp")+ylim(0,0.4)+
  geom_ribbon(aes(x=Start_date,ymin=(total_sympI_lo)/Total_size,ymax=(total_sympI_hi)/Total_size,fill=testing),alpha=0.1,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  scale_fill_manual(name="Surveillance\nTesting",breaks=c("Arrival","Quarantine","Both","None"),c("Arrival","Cocoon","Both","None"),values=testing_colors[c(2,3,4,1)])+
  scale_color_manual(name="Surveillance\nTesting",breaks=c("Arrival","Quarantine","Both","None"),labels=c("Arrival","Cocoon","Both","None"),values=testing_colors[c(2,3,4,1)])
ggsave("../figures/loAsymp/surveillance_testing_symp.pdf",plot_stp,width=6,height=4,pointsize=14,units="in")

plotC_st<-plot_grid(plot_at,plot_sti,plot_att,plot_qt,ncol=2,labels="AUTO")

ggsave("../figures/loAsymp/combo_surveillance_testing_impact.pdf",plotC_st,width=14,height=10,pointsize=14,units="in")

plotC_st_c1<-plot_grid(plot_at,plot_att,plot_qt,ncol=1,labels="AUTO")
plotC_st_c2<-plot_grid(plot_sti,plot_stp,ncol=1,labels=c("D","E"))
plotC_st_b<-plot_grid(plotC_st_c1,plotC_st_c2,ncol=2)

ggsave("../figures/loAsymp/combo_surveillance_testing_impact_2panel.pdf",plotC_st_b,width=8.5,height=8/1.4,pointsize=14,units="in")


##
## Masks
##
mask_colors=c("red2", "green3", "greenyellow","orchid3")
mask_colors=c("red2", "green3", "cyan2","orchid3")
pdf(file="../figures/loAsymp/masks_impact_infections.pdf",width=8,height=6,pointsize=14)
plot_mi<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","No Masks","High Masks"),],aes(x=Start_date,y=(total_inf_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion infected\n during training camp")+ylim(0,1)+
  geom_ribbon(aes(x=Start_date,ymin=(total_inf_lo)/Total_size,ymax=(total_inf_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  scale_fill_manual(name="Masking",breaks=c("No Masks","Both","High Masks"),labels=c("None","Baseline","100%"),values=mask_colors[1:3])+
  scale_color_manual(name="Masking",breaks=c("No Masks","Both","High Masks"),labels=c("None","Baseline","100%"),values=mask_colors[1:3])
plot_mi
dev.off()
pdf(file="../figures/loAsymp/masks_impact_symp_infections.pdf",width=8,height=6,pointsize=14)
plot_ms<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","No Masks","High Masks"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion with\nsymptomatic infection")+ylim(0,0.4)+
  geom_ribbon(aes(x=Start_date,ymin=(total_sympI_lo)/Total_size,ymax=(total_sympI_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  scale_fill_manual(name="Masking",breaks=c("No Masks","Both","High Masks"),labels=c("None","Baseline","100%"),values=mask_colors[1:3])+
  scale_color_manual(name="Masking",breaks=c("No Masks","Both","High Masks"),labels=c("None","Baseline","100%"),values=mask_colors[1:3])
plot_ms
dev.off()

pdf(file="../figures/loAsymp/intervention_impact_infections.pdf",width=8,height=6,pointsize=14)
plot_ii<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","No Masks","None"),],aes(x=Start_date,y=(total_inf_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion infected\n during training camp")+ylim(0,1)+
  geom_ribbon(aes(x=Start_date,ymin=(total_inf_lo)/Total_size,ymax=(total_inf_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  scale_fill_manual(name="Intervention",breaks=c("No Masks","None","Both"),labels=c("Testing","Masks","Both"),values=c(testing_colors[2],mask_colors[c(2,4)]))+
  scale_color_manual(name="Intervention",breaks=c("No Masks","None","Both"),labels=c("Testing","Masks","Both"),values=c(testing_colors[2],mask_colors[c(2,4)]))
plot_ii
dev.off()

##Vaccine

pdf(file="../figures/loAsymp/vaccine_impact_infections.pdf",width=8,height=6,pointsize=14)
plot_vi<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","No Vaccine","High Vaccine"),],aes(x=Start_date,y=(total_inf_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion infected\n during training camp")+ylim(0,1)+
  geom_ribbon(aes(x=Start_date,ymin=(total_inf_lo)/Total_size,ymax=(total_inf_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_vline(xintercept = as.Date("2021-05-22"),lty=3,col="black")+
  scale_fill_manual(name="Arrival\nVaccination",breaks=c("No Vaccine","Both","High Vaccine"),labels=c("0%","Baseline","100%"),values=mask_colors[1:3])+
  scale_color_manual(name="Arrival\nVaccination",breaks=c("No Vaccine","Both","High Vaccine"),labels=c("0%","Baseline","100%"),values=mask_colors[1:3])
plot_vi
dev.off()


pdf(file="../figures/loAsymp/vaccine_impact_symptomatic_infections.pdf",width=8,height=6,pointsize=14)
plot_vs<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","No Vaccine","High Vaccine"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion with\nsymptomatic infection")+ylim(0,0.4)+
  geom_ribbon(aes(x=Start_date,ymin=(total_sympI_lo)/Total_size,ymax=(total_sympI_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_vline(xintercept = as.Date("2021-05-22"),lty=3,col="black")+
  scale_fill_manual(name="Arrival\nVaccination",breaks=c("No Vaccine","Both","High Vaccine"),labels=c("0%","Baseline","100%"),values=mask_colors[1:3])+
  scale_color_manual(name="Arrival\nVaccination",breaks=c("No Vaccine","Both","High Vaccine"),labels=c("0%","Baseline","100%"),values=mask_colors[1:3])
plot_vs
dev.off()

#masks and early vaccination impact
plotC_vm<-plot_grid(plot_mi,plot_ms,plot_vi,plot_vs,ncol=2,labels="AUTO")

ggsave("../figures/loAsymp/combo_mask_vaccine_impact_symptomatic.pdf",plotC_vm,width=8,height=5,pointsize=14,units="in")

##
## Vaccinating Staff
##
pdf(file="../figures/loAsymp/vaccine_staff_impact_infections.pdf",width=8,height=6,pointsize=14)
plot_vSi<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Staff Vaccine"),],aes(x=Start_date,y=(total_inf_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals\n infected during training camp")+ylim(0,1)+
  geom_ribbon(aes(x=Start_date,ymin=(total_inf_lo)/Total_size,ymax=(total_inf_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_vline(xintercept = as.Date("2021-05-22"),lty=3,col="black")+
  geom_vline(xintercept = as.Date("2021-12-01"),lty=3,col="black")+
  scale_fill_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Vaccine"),labels=c("Voluntary","Required"))+
  scale_color_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Vaccine"),labels=c("Voluntary","Required"))
plot_vSi
dev.off()

pdf(file="../figures/loAsymp/vaccine_staff_impact_symptomatic_infections.pdf",width=8,height=6,pointsize=14)
plot_vSs<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Staff Vaccine"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals with symptomatic\n infections during training camp")+ylim(0,0.4)+
  geom_ribbon(aes(x=Start_date,ymin=(total_sympI_lo)/Total_size,ymax=(total_sympI_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_vline(xintercept = as.Date("2021-05-22"),lty=3,col="black")+
  geom_vline(xintercept = as.Date("2021-12-01"),lty=3,col="black")+
  scale_fill_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Vaccine"),labels=c("Voluntary","Required"))+
  scale_color_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Vaccine"),labels=c("Voluntary","Required"))
plot_vSs
dev.off()

pdf(file="../figures/loAsymp/vaccine_staff_impact_importations.pdf",width=8,height=6,pointsize=14)
plot_vS_imp<-ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Staff Vaccine"),],aes(x=Start_date,y=(imported_inf_median)/Staff_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of trainers and staff\n infected off base")+ylim(0,.2)+
  geom_ribbon(aes(x=Start_date,ymin=(imported_inf_lo)/Staff_size,ymax=(imported_inf_hi)/Staff_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_vline(xintercept = as.Date("2021-05-22"),lty=3,col="black")+
  geom_vline(xintercept = as.Date("2021-12-01"),lty=3,col="black")+
  scale_fill_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Vaccine"),labels=c("Voluntary","Required"))+
  scale_color_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Vaccine"),labels=c("Voluntary","Required"))
plot_vS_imp
dev.off()

#Look at early vaccine
ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Early Vaccine","High Vaccine"),],aes(x=Start_date,y=(total_inf_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals\n infected during training camp")+ylim(0,1)+
  geom_ribbon(aes(x=Start_date,ymin=(total_inf_lo)/Total_size,ymax=(total_inf_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)+
  geom_vline(xintercept = as.Date("2021-05-22"),lty=3,col="black")+
  scale_fill_discrete(name="Arrival\nVaccination",breaks=c("Early Vaccine","Both","High Vaccine"),labels=c("0%","Baseline","100%"))+
  scale_color_discrete(name="Arrival\nVaccination",breaks=c("Early Vaccine","Both","High Vaccine"),labels=c("0%","Baseline","100%"))

##
## Early vaccine more effective than testing at reducing symptomatic infeciotn, but not always for total infections
##
ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Early Vaccine","High Vaccine"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals\n infected during training camp")+ylim(0,0.3)
ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Early Vaccine","Early Vaccine Only"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals\n infected during training camp")+ylim(0,0.3)
ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Early Vaccine","Staff Vaccine"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals\n infected during training camp")+ylim(0,0.3)

#plotC_vm<-plot_grid(plot_mi,plot_vi,ncol=2,labels="AUTO")
plotC_vm2<-plot_grid(plot_mi,plot_vi,plot_vS_imp,plot_vSi,ncol=2,labels="AUTO")

#ggsave("../figures/loAsymp/combo_mask_vaccine_impact.pdf",plotC_vm,width=14,height=6,pointsize=14,units="in")
ggsave("../figures/loAsymp/combo_mask_vaccine_impact_staff.pdf",plotC_vm2,width=14,height=10,pointsize=14,units="in")

##
## Staff surveillance
##
ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Staff Tests","Staff Tests Daily Antigen","Staff Tests Weekly Antigen"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals with symptomatic\n infections during training camp")+ylim(0,0.4)+
  geom_ribbon(aes(x=Start_date,ymin=(total_sympI_lo)/Total_size,ymax=(total_sympI_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)
#scale_fill_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Tests"),labels=c("Voluntary","Tests"))+
#scale_color_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Tests"),labels=c("Voluntary","Tests"))
ggplot(combo_sum_test[combo_sum_test$testing %in% c("Both","Staff Tests Only","Staff Tests Daily Antigen"),],aes(x=Start_date,y=(total_sympI_median)/Total_size,color=testing))+geom_line(lwd=2)+theme_bw()+
  xlab("Arrival date")+ylab("Portion of individuals with symptomatic\n infections during training camp")+ylim(0,0.5)+
  geom_ribbon(aes(x=Start_date,ymin=(total_sympI_lo)/Total_size,ymax=(total_sympI_hi)/Total_size,fill=testing),alpha=0.2,color=NA)+
  geom_vline(xintercept = as.Date("2022-03-02"),lty=2)
#scale_fill_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Tests"),labels=c("Voluntary","Tests"))+
#scale_color_discrete(name="Staff\nVaccination",breaks=c("Both","Staff Tests"),labels=c("Voluntary","Tests"))

combo_dat$total_inf_pct=combo_dat$total_infections/combo_dat$Total_size
combo_dat$total_symp_inf_pct=(combo_dat$symptomatic_quarantine_infections+combo_dat$symptomatic_training_infections)/combo_dat$Total_size
combo_dat$import_inf_pct=combo_dat$imported_infections/(combo_dat$Total_size-combo_dat$Battalion_fill)
combo_dat_test=combo_dat[combo_dat$testing %in% c("Both","None"),]

###Boxplots
plot_pct_q<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(pct_inf_prior_2week*2,2)/2),y=quarantine_infections/Total_size))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(fill=col_q)+theme_bw()+
  scale_x_discrete(breaks=seq(0,0.08,.01))+ylim(0,0.5)+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+ylab("Portion of trainees\n infected during cocoon") #+geom_smooth(method="lm")
plot_pct_t<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(pct_inf_prior_2week*2,2)/2),y=total_inf_pct))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(fill=col_arrive,alpha=0.6)+theme_bw()+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.08,.01))+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")
plot_imp_t<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(import_inf_pct,2)),y=total_inf_pct))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(fill=col_imp,alpha=0.6)+theme_bw()+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.2,.02))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")

plot_pct_imp<-plot_grid(plot_pct_q,plot_pct_t,plot_imp_t,ncol=1,labels="AUTO")

##Violin plots
plot_pct_q_viol<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(pct_inf_prior_2week*2,2)/2),y=quarantine_infections/Total_size))+
  geom_jitter(alpha=0.04,height=0,width=0.1)+#+theme_bw()+geom_smooth()
  geom_violin(fill=col_q,alpha=0.6,draw_quantiles = c(0.25, 0.5, 0.75))+theme_bw()+
  scale_x_discrete(breaks=seq(0,0.08,.01))+ylim(0,0.5)+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+ylab("Portion of trainees\n infected during cocoon") #+geom_smooth(method="lm")
plot_pct_t_viol<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(pct_inf_prior_2week*2,2)/2),y=total_inf_pct))+
  geom_jitter(alpha=0.02,height=0,width=0.1)+#theme_bw()+geom_smooth()
  geom_violin(fill=col_arrive,alpha=0.6,draw_quantiles = c(0.25, 0.5, 0.75))+theme_bw()+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.08,.01))+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")
plot_imp_t_viol<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(import_inf_pct,2)),y=total_inf_pct))+
  geom_jitter(alpha=0.05,height=0,width=0.1)+#theme_bw()+geom_smooth()
  geom_violin(fill=col_imp,alpha=0.6,draw_quantiles = c(0.25, 0.5, 0.75))+theme_bw()+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.2,.02))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")

plot_pct_imp_viol<-plot_grid(plot_pct_q_viol,plot_pct_t_viol,plot_imp_t_viol,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/infections_by_arrival_import_violin.pdf",plot_pct_imp_viol,width=4.5,height=7,pointsize=14,units="in")

plot_pct_imp_poster<-plot_grid(plot_pct_t,plot_imp_t,ncol=2,labels="AUTO")

ggsave("../figures/loAsymp/infections_by_arrival_import_poster.pdf",plot_pct_imp_poster,width=12,height=6,pointsize=14,units="in")

###Infected when no arrival testing
combo_dat$date_ind=ifelse(combo_dat$Start_date<as.Date("2021-07-01"),0,1)
#
#Compare dates
#
combo_dat$period=ifelse(combo_dat$Start_date<"2021-07-01",0,1)
plot_pct_date<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(pct_inf_prior_2week*2,2)/2),y=total_inf_pct,fill=as.factor(period)))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot()+theme_bw()+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.08,.01))+
  scale_fill_discrete(name="Period",breaks=c(0,1),labels=c("Before 7/1/21","After 7/1/21"))+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")
plot_imp_date<-ggplot(combo_dat[combo_dat$testing=="Both",],aes(x=as.factor(round(import_inf_pct,2)),y=total_inf_pct,fill=as.factor(period)))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot()+theme_bw()+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.2,.02))+
  scale_fill_discrete(name="Period",breaks=c(0,1),labels=c("Before 7/1/21","After 7/1/21"))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")
plot_pct_imp_date<-plot_grid(plot_pct_date,plot_imp_date,ncol=2,labels="AUTO")

ggsave("../figures/loAsymp/infections_by_arrival_import_waves.pdf",plot_pct_imp_date,width=14,height=6,pointsize=14,units="in")

##
## Impact of R0 on outbreak size
##
R0_val=median(R0_samples$R0)
R0_dat=var_dat
R0_dat$R0=R0_dat$other_pct*R0_val+
  R0_dat$alpha_pct*R0_val*var_p$trans_mult[var_p$variant=="alpha"]+
  R0_dat$gamma_pct*R0_val*var_p$trans_mult[var_p$variant=="gamma"]+
  R0_dat$delta_pct*R0_val*var_p$trans_mult[var_p$variant=="delta"]+
  R0_dat$omicron_pct*R0_val*var_p$trans_mult[var_p$variant=="omicron"]
#Adjust variant week dates so they match model weeks
tdates=unique(combo_dat$Start_date)
R0_dat$day=as.numeric(R0_dat$week-R0_dat$week[1])

R0.gam=gam(R0~s(day,k=nrow(R0_dat)),data=R0_dat)
R0.day=predict(R0.gam,newdata=data.frame(day=0:max(R0_dat$day)))
R0.df=data.frame(day=0:max(R0_dat$day),R0=R0.day)
R0.df$date=min(R0_dat$week)+R0.df$day
##Add weekly R0 value
combo_dat = combo_dat %>% left_join(R0.df,by=c("Start_date"="date"))

plot_Inf_R0<-ggplot(combo_dat[combo_dat$testing %in% c("Both"),],
                    aes(x=round(R0,0),y=total_inf_pct,group=round(R0,0)))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(outlier.alpha = 0.2,fill="lightblue")+theme(legend.position = "top")+ylim(0,1)+  
  scale_x_continuous(breaks=seq(5,20,by=5),labels=seq(5,20,by=5))+xlim(5,20)+
  xlab(expression(paste("Basic reproduction number (",R[0],")")))+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")

plot_Inf_Imp_R0<-ggplot(combo_dat[combo_dat$testing %in% c("Both"),],
                        aes(x=as.factor(round(import_inf_pct/2,2)*2),y=total_inf_pct,fill=as.factor(round(R0/5,0)*5)))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(outlier.alpha = 0.2)+ylim(0,1)+
  scale_x_discrete(breaks=seq(0,0.21,.04))+
  theme(legend.position = c(0.86,0.27),legend.title = element_text(size=7),legend.text = element_text(size=7))+  
  scale_fill_discrete(name=expression(R[0]),breaks=c(5,10,15,20),labels=c("3-7.9","8-12.9","13-17.9","18-22.9"))+
  #scale_fill_discrete(name=expression(R[0]),breaks=c(12,18,24,30),labels=c("9-14.9","15-20.9","21-26.9","27-32.9"))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")

plot_R0<-plot_grid(plot_Inf_R0,plot_Inf_Imp_R0,ncol=2,labels="AUTO")

ggsave("../figures/loAsymp/infections_by_R0.pdf",plot_R0,width=8,height=3,pointsize=14,units="in")


plot_R0_vert<-plot_grid(plot_Inf_R0,plot_Inf_Imp_R0,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/infections_by_R0_vertical.pdf",plot_R0_vert,width=4.5,height=6,pointsize=14,units="in")

##
## Statistical test of intervention impact
##
impact_dat<-combo_dat %>% group_by(Start_date,testing) %>% summarize(total_infections=mean(total_infections,na.rm = T),
                                                                     total_symp_infections=mean(symptomatic_quarantine_infections+
                                                                                                  symptomatic_training_infections,na.rm=T),
                                                                     pct_inf_prior_2week=mean(pct_inf_prior_2week),
                                                                     import_inf_pct=mean(import_inf_pct))
impact_dat_pv=impact_dat[impact_dat$Start_date>"2021-05-31",]

##
##T-tests and Mann-Whitney tests for impact of different interventions
##

##Arrival and cocoon surveillance testing
at.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Quarantine"],
                 impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95) #,alternative="l")
at.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Quarantine"])
at.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Quarantine"])
at.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Quarantine"])

at.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Quarantine"],
                  impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95)
at.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Quarantine"])
at.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Quarantine"])
at.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Quarantine"])

ct.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Arrival"],
                 impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95) #,alternative="l")
ct.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Arrival"])
ct.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Arrival"])
ct.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Arrival"])

ct.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Arrival"],
                  impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95)
ct.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Arrival"])
ct.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Arrival"])
ct.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Arrival"])

##Antigen testing
dat.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
             impact_dat_pv$total_infections[impact_dat_pv$testing=="Staff Tests Daily Antigen"],paired=T)
dat.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
                  impact_dat_pv$total_infections[impact_dat_pv$testing=="Staff Tests Daily Antigen"],paired=T,conf.int = T,conf.level=0.95)
dat.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
dat.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
dat.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
dat.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
dat.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
dat.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])

wat.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
             impact_dat_pv$total_infections[impact_dat_pv$testing=="Staff Tests Weekly Antigen"],paired=T) #,alternative="l")
wat.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
wat.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
wat.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
wat.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
                  impact_dat_pv$total_infections[impact_dat_pv$testing=="Staff Tests Weekly Antigen"],paired=T,conf.int = T,conf.level=0.95) #,alternative="l")
wat.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
wat.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
wat.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])

dat.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
                   impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Staff Tests Daily Antigen"],paired=T,conf.int = T,conf.level=0.95)
dat.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
dat.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
dat.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
wat.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
                   impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Staff Tests Weekly Antigen"],paired=T,conf.int = T,conf.level=0.95) #,alternative="l")
wat.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
wat.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
wat.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

###Arrival vaccination

v.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Vaccine"],
                impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95)
v.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Vaccine"])
v.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Vaccine"])
v.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Vaccine"])
v.wt=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Vaccine"],
                 impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95)
v.wt$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Vaccine"])
v.wt$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Vaccine"])
v.wt$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Vaccine"])


hv.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
                 impact_dat_pv$total_infections[impact_dat_pv$testing=="High Vaccine"],paired=T,conf.int = T,conf.level=0.95)
hv.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
hv.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
hv.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
hv.wt=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
                  impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="High Vaccine"],paired=T,conf.int = T,conf.level=0.95)
hv.wt$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
hv.wt$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
hv.wt$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

ev.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
                 impact_dat_pv$total_infections[impact_dat_pv$testing=="Early Vaccine"],paired=T,conf.int = T,conf.level=0.95)
ev.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
ev.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
ev.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
ev.wt=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
                  impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Early Vaccine"],paired=T,conf.int = T,conf.level=0.95)
ev.wt$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
ev.wt$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
ev.wt$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

ev.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
            impact_dat_pv$total_infections[impact_dat_pv$testing=="Early Vaccine"],paired=T)
ev.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
ev.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
ev.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
ev.st=t.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
             impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Early Vaccine"],paired=T)
ev.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
ev.st$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
ev.st$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

##Baseline mask
mask.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"],
                    impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95)
mask.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"])
mask.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"])
mask.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"])
mask.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"],
                   impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],paired=T,conf.int = T,conf.level=0.95)
mask.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"])
mask.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"])
mask.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"])

##High mask
hmask.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"],
                    impact_dat_pv$total_infections[impact_dat_pv$testing=="High Masks"],paired=T,conf.int = T,conf.level=0.95)
hmask.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"])
hmask.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"])
hmask.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks"])
hmask.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"],
                     impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="High Masks"],paired=T,conf.int = T,conf.level=0.95)
hmask.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"])
hmask.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"])
hmask.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks"])

sv.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
            impact_dat_pv$total_infections[impact_dat_pv$testing=="Staff Vaccine"],paired=T)
sv.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
sv.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
sv.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
sv.st=t.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
             impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Staff Vaccine"],paired=T)
sv.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
sv.st$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
sv.st$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
sv.w=wilcox.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
                 impact_dat_pv$total_infections[impact_dat_pv$testing=="Staff Vaccine"],paired=T,conf.int = T,conf.level=0.95)
sv.w$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
sv.w$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
sv.w$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
sv.sw=wilcox.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
                  impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Staff Vaccine"],paired=T,conf.int = T,conf.level=0.95)
sv.sw$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
sv.sw$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
sv.sw$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

evst.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
              impact_dat_pv$total_infections[impact_dat_pv$testing=="Early Vaccine and Staff Tests Daily Antigen"],paired=T)
evst.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
evst.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
evst.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
evst.st=t.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
               impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Early Vaccine and Staff Tests Daily Antigen"],paired=T)
evst.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
evst.st$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
evst.st$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
ev.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])+dat.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

all.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"],
             impact_dat_pv$total_infections[impact_dat_pv$testing=="Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine"],paired=T)
all.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
all.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
all.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="Both"])
all.st=t.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"],
              impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine"],paired=T)
all.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
all.st$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])
all.st$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Both"])

allP.t=t.test(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks or Tests"],
              impact_dat_pv$total_infections[impact_dat_pv$testing=="Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine"],paired=T)
allP.t$estimate/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks or Tests"])
allP.t$conf.int[1]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks or Tests"])
allP.t$conf.int[2]/mean(impact_dat_pv$total_infections[impact_dat_pv$testing=="No Masks or Tests"])
allP.st=t.test(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks or Tests"],
               impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine"],paired=T)
allP.st$estimate/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks or Tests"])
allP.st$conf.int[1]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks or Tests"])
allP.st$conf.int[2]/mean(impact_dat_pv$total_symp_infections[impact_dat_pv$testing=="No Masks or Tests"])

##Generate summaries of percent reductions
impact_TS = impact_dat[,c(1:2,4:5)]
impact_TS_w = impact_TS %>% pivot_wider(names_from=testing,values_from=total_symp_infections)
impact_TI = impact_dat[,c(1:3,5)]
impact_TI_w = impact_TI %>% pivot_wider(names_from=testing,values_from=total_infections)
impact_TI_w$diff_test=impact_TI_w$None-impact_TI_w$Both
impact_TI_w$diff_mask=impact_TI_w$'No Masks'-impact_TI_w$'High Masks'
impact_TI_w$diff_test_pct=(impact_TI_w$None-impact_TI_w$Both)/impact_TI_w$None
impact_TI_w$diff_mask_pct=(impact_TI_w$'No Masks'-impact_TI_w$'High Masks')/impact_TI_w$'No Masks'
impact_TI_w$diff_staff_vacc=impact_TI_w$'Both'-impact_TI_w$'Staff Vaccine'
impact_TI_w$diff_staff_vacc_pct=(impact_TI_w$'Both'-impact_TI_w$'Staff Vaccine')/impact_TI_w$'Both'
impact_TI_w$diff_vacc=impact_TI_w$'No Vaccine'-impact_TI_w$'High Vaccine'
impact_TI_w$diff_vacc_pct=(impact_TI_w$'No Vaccine'-impact_TI_w$'High Vaccine')/impact_TI_w$'No Vaccine'
impact_TI_w$diff_early_vacc=impact_TI_w$'Both'-impact_TI_w$'Early Vaccine'
impact_TI_w$diff_early_vacc_pct=(impact_TI_w$'Both'-impact_TI_w$'Early Vaccine')/impact_TI_w$'Both'
impact_TI_w$diff_staff_daily_test=impact_TI_w$'Both'-impact_TI_w$'Staff Tests Daily Antigen'
impact_TI_w$diff_staff_daily_test_pct=(impact_TI_w$'Both'-impact_TI_w$'Staff Tests Daily Antigen')/impact_TI_w$'Both'
impact_TI_w$diff_staff_weekly_test=impact_TI_w$'Both'-impact_TI_w$'Staff Tests'
impact_TI_w$diff_staff_weekly_test_pct=(impact_TI_w$'Both'-impact_TI_w$'Staff Tests')/impact_TI_w$'Both'
impact_TI_w$diff_EV_ST=impact_TI_w$'Both'-impact_TI_w$'Early Vaccine and Staff Tests Daily Antigen'
impact_TI_w$diff_EV_ST_pct=(impact_TI_w$'Both'-impact_TI_w$'Early Vaccine and Staff Tests Daily Antigen')/impact_TI_w$'Both'
impact_TI_w$diff_EV_ST_SV=impact_TI_w$'Both'-impact_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine'
impact_TI_w$diff_EV_ST_SV_pct=(impact_TI_w$'Both'-impact_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/impact_TI_w$'Both'
impact_TI_w$diff_ALL=(impact_TI_w$'No Masks or Tests'-impact_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')
impact_TI_w$diff_ALL_pct=(impact_TI_w$'No Masks or Tests'-impact_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/impact_TI_w$'No Masks or Tests'

impact_TS_w$diff_test=impact_TS_w$None-impact_TS_w$Both
impact_TS_w$diff_mask=impact_TS_w$'No Masks'-impact_TS_w$'High Masks'
impact_TS_w$diff_test_pct=(impact_TS_w$None-impact_TS_w$Both)/impact_TS_w$None
impact_TS_w$diff_mask_pct=(impact_TS_w$'No Masks'-impact_TS_w$'High Masks')/impact_TS_w$'No Masks'
impact_TS_w$diff_staff_vacc=impact_TS_w$'Both'-impact_TS_w$'Staff Vaccine'
impact_TS_w$diff_staff_vacc_pct=(impact_TS_w$'Both'-impact_TS_w$'Staff Vaccine')/impact_TS_w$'Both'
impact_TS_w$diff_vacc=impact_TS_w$'No Vaccine'-impact_TS_w$'High Vaccine'
impact_TS_w$diff_vacc_pct=(impact_TS_w$'No Vaccine'-impact_TS_w$'High Vaccine')/impact_TS_w$'No Vaccine'
impact_TS_w$diff_early_vacc=impact_TS_w$'Both'-impact_TS_w$'Early Vaccine'
impact_TS_w$diff_early_vacc_pct=(impact_TS_w$'Both'-impact_TS_w$'Early Vaccine')/impact_TS_w$'Both'
impact_TS_w$diff_staff_daily_test=impact_TS_w$'Both'-impact_TS_w$'Staff Tests Daily Antigen'
impact_TS_w$diff_staff_daily_test_pct=(impact_TS_w$'Both'-impact_TS_w$'Staff Tests Daily Antigen')/impact_TS_w$'Both'
impact_TS_w$diff_staff_weekly_test=impact_TS_w$'Both'-impact_TS_w$'Staff Tests'
impact_TS_w$diff_staff_weekly_test_pct=(impact_TS_w$'Both'-impact_TS_w$'Staff Tests')/impact_TS_w$'Both'
impact_TS_w$diff_EV_ST=impact_TS_w$'Both'-impact_TS_w$'Early Vaccine and Staff Tests Daily Antigen'
impact_TS_w$diff_EV_ST_pct=(impact_TS_w$'Both'-impact_TS_w$'Early Vaccine and Staff Tests Daily Antigen')/impact_TS_w$'Both'
impact_TS_w$diff_EV_ST_SV=impact_TS_w$'Both'-impact_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine'
impact_TS_w$diff_EV_ST_SV_pct=(impact_TS_w$'Both'-impact_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/impact_TS_w$'Both'
impact_TS_w$diff_ALL=(impact_TS_w$'No Masks or Tests'-impact_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')
impact_TS_w$diff_ALL_pct=(impact_TS_w$'No Masks or Tests'-impact_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/impact_TS_w$'No Masks or Tests'

impact_TI_l=impact_TI_w %>% select(Start_date,pct_inf_prior_2week,diff_test,diff_mask,diff_vacc,diff_mask_pct,
                                   diff_test_pct,diff_vacc,diff_vacc_pct,diff_staff_vacc,diff_staff_vacc_pct,diff_early_vacc,
                                   diff_early_vacc_pct,diff_staff_daily_test,diff_staff_daily_test_pct,diff_staff_weekly_test,
                                   diff_staff_weekly_test_pct,diff_EV_ST,diff_EV_ST_pct,diff_EV_ST_SV,diff_EV_ST_SV_pct,diff_ALL,diff_ALL_pct) %>% 
  pivot_longer(cols=starts_with("diff"),names_to="Intervention",values_to="impact_pct")
impact_TS_l=impact_TS_w %>% select(Start_date,pct_inf_prior_2week,diff_test,diff_mask,diff_vacc,diff_mask_pct,
                                   diff_test_pct,diff_vacc,diff_vacc_pct,diff_staff_vacc,diff_staff_vacc_pct,diff_early_vacc,
                                   diff_early_vacc_pct,diff_staff_daily_test,diff_staff_daily_test_pct,diff_staff_weekly_test,
                                   diff_staff_weekly_test_pct,diff_EV_ST,diff_EV_ST_pct,diff_EV_ST_SV,diff_EV_ST_SV_pct,diff_ALL,diff_ALL_pct) %>% 
  pivot_longer(cols=starts_with("diff"),names_to="Intervention",values_to="impact_pct")

###Add import_inf_pct (from baseline scenario)
import_dat<-combo_dat %>% filter(testing=="Both") %>% group_by(Start_date) %>% summarize(import_inf_pct=mean(import_inf_pct))

###Combo plot
impact_TS_l$measure="symptomatic"
impact_TI_l$measure="infection"
impact_l_c=impact_TI_l %>% bind_rows(impact_TS_l)
impact_l_c=impact_l_c %>% left_join(import_dat,by="Start_date")
impact_l_pv=impact_l_c[impact_l_c$Start_date>"2021-05-31",]
impact_TS_l_pv=impact_TS_l[impact_TS_l$Start_date>"2021-05-31",]
impact_TI_l_pv=impact_TI_l[impact_TI_l$Start_date>"2021-05-31",]
impact_TS_l_pv=impact_TS_l_pv %>% left_join(import_dat,by="Start_date")
impact_TI_l_pv=impact_TI_l_pv %>% left_join(import_dat,by="Start_date")
##Summarize impact
impact_sum=impact_l_pv %>% group_by(Intervention,measure) %>% summarize(impact_pct_m=mean(impact_pct),
                                                                        impact_pct_med=median(impact_pct),
                                                                        pct_lo=quantile(impact_pct,.025),
                                                                        pct_hi=quantile(impact_pct,.975))
impact_sum=impact_sum[grepl("pct",impact_sum$Intervention),]

plot_impact_sum<-ggplot(impact_l_pv[impact_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                                    "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_staff_weekly_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),],
                        aes(x=as.factor(measure),y=impact_pct*100,fill=Intervention))+geom_boxplot(outlier.alpha = 0.3)+theme_bw()+
  coord_cartesian(ylim=c(-10,100))+scale_x_discrete(breaks=c("infection","symptomatic"),labels=c("Infections","Symptomatic infections"))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                                   "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_staff_weekly_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),
                    labels=c("Masks","Screening","Arrival\nvaccination","Staff\nvaccination","Early\nvaccination","Staff -\ndaily tests","Staff -\nweekly tests",
                             "Early vacc. +\nstaff tests","Early vacc. +\nstaff tests +\nstaff vacc.","All interventions"))+
  xlab("")+  ylab("Intervention impact (% reduction)") +geom_hline(yintercept=0,lty=2)

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_summary.pdf",plot_impact_sum,width=6,height=4,pointsize=14,units="in")

##Try interventions on y-axis - with colors being symptomatic vs. asymptomatic
plot_impact_sum_alt<-ggplot(impact_l_pv[impact_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                                        "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_staff_weekly_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),],
                            aes(x=as.factor(Intervention),y=impact_pct*100,fill=as.factor(measure)))+geom_boxplot(outlier.alpha = 0.3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  coord_cartesian(ylim=c(-10,100))+scale_fill_brewer(palette="Pastel1",name="",breaks=c("infection","symptomatic"),labels=c("Infections","Symptomatic\ninfections"))+
  scale_x_discrete(breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                            "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_staff_weekly_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),
                   labels=c("Masks","Screening","Arrival\nvaccination","Staff\nvaccination","Early\nvaccination","Staff -\ndaily tests","Staff -\nweekly tests",
                            "Early vacc. +\nstaff tests","Early vacc. +\nstaff tests +\nstaff vacc.","All interventions"))+
  xlab("")+  ylab("Intervention impact (% reduction)") +geom_hline(yintercept=0,lty=2)

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_summary_alt.pdf",plot_impact_sum_alt,width=6,height=4,pointsize=14,units="in")



plot_impact_sum_poster<-ggplot(impact_l_pv[impact_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct"),],
                               aes(x=as.factor(measure),y=impact_pct*100,fill=Intervention))+geom_boxplot(outlier.alpha = 0.3)+theme_bw()+
  coord_cartesian(ylim=c(-10,60))+scale_x_discrete(breaks=c("infection","symptomatic"),labels=c("Infections","Symptomatic infections"))+
  scale_fill_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                   "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_staff_weekly_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),
                      labels=c("Masks","Screening","Arrival\nvaccination","Staff\nvaccination","Early\nvaccination","Staff -\ndaily tests","Staff -\nweekly tests",
                               "Early vacc. +\nstaff tests","Early vacc. +\nstaff tests +\nstaff vacc.","All interventions"))+
  xlab("")+  ylab("Intervention impact (% reduction)") +geom_hline(yintercept=0,lty=2)

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_summary_poster.pdf",plot_impact_sum_poster,width=6,height=4,pointsize=14,units="in")


plot_impact_arr<-ggplot(impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),],aes(x=as.factor(round(pct_inf_prior_2week,2)),y=impact_pct*100,fill=(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+coord_cartesian(ylim=c(-5,100))+
  scale_x_discrete(breaks=seq(0,0.08,.01))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                    labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+
  ylab("Intervention impact\n(% reduction in symptomatic infections)")+
  geom_hline(aes(yintercept=tsm,color="diff_early_vacc_pct"),lty=2)+
  geom_hline(aes(yintercept=msm,color="diff_mask_pct"),lty=2)+
  geom_hline(aes(yintercept=evsm,color="diff_staff_daily_test_pct"),lty=2)+
  geom_hline(aes(yintercept=asm,color="diff_ALL_pct"),lty=2)+
  geom_hline(aes(yintercept=stsm,color="diff_test_pct"),lty=2)+guides(color = FALSE)


plot_int_imp<-ggplot(impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),],aes(x=as.factor(round(import_inf_pct,2)),y=impact_pct*100,fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+coord_cartesian(ylim=c(-5,100))+
  scale_x_discrete(breaks=seq(0,0.22,.02))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                    labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  xlab("Portion of staff infected off base")+
  ylab("Intervention impact\n(% reduction in symptomatic infections)")+
  geom_hline(aes(yintercept=tsm,color="diff_early_vacc_pct"),lty=2)+
  geom_hline(aes(yintercept=msm,color="diff_mask_pct"),lty=2)+
  geom_hline(aes(yintercept=evsm,color="diff_staff_daily_test_pct"),lty=2)+
  geom_hline(aes(yintercept=asm,color="diff_ALL_pct"),lty=2)+
  geom_hline(aes(yintercept=stsm,color="diff_test_pct"),lty=2)+guides(color = FALSE)

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_imports.pdf",plot_int_imp,width=6,height=4,pointsize=14,units="in")

plot_int_imp_poster<-ggplot(impact_TI_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct"),],aes(x=as.factor(round(import_inf_pct,2)),y=impact_pct*100,fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+coord_cartesian(ylim=c(-5,60))+
  scale_x_discrete(breaks=seq(0,0.22,.02))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct"),
                    labels=c("Masks","Arrival\nscreening","Arrival\nvaccination","Staff\nvaccination"))+
  xlab("Portion of staff infected off base")+
  ylab("Intervention impact\n(% reduction in symptomatic infections)")

plot_int_imp_poster<-ggplot(impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct"),],
                            aes(x=round(import_inf_pct,2),y=impact_pct*100,color=as.factor(Intervention),fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  theme_bw()+coord_cartesian(ylim=c(-5,60))+
  geom_jitter(size=2,alpha=0.6)+geom_smooth(method="lm",alpha=0.3)+
  #scale_x_discrete(breaks=seq(0,0.22,.02))+
  scale_color_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct"),
                       labels=c("Masks","Arrival\nscreening","Arrival\nvaccination","Staff\nvaccination"))+
  scale_fill_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct"),
                      labels=c("Masks","Arrival\nscreening","Arrival\nvaccination","Staff\nvaccination"))+
  xlab("Portion of staff infected off base")+
  ylab("Intervention impact\n(% reduction in symptomatic infections)")


ggsave("../figures/loAsymp/interventions_impact_post-vaccination_imports_poster.pdf",plot_int_imp_poster,width=6,height=4,pointsize=14,units="in")

# plot_pct_inter<-plot_grid(plot_int_arr,plot_int_imp,ncol=1,labels="AUTO")
# 
# ggsave("../figures/loAsymp/interventions_impact_post-vaccination_by_arrival_import.pdf",plot_pct_inter,width=8,height=10,pointsize=14,units="in")

impact_pct_arrival<-ggplot(impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),],
                           aes(x=pct_inf_prior_2week,y=impact_pct*100,color=as.factor(Intervention),fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  #geom_boxplot(outlier.alpha = 0.2)+
  geom_point(size=1,alpha=0.6)+geom_smooth(alpha=0.3,method="lm")+
  theme_bw()+  coord_cartesian(ylim=c(-5,100))+
  #scale_x_discrete(breaks=seq(0,0.08,.01))+
  scale_fill_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                      labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  scale_color_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                       labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  xlab("Portion of trainees infected prior to arrival")+
  ylab("Intervention impact\n(% reduction)")


impact_pct_import<-ggplot(impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),],
                          aes(x=import_inf_pct,y=impact_pct*100,color=as.factor(Intervention),fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  #geom_boxplot(outlier.alpha = 0.2)+
  geom_point(size=1,alpha=0.6)+geom_smooth(alpha=0.3,method="lm")+
  theme_bw()+  coord_cartesian(ylim=c(-5,100))+
  #scale_x_discrete(breaks=seq(0,0.08,.01))+
  scale_fill_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                      labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  scale_color_discrete(name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                       labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  xlab("Portion of staff infected off base")+
  ylab("Intervention impact\n(% reduction)")

plot_pct_inter<-plot_grid(impact_pct_arrival,impact_pct_import,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_by_arrival_import_UPD.pdf",plot_pct_inter,width=6,height=7.5,pointsize=14,units="in")


##Linear regression - test significance of impact of arrivals and imports
summary(lm(impact_pct~import_inf_pct,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask"),]))
summary(lm(impact_pct~import_inf_pct,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_early_vacc"),]))
summary(lm(impact_pct~import_inf_pct,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_staff_vacc"),]))
summary(lm(impact_pct~import_inf_pct,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_staff_daily_test"),]))
summary(lm(impact_pct~import_inf_pct,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_test"),]))
summary(lm(impact_pct~import_inf_pct,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_ALL"),]))

summary(lm(impact_pct~pct_inf_prior_2week,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_mask"),]))
summary(lm(impact_pct~pct_inf_prior_2week,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_early_vacc"),]))
summary(lm(impact_pct~pct_inf_prior_2week,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_staff_vacc"),]))
summary(lm(impact_pct~pct_inf_prior_2week,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_staff_daily_test"),]))
summary(lm(impact_pct~pct_inf_prior_2week,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_test"),]))
summary(lm(impact_pct~pct_inf_prior_2week,
           data=impact_TS_l_pv[impact_TS_l_pv$Intervention %in% c("diff_ALL"),]))


### Plot ...
plot_test_1<-ggplot(combo_dat_test,aes(x=as.factor(round(import_inf_pct,2)),y=total_infections/Total_size,fill=testing))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+ylim(0,1)+scale_fill_brewer(palette="Accent",name="Arrival\nTesting",breaks=c("None","Both"),labels=c("No","Yes"))+
  scale_x_discrete(breaks=seq(0,0.2,.02))+ #scale_fill_discrete(name="Arrival Testing",breaks=c("None","Both"),labels=c("No","Yes"))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")
plot_test_1s<-ggplot(combo_dat_test,aes(x=as.factor(round(import_inf_pct,2)),y=(symptomatic_quarantine_infections+symptomatic_training_infections)/Total_size,fill=testing))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+ylim(0,0.4)+scale_fill_brewer(palette="Accent",name="Arrival\nTesting",breaks=c("None","Both"),labels=c("No","Yes"))+
  scale_x_discrete(breaks=seq(0,0.2,.02))+ #scale_fill_discrete(name="Arrival Testing",breaks=c("None","Both"),labels=c("No","Yes"))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")



##Calcultae diff between testing and no testing -Infections
dat_TI = combo_dat[,c(2:3,16,34,35)]
dat_TI_b=combo_dat[combo_dat$testing=="Both",c("rep","Start_date","import_inf_pct")]
dat_TI_w = dat_TI %>% pivot_wider(names_from=testing,values_from=total_inf_pct)
dat_TI_w$diff_test=dat_TI_w$None-dat_TI_w$Both
dat_TI_w$diff_mask=dat_TI_w$'No Masks'-dat_TI_w$'High Masks'
dat_TI_w$diff_test_pct=(dat_TI_w$None-dat_TI_w$Both)/dat_TI_w$None
dat_TI_w$diff_mask_pct=(dat_TI_w$'No Masks'-dat_TI_w$'High Masks')/dat_TI_w$'No Masks'
dat_TI_w$diff_staff_vacc=dat_TI_w$'Both'-dat_TI_w$'Staff Vaccine'
dat_TI_w$diff_staff_vacc_pct=(dat_TI_w$'Both'-dat_TI_w$'Staff Vaccine')/dat_TI_w$'Both'
dat_TI_w$diff_vacc=dat_TI_w$'No Vaccine'-dat_TI_w$'High Vaccine'
dat_TI_w$diff_vacc_pct=(dat_TI_w$'No Vaccine'-dat_TI_w$'High Vaccine')/dat_TI_w$'No Vaccine'
dat_TI_w$diff_early_vacc=dat_TI_w$'No Vaccine'-dat_TI_w$'Early Vaccine'
dat_TI_w$diff_early_vacc_pct=(dat_TI_w$'No Vaccine'-dat_TI_w$'Early Vaccine')/dat_TI_w$'No Vaccine'
dat_TI_w$diff_staff_daily_test=dat_TI_w$'Both'-dat_TI_w$'Staff Tests Daily Antigen'
dat_TI_w$diff_staff_daily_test_pct=(dat_TI_w$'Both'-dat_TI_w$'Staff Tests Daily Antigen')/dat_TI_w$'Both'
dat_TI_w$diff_staff_weekly_test=dat_TI_w$'Both'-dat_TI_w$'Staff Tests'
dat_TI_w$diff_staff_weekly_test_pct=(dat_TI_w$'Both'-dat_TI_w$'Staff Tests')/dat_TI_w$'Both'
dat_TI_w$diff_EV_ST=dat_TI_w$'No Vaccine'-dat_TI_w$'Early Vaccine and Staff Tests Daily Antigen'
dat_TI_w$diff_EV_ST_pct=(dat_TI_w$'No Vaccine'-dat_TI_w$'Early Vaccine and Staff Tests Daily Antigen')/dat_TI_w$'No Vaccine'
dat_TI_w$diff_EV_ST_SV=dat_TI_w$'No Vaccine'-dat_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine'
dat_TI_w$diff_EV_ST_SV_pct=(dat_TI_w$'No Vaccine'-dat_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/dat_TI_w$'No Vaccine'
dat_TI_w$diff_ALL_pct=(dat_TI_w$'No Masks or Tests'-dat_TI_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/dat_TI_w$'No Masks or Tests'

dat_TI_w = dat_TI_w %>% left_join(dat_TI_b,by=c('rep','Start_date'))

##Calcultae diff between testing and no testing - Symptomatic Infections
dat_TS = combo_dat[,c(2:3,16,34,36)]
dat_TS_b=combo_dat[combo_dat$testing=="Both",c("rep","Start_date","import_inf_pct")]
dat_TS_w = dat_TS %>% pivot_wider(names_from=testing,values_from=total_symp_inf_pct)
dat_TS_w$diff_test=dat_TS_w$None-dat_TS_w$Both
dat_TS_w$diff_test_pct=(dat_TS_w$None-dat_TS_w$Both)/dat_TS_w$None
dat_TS_w$diff_mask=dat_TS_w$'No Masks'-dat_TS_w$'High Masks'
dat_TS_w$diff_mask_pct=(dat_TS_w$'No Masks'-dat_TS_w$'High Masks')/dat_TS_w$'No Masks'
dat_TS_w$diff_staff_vacc=dat_TS_w$'Both'-dat_TS_w$'Staff Vaccine'
dat_TS_w$diff_staff_vacc_pct=(dat_TS_w$'Both'-dat_TS_w$'Staff Vaccine')/dat_TS_w$'Both'
dat_TS_w$diff_vacc=dat_TS_w$'No Vaccine'-dat_TS_w$'High Vaccine'
dat_TS_w$diff_vacc_pct=(dat_TS_w$'No Vaccine'-dat_TS_w$'High Vaccine')/dat_TS_w$'No Vaccine'
dat_TS_w$diff_early_vacc=dat_TS_w$'No Vaccine'-dat_TS_w$'Early Vaccine'
dat_TS_w$diff_early_vacc_pct=(dat_TS_w$'No Vaccine'-dat_TS_w$'Early Vaccine')/dat_TS_w$'No Vaccine'
dat_TS_w$diff_staff_daily_test=dat_TS_w$'Both'-dat_TS_w$'Staff Tests Daily Antigen'
dat_TS_w$diff_staff_daily_test_pct=(dat_TS_w$'Both'-dat_TS_w$'Staff Tests Daily Antigen')/dat_TS_w$'Both'
dat_TS_w$diff_staff_weekly_test=dat_TS_w$'Both'-dat_TS_w$'Staff Tests'
dat_TS_w$diff_staff_weekly_test_pct=(dat_TS_w$'Both'-dat_TS_w$'Staff Tests')/dat_TS_w$'Both'
dat_TS_w$diff_EV_ST=dat_TS_w$'No Vaccine'-dat_TS_w$'Early Vaccine and Staff Tests Daily Antigen'
dat_TS_w$diff_EV_ST_pct=(dat_TS_w$'No Vaccine'-dat_TS_w$'Early Vaccine and Staff Tests Daily Antigen')/dat_TS_w$'No Vaccine'
dat_TS_w$diff_EV_ST_SV=dat_TS_w$'No Vaccine'-dat_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine'
dat_TS_w$diff_EV_ST_SV_pct=(dat_TS_w$'No Vaccine'-dat_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/dat_TS_w$'No Vaccine'
dat_TS_w$diff_ALL_pct=(dat_TS_w$'No Masks or Tests'-dat_TS_w$'Early Vaccine and Staff Tests Daily Antigen and Staff Vaccine')/dat_TS_w$'No Masks or Tests'

dat_TS_w = dat_TS_w %>% left_join(dat_TS_b,by=c('rep','Start_date'))


#Tests
##PLOT A - IMPACT OF TESTING
plot_test_1<-ggplot(combo_dat_test,aes(x=as.factor(round(import_inf_pct,2)),y=total_infections/Total_size,fill=testing))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+ylim(0,1)+scale_fill_brewer(palette="Accent",name="Arrival\nTesting",breaks=c("None","Both"),labels=c("No","Yes"))+
  scale_x_discrete(breaks=seq(0,0.2,.02))+ #scale_fill_discrete(name="Arrival Testing",breaks=c("None","Both"),labels=c("No","Yes"))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")
plot_test_1s<-ggplot(combo_dat_test,aes(x=as.factor(round(import_inf_pct,2)),y=(symptomatic_quarantine_infections+symptomatic_training_infections)/Total_size,fill=testing))+
  #geom_jitter(alpha=0.1)+theme_bw()+geom_smooth()
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+ylim(0,0.4)+scale_fill_brewer(palette="Accent",name="Arrival\nTesting",breaks=c("None","Both"),labels=c("No","Yes"))+
  scale_x_discrete(breaks=seq(0,0.2,.02))+ #scale_fill_discrete(name="Arrival Testing",breaks=c("None","Both"),labels=c("No","Yes"))+
  xlab("Portion of staff infected off base")+ylab("Portion of individuals\n infected during training camp") #+geom_smooth(method="lm")

dat_TI_l=dat_TI_w %>% select(rep,Start_date,pct_inf_prior_2week,import_inf_pct,diff_test,diff_mask,diff_vacc,diff_mask_pct,
                             diff_test_pct,diff_vacc,diff_vacc_pct,diff_staff_vacc,diff_staff_vacc_pct,diff_early_vacc,
                             diff_early_vacc_pct,diff_staff_daily_test,diff_staff_daily_test_pct,diff_staff_weekly_test,
                             diff_staff_weekly_test_pct,diff_EV_ST,diff_EV_ST_pct,diff_EV_ST_SV,diff_EV_ST_SV_pct,diff_ALL_pct) %>% 
  pivot_longer(cols=starts_with("diff"),names_to="Intervention",values_to="impact_pct")
dat_TS_l=dat_TS_w %>% select(rep,Start_date,pct_inf_prior_2week,import_inf_pct,diff_test,diff_mask,diff_vacc,diff_mask_pct,
                             diff_test_pct,diff_vacc,diff_vacc_pct,diff_staff_vacc,diff_staff_vacc_pct,diff_early_vacc,
                             diff_early_vacc_pct,diff_staff_daily_test,diff_staff_daily_test_pct,diff_staff_weekly_test,
                             diff_staff_weekly_test_pct,diff_EV_ST,diff_EV_ST_pct,diff_EV_ST_SV,diff_EV_ST_SV_pct,diff_ALL_pct) %>% 
  pivot_longer(cols=starts_with("diff"),names_to="Intervention",values_to="impact_pct")

#IMPACT TESTING AND MASKING BASED ON PCT INFECTED
plot_test_3<-ggplot(dat_TI_l[dat_TI_l$Intervention %in% c("diff_mask","diff_test"),],aes(x=as.factor(round(pct_inf_prior_2week*2,2)/2),y=impact_pct,fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+
  scale_x_discrete(breaks=seq(0,0.08,.01))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask","diff_test"),labels=c("Masks","Testing"))+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+ylab("Intervention impact\n(difference in portion infected)") #+geom_smooth(method="lm")
#IMPACT OF TESTING AND MASKING BASED ON IMPORTS
#Masks
plot_test_2<-ggplot(dat_TI_l[dat_TI_l$Intervention %in% c("diff_mask","diff_test"),],aes(x=as.factor(round(import_inf_pct,2)),y=impact_pct,fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+
  scale_x_discrete(breaks=seq(0,0.2,.02))+ #scale_fill_discrete(name="Arrival Testing",breaks=c("None","Both"),labels=c("No","Yes"))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask","diff_test"),labels=c("Masks","Testing"))+
  xlab("Portion of staff infected off base")+ylab("Intervention impact\n(difference in portion infected)") #+geom_smooth(method="lm")

plot_intervention<-plot_grid(plot_test_1,plot_test_2,plot_test_3,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/infections_by_arrival_import_impact.pdf",plot_intervention,width=8,height=14,pointsize=14,units="in")

ggplot(dat_TI_l[dat_TI_l$Intervention %in% c("diff_mask","diff_test","diff_vacc","diff_staff_vacc"),],aes(x=as.factor(round(import_inf_pct,2)),y=impact_pct,fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+
  scale_x_discrete(breaks=seq(0,0.2,.02))+ #scale_fill_discrete(name="Arrival Testing",breaks=c("None","Both"),labels=c("No","Yes"))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask","diff_test"),labels=c("Masks","Testing"))+
  xlab("Portion of staff infected off base")+ylab("Intervention impact\n(difference in portion infected)") #+geom_smooth(method="lm")

###Combo plot
dat_TS_l$measure="symptomatic"
dat_TI_l$measure="infection"
dat_l_c=dat_TI_l %>% bind_rows(dat_TS_l)
dat_l_pv=dat_l_c[dat_l_c$Start_date>"2021-05-31",]
dat_TS_l_pv=dat_TS_l[dat_TS_l$Start_date>"2021-05-31",]
dat_TI_l_pv=dat_TI_l[dat_TI_l$Start_date>"2021-05-31",]

##Assess mean impact per week
dat_l_c$impact_pct=ifelse(is.infinite(dat_l_c$impact_pct),NA,dat_l_c$impact_pct)
dat_l_c_sum = dat_l_c %>% group_by(Start_date,Intervention,measure) %>% 
  summarize(impact_sum=sum(impact_pct,na.rm=T)/200,
            impact_mean=mean(impact_pct,na.rm=T),
            impact_median=median(impact_pct,na.rm=T),
            import_inf_pct=mean(import_inf_pct ,na.rm=T),
            pct_inf_prior_2week =mean(pct_inf_prior_2week ,na.rm=T))


dat_l_pv=dat_l_c[dat_l_c$Start_date>"2021-05-31",]
dat_l_pv_sum = dat_l_pv %>% group_by(Start_date,Intervention,measure) %>% 
  summarize(impact_sum=sum(impact_pct,na.rm=T)/200,
            impact_mean=mean(impact_pct,na.rm=T),
            impact_median=median(impact_pct,na.rm=T),
            import_inf_pct=mean(import_inf_pct ,na.rm=T),
            pct_inf_prior_2week =mean(pct_inf_prior_2week ,na.rm=T))

msm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_mask_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)
tsm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_test_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)
vsm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_vacc_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)
svsm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_staff_vacc_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)
evsm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_early_vacc_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)
stsm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_staff_daily_test_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)
asm=median(dat_TS_l_pv$impact_pct[dat_TS_l_pv$Intervention=="diff_ALL_pct"&!is.infinite(dat_TS_l_pv$impact_pct)]*100,na.rm=T)

dat_l_sum =dat_l_pv[!is.infinite(dat_l_pv$impact_pct),] %>% group_by(Intervention,measure) %>% 
  summarize(impact_m=mean(impact_pct,na.rm=T),
            impact_med=median(impact_pct,na.rm=T),
            impact_lo=quantile(impact_pct,c(0.025),na.rm=T),impact_hi=quantile(impact_pct,c(0.975),na.rm=T))

plot_int_sum<-ggplot(dat_l_pv[dat_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                           "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),],
                     aes(x=as.factor(measure),y=impact_pct*100,fill=Intervention))+geom_boxplot(outlier.alpha = 0.1)+theme_bw()+
  coord_cartesian(ylim=c(-100,100))+scale_x_discrete(breaks=c("infection","symptomatic"),labels=c("Infections","Symptomatic infections"))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                                   "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),
                    labels=c("Masks","Screening","Arrival\nvaccination","Staff\nvaccination","Early\nvaccination","Staff -\ndaily tests",
                             "Early vacc. +\nstaff tests","Early vacc. +\nstaff tests +\nstaff vacc.","All interventions"))+
  xlab("")+  ylab("Intervention impact (% reduction)") #+geom_smooth(method="lm")

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_update.pdf",plot_int_sum,width=6,height=4,pointsize=14,units="in")

###
plot_int_sum2<-ggplot(dat_l_pv_sum[dat_l_pv_sum$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                                    "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),],
                      aes(x=as.factor(measure),y=impact_pct*100,fill=Intervention))+geom_boxplot(outlier.alpha = 0.1)+theme_bw()+
  coord_cartesian(ylim=c(-100,100))+scale_x_discrete(breaks=c("infection","symptomatic"),labels=c("Infections","Symptomatic infections"))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_vacc_pct","diff_staff_vacc_pct",
                                                                   "diff_early_vacc_pct","diff_staff_daily_test_pct","diff_EV_ST_pct","diff_EV_ST_SV_pct","diff_ALL_pct"),
                    labels=c("Masks","Screening","Arrival\nvaccination","Staff\nvaccination","Early\nvaccination","Staff -\ndaily tests",
                             "Early vacc. +\nstaff tests","Early vacc. +\nstaff tests +\nstaff vacc.","All interventions"))+
  xlab("")+  ylab("Intervention impact (% reduction)") #+geom_smooth(method="lm")

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_summary.pdf",plot_int_sum2,width=6,height=4,pointsize=14,units="in")

###
##For some reason line colors don't match the scale fills - so lines have been assigned color that doesn't match name but is correct value
plot_int_arr<-ggplot(dat_TS_l_pv[dat_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),],aes(x=as.factor(round(pct_inf_prior_2week,2)),y=impact_pct*100,fill=(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+coord_cartesian(ylim=c(-120,100))+
  scale_x_discrete(breaks=seq(0,0.08,.01))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                    labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  xlab("Portion of trainees infected prior to arrival (2 weeks)")+
  ylab("Intervention impact\n(% reduction in symptomatic infections)")+
  geom_hline(aes(yintercept=tsm,color="diff_early_vacc_pct"),lty=2)+
  geom_hline(aes(yintercept=msm,color="diff_mask_pct"),lty=2)+
  geom_hline(aes(yintercept=evsm,color="diff_staff_daily_test_pct"),lty=2)+
  geom_hline(aes(yintercept=asm,color="diff_ALL_pct"),lty=2)+
  geom_hline(aes(yintercept=stsm,color="diff_test_pct"),lty=2)+guides(color = FALSE)


plot_int_imp<-ggplot(dat_TS_l_pv[dat_TS_l_pv$Intervention %in% c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),],aes(x=as.factor(round(import_inf_pct/2,2)*2),y=impact_pct*100,fill=as.factor(Intervention)))+geom_hline(yintercept=0,lty=2)+
  geom_boxplot(outlier.alpha = 0.2)+theme_bw()+coord_cartesian(ylim=c(-120,100))+
  scale_x_discrete(breaks=seq(0,0.22,.02))+
  scale_fill_brewer(palette="Pastel1",name="Intervention",breaks=c("diff_mask_pct","diff_test_pct","diff_ALL_pct","diff_staff_daily_test_pct","diff_early_vacc_pct"),
                    labels=c("Masks","Arrival\nscreening","All","Staff\ntesting","Early\nvaccination"))+
  xlab("Portion of staff infected off base")+
  ylab("Intervention impact\n(% reduction in symptomatic infections)")+
  geom_hline(aes(yintercept=tsm,color="diff_early_vacc_pct"),lty=2)+
  geom_hline(aes(yintercept=msm,color="diff_mask_pct"),lty=2)+
  geom_hline(aes(yintercept=evsm,color="diff_staff_daily_test_pct"),lty=2)+
  geom_hline(aes(yintercept=asm,color="diff_ALL_pct"),lty=2)+
  geom_hline(aes(yintercept=stsm,color="diff_test_pct"),lty=2)+guides(color = FALSE)

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_imports.pdf",plot_int_imp,width=6,height=4,pointsize=14,units="in")

plot_pct_inter<-plot_grid(plot_int_arr,plot_int_imp,ncol=1,labels="AUTO")

ggsave("../figures/loAsymp/interventions_impact_post-vaccination_by_arrival_import.pdf",plot_pct_inter,width=8,height=10,pointsize=14,units="in")


pdf(file="../figures/loAsymp/masking_impact_by_imports.pdf",width=8,height=6,pointsize=14)
ggplot(combo_dat[combo_dat$testing %in% c("Both","No Masks"),],aes(x=pct_inf_prior_2week,y=(quarantine_infections+training_infections)/Total_size,color=testing))+
  geom_jitter(alpha=0.1)+theme_bw()+
  xlab("Infected past 2 weeks")+ylab("Portion of recruits infected during training camp")+geom_smooth(method="lm") 
dev.off()

