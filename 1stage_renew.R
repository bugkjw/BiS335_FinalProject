library(knitr)
library(survival)
library(dplyr)
clin <- readRDS("clinical.Rds")
mut <- readRDS("mutation.Rds")
gex <- readRDS("expression.Rds")

data_st1<-clin[clin$stage == sort(unique(clin$stage))[1],]
data_st2<-clin[clin$stage == sort(unique(clin$stage))[2],]
data_st3<-clin[clin$stage == sort(unique(clin$stage))[3],]

data_st1<-data_st1[order(data_st1$survival_time),]
data_st2<-data_st2[order(data_st2$survival_time),]
data_st3<-data_st3[order(data_st3$survival_time),]

Median <- matrix(nrow=1, ncol=length(unique(clin$stage)))
colnames(Median)<-c("stage i", "stage ii", "stage iii")
rownames(Median)<-c("survival time")

Quantile <- matrix(nrow=1, ncol=length(unique(clin$stage)))
rownames(Quantile)<-c("survival time")
colnames(Quantile)<-c("stage i", "stage ii", "stage iii")

data_st1_all1 <- filter(data_st1, sample_id%in%intersect(mut[,1],colnames(gex)))
data_st1_all2 <- filter(data_st2, sample_id%in%intersect(mut[,1],colnames(gex)))
data_st1_all3 <- filter(data_st3, sample_id%in%intersect(mut[,1],colnames(gex)))
clin_all<-rbind(data_st1_all1,data_st1_all2,data_st1_all3)

p_val<-matrix(nrow=1, ncol=2)

#stage 1
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_st1$sample_id)+1-rev(which(data_st1$vital_status==1))
deaths<-(c(0,rev(which(data_st1$vital_status==1)))-c(rev(which(data_st1$vital_status==1)),0))[-c(1)]

data_st1$Number_alive<-c(rep(number_alive[1],length(data_st1$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_st1$Deaths<-c(rep(0,length(data_st1$sample_id)))
data_st1$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_st1$prop_surv<-rev((data_st1$Number_alive-data_st1$Deaths)/data_st1$Number_alive)
data_st1$prop_surv[is.na(data_st1$prop_surv)]<-1
data_st1$cum_prop_surv<-cumprod(data_st1$prop_surv)
#median survival time
if(is.na(which(data_st1$cum_prop_surv<=0.75)[1])){
  Median[1]<-c("No median")
  Quantile[1]<-c("No Quantile")}else
    if(is.na(which(data_st1$cum_prop_surv<=0.5)[1])){
      Median[1]<-c("No median")
      Quantile[1]<-data_st1$survival_time[which(data_st1$cum_prop_surv<=0.75)[1]]
    }else
      if(data_st1$vital_status[which(data_st1$cum_prop_surv<=0.5)[1]]==1){
        Median[1]<-data_st1$survival_time[which(data_st1$cum_prop_surv<=0.5)[1]]
        Quantile[1]<-c("Don't need")}else 
          if(data_st1$vital_status[which(data_st1$cum_prop_surv<=0.75)[1]]==1){
            Quantile[1]<-data_st1$survival_time[which(data_st1$cum_prop_surv<=0.75)[1]]}else{
              Median[1]<-c("No median")
              Quantile[1]<-c("No Quantile")}

#stage 2
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_st2$sample_id)+1-rev(which(data_st2$vital_status==1))
deaths<-(c(0,rev(which(data_st2$vital_status==1)))-c(rev(which(data_st2$vital_status==1)),0))[-c(1)]

data_st2$Number_alive<-c(rep(number_alive[1],length(data_st2$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_st2$Deaths<-c(rep(0,length(data_st2$sample_id)))
data_st2$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_st2$prop_surv<-rev((data_st2$Number_alive-data_st2$Deaths)/data_st2$Number_alive)
data_st2$prop_surv[is.na(data_st2$prop_surv)]<-1
data_st2$cum_prop_surv<-cumprod(data_st2$prop_surv)
#median survival time
if(is.na(which(data_st2$cum_prop_surv<=0.75)[1])){
  Median[2]<-c("No median")
  Quantile[2]<-c("No Quantile")}else
    if(is.na(which(data_st2$cum_prop_surv<=0.5)[1])){
      Median[2]<-c("No median")
      Quantile[2]<-data_st2$survival_time[which(data_st2$cum_prop_surv<=0.75)[1]]
    }else
      if(data_st2$vital_status[which(data_st2$cum_prop_surv<=0.5)[1]]==1){
        Median[2]<-data_st2$survival_time[which(data_st2$cum_prop_surv<=0.5)[1]]
        Quantile[2]<-c("Don't need")}else 
          if(data_st2$vital_status[which(data_st2$cum_prop_surv<=0.75)[1]]==1){
            Quantile[2]<-data_st2$survival_time[which(data_st2$cum_prop_surv<=0.75)[1]]}else{
              Median[2]<-c("No median")
              Quantile[2]<-c("No Quantile")}
###########################

#stage 3
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_st3$sample_id)+1-rev(which(data_st3$vital_status==1))
deaths<-(c(0,rev(which(data_st3$vital_status==1)))-c(rev(which(data_st3$vital_status==1)),0))[-c(1)]

data_st3$Number_alive<-c(rep(number_alive[1],length(data_st3$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_st3$Deaths<-c(rep(0,length(data_st3$sample_id)))
data_st3$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_st3$prop_surv<-rev((data_st3$Number_alive-data_st3$Deaths)/data_st3$Number_alive)
data_st3$prop_surv[is.na(data_st3$prop_surv)]<-1
data_st3$cum_prop_surv<-cumprod(data_st3$prop_surv)
#median survival time
if(is.na(which(data_st3$cum_prop_surv<=0.75)[1])){
  Median[3]<-c("No median")
  Quantile[3]<-c("No Quantile")}else
    if(is.na(which(data_st3$cum_prop_surv<=0.5)[1])){
      Median[3]<-c("No median")
      Quantile[3]<-data_st3$survival_time[which(data_st3$cum_prop_surv<=0.75)[1]]
    }else
      if(data_st3$vital_status[which(data_st3$cum_prop_surv<=0.5)[1]]==1){
        Median[3]<-data_st3$survival_time[which(data_st3$cum_prop_surv<=0.5)[1]]
        Quantile[3]<-c("Don't need")}else 
          if(data_st3$vital_status[which(data_st3$cum_prop_surv<=0.75)[1]]==1){
            Quantile[3]<-data_st3$survival_time[which(data_st3$cum_prop_surv<=0.75)[1]]}else{
              Median[3]<-c("No median")
              Quantile[3]<-c("No Quantile")}
###########################
Median
Quantile

#graph
library(survminer)
ev<-1*(clin$vital_status==1)
fut<-as.numeric(clin$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin$stage,data=clin)
ggsurvplot(fit,data=clin,pval=TRUE)
fit<-survfit(su~clin$stage,data=clin)
#find p-value
survdiff(su~clin$stage,data=clin)
sdf<-survdiff(su~clin$stage, data=clin)
p_val[1]<-1-pchisq(sdf$chisq, length(sdf$n)-1)

#With all data
ev<-1*(clin_all$vital_status==1)
fut<-as.numeric(clin_all$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin_all$stage,data=clin_all)
ggsurvplot(fit,data=clin_all,pval=TRUE)
fit<-survfit(su~clin_all$stage,data=clin_all)
#find p-value
survdiff(su~clin_all$stage,data=clin_all)
sdf<-survdiff(su~clin_all$stage,data=clin_all)
p_val[2]<-1-pchisq(sdf$chisq, length(sdf$n)-1)


