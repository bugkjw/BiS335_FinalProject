library(knitr)
install.packages("survminer")
clin <- readRDS("C:/Users/Shin Dongho/Desktop/Finalterm-Project/clinical.Rds")
mut <- readRDS("C:/Users/Shin Dongho/Desktop/Finalterm-Project/mutation.Rds")
expre <- readRDS("C:/Users/Shin Dongho/Desktop/Finalterm-Project/expression.Rds")

data_st1<-clin[clin$stage == sort(unique(clin$stage))[1],]
data_st2<-clin[clin$stage == sort(unique(clin$stage))[2],]
data_st3<-clin[clin$stage == sort(unique(clin$stage))[3],]

data_st1<-data_st1[order(data_st1$survival_time),]
data_st2<-data_st2[order(data_st2$survival_time),]
data_st3<-data_st3[order(data_st3$survival_time),]

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
data_st1$survival_time[which(data_st1$cum_prop_surv<=0.5&data_st1$cum_prop_surv!=0)[1]-1]
###########################

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
data_st2$survival_time[which(data_st2$cum_prop_surv<=0.5&data_st2$cum_prop_surv!=0)[1]-1]
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
data_st3$survival_time[which(data_st3$cum_prop_surv<=0.5&data_st3$cum_prop_surv!=0)[1]-1]
###########################

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
