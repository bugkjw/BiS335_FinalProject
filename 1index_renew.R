library(knitr)
library(survival)
library(dplyr)
clin <- readRDS("clinical.Rds")
mut <- readRDS("mutation.Rds")
gex <- readRDS("expression.Rds")

data_ind1<-na.omit(clin[clin$survival_index == sort(unique(clin$survival_index))[1],])
data_ind2<-na.omit(clin[clin$survival_index == sort(unique(clin$survival_index))[2],])
data_ind3<-na.omit(clin[clin$survival_index == sort(unique(clin$survival_index))[3],])
data_ind4<-na.omit(clin[clin$survival_index == sort(unique(clin$survival_index))[4],])

data_ind1<-data_ind1[order(data_ind1$survival_time),]
data_ind2<-data_ind2[order(data_ind2$survival_time),]
data_ind3<-data_ind3[order(data_ind3$survival_time),]
data_ind4<-data_ind4[order(data_ind4$survival_time),]


Median <- matrix(nrow=1, ncol=4)
colnames(Median)<-c(1, 2, 3, 4)
rownames(Median)<-c("survival time")

Quantile <- matrix(nrow=1, ncol=4)
rownames(Quantile)<-c("survival time")
colnames(Quantile)<-c(1,2,3,4)

data_ind1_all1 <- filter(data_ind1, sample_id%in%intersect(mut[,1],colnames(gex)))
data_ind1_all2 <- filter(data_ind2, sample_id%in%intersect(mut[,1],colnames(gex)))
data_ind1_all3 <- filter(data_ind3, sample_id%in%intersect(mut[,1],colnames(gex)))
data_ind1_all4 <- filter(data_ind4, sample_id%in%intersect(mut[,1],colnames(gex)))
clin_all<-rbind(data_ind1_all1,data_ind1_all2,data_ind1_all3,data_ind1_all4)

p_val<-matrix(nrow=1, ncol=2)

#survival_index 1
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_ind1$sample_id)+1-rev(which(data_ind1$vital_status==1))
deaths<-(c(0,rev(which(data_ind1$vital_status==1)))-c(rev(which(data_ind1$vital_status==1)),0))[-c(1)]

data_ind1$Number_alive<-c(rep(number_alive[1],length(data_ind1$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_ind1$Deaths<-c(rep(0,length(data_ind1$sample_id)))
data_ind1$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_ind1$prop_surv<-rev((data_ind1$Number_alive-data_ind1$Deaths)/data_ind1$Number_alive)
data_ind1$prop_surv[is.na(data_ind1$prop_surv)]<-1
data_ind1$cum_prop_surv<-cumprod(data_ind1$prop_surv)
#median survival time
if(is.na(which(data_ind1$cum_prop_surv<=0.75)[1])){
  Median[1]<-c("No median")
  Quantile[1]<-c("No Quantile")}else
    if(is.na(which(data_ind1$cum_prop_surv<=0.5)[1])){
      Median[1]<-c("No median")
      Quantile[1]<-data_ind1$survival_time[which(data_ind1$cum_prop_surv<=0.75)[1]]
    }else
      if(data_ind1$vital_status[which(data_ind1$cum_prop_surv<=0.5)[1]]==1){
        Median[1]<-data_ind1$survival_time[which(data_ind1$cum_prop_surv<=0.5)[1]]
        Quantile[1]<-c("Don't need")}else 
          if(data_ind1$vital_status[which(data_ind1$cum_prop_surv<=0.75)[1]]==1){
            Quantile[1]<-data_ind1$survival_time[which(data_ind1$cum_prop_surv<=0.75)[1]]}else{
              Median[1]<-c("No median")
              Quantile[1]<-c("No Quantile")}
###########################

#survival_index 2
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_ind2$sample_id)+1-rev(which(data_ind2$vital_status==1))
deaths<-(c(0,rev(which(data_ind2$vital_status==1)))-c(rev(which(data_ind2$vital_status==1)),0))[-c(1)]

data_ind2$Number_alive<-c(rep(number_alive[1],length(data_ind2$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_ind2$Deaths<-c(rep(0,length(data_ind2$sample_id)))
data_ind2$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_ind2$prop_surv<-rev((data_ind2$Number_alive-data_ind2$Deaths)/data_ind2$Number_alive)
data_ind2$prop_surv[is.na(data_ind2$prop_surv)]<-1
data_ind2$cum_prop_surv<-cumprod(data_ind2$prop_surv)
#median survival time
if(is.na(which(data_ind2$cum_prop_surv<=0.75)[1])){
  Median[2]<-c("No median")
  Quantile[2]<-c("No Quantile")}else
    if(is.na(which(data_ind2$cum_prop_surv<=0.5)[1])){
      Median[2]<-c("No median")
      Quantile[2]<-data_ind2$survival_time[which(data_ind2$cum_prop_surv<=0.75)[1]]
    }else
      if(data_ind2$vital_status[which(data_ind2$cum_prop_surv<=0.5)[1]]==1){
        Median[2]<-data_ind2$survival_time[which(data_ind2$cum_prop_surv<=0.5)[1]]
        Quantile[2]<-c("Don't need")}else 
          if(data_ind2$vital_status[which(data_ind2$cum_prop_surv<=0.75)[1]]==1){
            Quantile[2]<-data_ind2$survival_time[which(data_ind2$cum_prop_surv<=0.75)[1]]}else{
              Median[2]<-c("No median")
              Quantile[2]<-c("No Quantile")}
###########################

#survival_index 3
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_ind3$sample_id)+1-rev(which(data_ind3$vital_status==1))
deaths<-(c(0,rev(which(data_ind3$vital_status==1)))-c(rev(which(data_ind3$vital_status==1)),0))[-c(1)]

data_ind3$Number_alive<-c(rep(number_alive[1],length(data_ind3$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_ind3$Deaths<-c(rep(0,length(data_ind3$sample_id)))
data_ind3$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_ind3$prop_surv<-rev((data_ind3$Number_alive-data_ind3$Deaths)/data_ind3$Number_alive)
data_ind3$prop_surv[is.na(data_ind3$prop_surv)]<-1
data_ind3$cum_prop_surv<-cumprod(data_ind3$prop_surv)
#median survival time
if(is.na(which(data_ind3$cum_prop_surv<=0.75)[1])){
  Median[3]<-c("No median")
  Quantile[3]<-c("No Quantile")}else
    if(is.na(which(data_ind3$cum_prop_surv<=0.5)[1])){
      Median[3]<-c("No median")
      Quantile[3]<-data_ind3$survival_time[which(data_ind3$cum_prop_surv<=0.75)[1]]
    }else
      if(data_ind3$vital_status[which(data_ind3$cum_prop_surv<=0.5)[1]]==1){
        Median[3]<-data_ind3$survival_time[which(data_ind3$cum_prop_surv<=0.5)[1]]
        Quantile[3]<-c("Don't need")}else 
          if(data_ind3$vital_status[which(data_ind3$cum_prop_surv<=0.75)[1]]==1){
            Quantile[3]<-data_ind3$survival_time[which(data_ind3$cum_prop_surv<=0.75)[1]]}else{
              Median[3]<-c("No median")
              Quantile[3]<-c("No Quantile")}
###########################


#survival_index 4
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_ind4$sample_id)+1-rev(which(data_ind4$vital_status==1))
deaths<-(c(0,rev(which(data_ind4$vital_status==1)))-c(rev(which(data_ind4$vital_status==1)),0))[-c(1)]

data_ind4$Number_alive<-c(rep(number_alive[1],length(data_ind4$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_ind4$Deaths<-c(rep(0,length(data_ind4$sample_id)))
data_ind4$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_ind4$prop_surv<-rev((data_ind4$Number_alive-data_ind4$Deaths)/data_ind4$Number_alive)
data_ind4$prop_surv[is.na(data_ind4$prop_surv)]<-1
data_ind4$cum_prop_surv<-cumprod(data_ind4$prop_surv)
#median survival time
if(is.na(which(data_ind4$cum_prop_surv<=0.75)[1])){
  Median[4]<-c("No median")
  Quantile[4]<-c("No Quantile")}else
    if(is.na(which(data_ind4$cum_prop_surv<=0.5)[1])){
      Median[4]<-c("No median")
      Quantile[4]<-data_ind4$survival_time[which(data_ind4$cum_prop_surv<=0.75)[1]]
    }else
      if(data_ind4$vital_status[which(data_ind4$cum_prop_surv<=0.5)[1]]==1){
        Median[4]<-data_ind4$survival_time[which(data_ind4$cum_prop_surv<=0.5)[1]]
        Quantile[4]<-c("Don't need")}else 
          if(data_ind4$vital_status[which(data_ind4$cum_prop_surv<=0.75)[1]]==1){
            Quantile[4]<-data_ind4$survival_time[which(data_ind4$cum_prop_surv<=0.75)[1]]}else{
              Median[4]<-c("No median")
              Quantile[4]<-c("No Quantile")}
###########################
Median
Quantile


#graph
library(survminer)
ev<-1*(clin$vital_status==1)
fut<-as.numeric(clin$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin$survival_index,data=clin)
ggsurvplot(fit,data=clin,pval=TRUE)
fit<-survfit(su~clin$survival_index,data=clin)
#find p-value
survdiff(su~clin$survival_index,data=clin)
sdf<-survdiff(su~clin$survival_index,data=clin)
p_val[1]<-1-pchisq(sdf$chisq, length(sdf$n)-1)


#With all data
ev<-1*(clin_all$vital_status==1)
fut<-as.numeric(clin_all$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin_all$survival_index,data=clin_all)
ggsurvplot(fit,data=clin_all,pval=TRUE)
fit<-survfit(su~clin_all$survival_index,data=clin_all)
#find p-value
survdiff(su~clin_all$survival_index,data=clin_all)
sdf<-survdiff(su~clin_all$survival_index,data=clin_all)
p_val[2]<-1-pchisq(sdf$chisq, length(sdf$n)-1)