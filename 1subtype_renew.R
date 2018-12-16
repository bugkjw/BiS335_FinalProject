#install.package("knitr")
#install.package("survival")
library(knitr)
library(survival)
library(dplyr)
clin <- readRDS("clinical.Rds")
mut <- readRDS("mutation.Rds")
gex <- readRDS("expression.Rds")

data_sub1<-na.omit(clin[clin$subtype == sort(unique(clin$subtype))[1],])
data_sub2<-na.omit(clin[clin$subtype == sort(unique(clin$subtype))[2],])
data_sub3<-na.omit(clin[clin$subtype == sort(unique(clin$subtype))[3],])
data_sub4<-na.omit(clin[clin$subtype == sort(unique(clin$subtype))[4],])
data_sub5<-na.omit(clin[clin$subtype == sort(unique(clin$subtype))[5],])

data_sub1<-data_sub1[order(data_sub1$survival_time),]
data_sub2<-data_sub2[order(data_sub2$survival_time),]
data_sub3<-data_sub3[order(data_sub3$survival_time),]
data_sub4<-data_sub4[order(data_sub4$survival_time),]
data_sub5<-data_sub5[order(data_sub5$survival_time),]

Median <- matrix(nrow=1, ncol=5)
colnames(Median)<-c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like")
rownames(Median)<-c("survival time")
Quantile <- matrix(nrow=1, ncol=5)
rownames(Quantile)<-c("survival time")
colnames(Quantile)<-c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like")

data_sub1_all1 <- filter(data_sub1, sample_id%in%intersect(mut[,1],colnames(gex)))
data_sub1_all2 <- filter(data_sub2, sample_id%in%intersect(mut[,1],colnames(gex)))
data_sub1_all3 <- filter(data_sub3, sample_id%in%intersect(mut[,1],colnames(gex)))
data_sub1_all4 <- filter(data_sub4, sample_id%in%intersect(mut[,1],colnames(gex)))
data_sub1_all5 <- filter(data_sub5, sample_id%in%intersect(mut[,1],colnames(gex)))
clin_all<-rbind(data_sub1_all1,data_sub1_all2,data_sub1_all3,data_sub1_all4,data_sub1_all5)

p_val<-matrix(nrow=1, ncol=2)

#subtype 1(Basal-like)
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_sub1$sample_id)+1-rev(which(data_sub1$vital_status==1))
deaths<-(c(0,rev(which(data_sub1$vital_status==1)))-c(rev(which(data_sub1$vital_status==1)),0))[-c(1)]

data_sub1$Number_alive<-c(rep(number_alive[1],length(data_sub1$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_sub1$Deaths<-c(rep(0,length(data_sub1$sample_id)))
data_sub1$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_sub1$prop_surv<-rev((data_sub1$Number_alive-data_sub1$Deaths)/data_sub1$Number_alive)
data_sub1$prop_surv[is.na(data_sub1$prop_surv)]<-1
data_sub1$cum_prop_surv<-cumprod(data_sub1$prop_surv)
#Median survival time
if(is.na(which(data_sub1$cum_prop_surv<=0.75)[1])){
  Median[1]<-c("No median")
  Quantile[1]<-c("No Quantile")}else
    if(is.na(which(data_sub1$cum_prop_surv<=0.5)[1])){
      Median[1]<-c("No median")
      Quantile[1]<-data_sub1$survival_time[which(data_sub1$cum_prop_surv<=0.75)[1]]
    }else
      if(data_sub1$vital_status[which(data_sub1$cum_prop_surv<=0.5)[1]]==1){
        Median[1]<-data_sub1$survival_time[which(data_sub1$cum_prop_surv<=0.5)[1]]
        Quantile[1]<-c("Don't need")}else 
          if(data_sub1$vital_status[which(data_sub1$cum_prop_surv<=0.75)[1]]==1){
            Quantile[1]<-data_sub1$survival_time[which(data_sub1$cum_prop_surv<=0.75)[1]]}else{
              Median[1]<-c("No median")
              Quantile[1]<-c("No Quantile")}
###########################

#subtype 2(HER2-enriched)
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_sub2$sample_id)+1-rev(which(data_sub2$vital_status==1))
deaths<-(c(0,rev(which(data_sub2$vital_status==1)))-c(rev(which(data_sub2$vital_status==1)),0))[-c(1)]

data_sub2$Number_alive<-c(rep(number_alive[1],length(data_sub2$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_sub2$Deaths<-c(rep(0,length(data_sub2$sample_id)))
data_sub2$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_sub2$prop_surv<-rev((data_sub2$Number_alive-data_sub2$Deaths)/data_sub2$Number_alive)
data_sub2$prop_surv[is.na(data_sub2$prop_surv)]<-1
data_sub2$cum_prop_surv<-cumprod(data_sub2$prop_surv)
#Median survival time
if(is.na(which(data_sub2$cum_prop_surv<=0.75)[1])){
  Median[2]<-c("No median")
  Quantile[2]<-c("No Quantile")}else
    if(is.na(which(data_sub2$cum_prop_surv<=0.5)[1])){
      Median[2]<-c("No median")
      Quantile[2]<-data_sub2$survival_time[which(data_sub2$cum_prop_surv<=0.75)[1]]
    }else
      if(data_sub2$vital_status[which(data_sub2$cum_prop_surv<=0.5)[1]]==1){
        Median[2]<-data_sub2$survival_time[which(data_sub2$cum_prop_surv<=0.5)[1]]
        Quantile[2]<-c("Don't need")}else 
          if(data_sub2$vital_status[which(data_sub2$cum_prop_surv<=0.75)[1]]==1){
            Quantile[2]<-data_sub2$survival_time[which(data_sub2$cum_prop_surv<=0.75)[1]]}else{
              Median[2]<-c("No median")
              Quantile[2]<-c("No Quantile")}
###########################

#subtype 3(Luminal A)
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_sub3$sample_id)+1-rev(which(data_sub3$vital_status==1))
deaths<-(c(0,rev(which(data_sub3$vital_status==1)))-c(rev(which(data_sub3$vital_status==1)),0))[-c(1)]

data_sub3$Number_alive<-c(rep(number_alive[1],length(data_sub3$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_sub3$Deaths<-c(rep(0,length(data_sub3$sample_id)))
data_sub3$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_sub3$prop_surv<-rev((data_sub3$Number_alive-data_sub3$Deaths)/data_sub3$Number_alive)
data_sub3$prop_surv[is.na(data_sub3$prop_surv)]<-1
data_sub3$cum_prop_surv<-cumprod(data_sub3$prop_surv)
#Median survival time
if(is.na(which(data_sub3$cum_prop_surv<=0.75)[1])){
  Median[3]<-c("No median")
  Quantile[3]<-c("No Quantile")}else
    if(is.na(which(data_sub3$cum_prop_surv<=0.5)[1])){
      Median[3]<-c("No median")
      Quantile[3]<-data_sub3$survival_time[which(data_sub3$cum_prop_surv<=0.75)[1]]
    }else
      if(data_sub3$vital_status[which(data_sub3$cum_prop_surv<=0.5)[1]]==1){
        Median[3]<-data_sub3$survival_time[which(data_sub3$cum_prop_surv<=0.5)[1]]
        Quantile[3]<-c("Don't need")}else 
          if(data_sub3$vital_status[which(data_sub3$cum_prop_surv<=0.75)[1]]==1){
            Quantile[3]<-data_sub3$survival_time[which(data_sub3$cum_prop_surv<=0.75)[1]]}else{
              Median[3]<-c("No median")
              Quantile[3]<-c("No Quantile")}
###########################


#subtype 4(Luminal B)
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_sub4$sample_id)+1-rev(which(data_sub4$vital_status==1))
deaths<-(c(0,rev(which(data_sub4$vital_status==1)))-c(rev(which(data_sub4$vital_status==1)),0))[-c(1)]

data_sub4$Number_alive<-c(rep(number_alive[1],length(data_sub4$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_sub4$Deaths<-c(rep(0,length(data_sub4$sample_id)))
data_sub4$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_sub4$prop_surv<-rev((data_sub4$Number_alive-data_sub4$Deaths)/data_sub4$Number_alive)
data_sub4$prop_surv[is.na(data_sub4$prop_surv)]<-1
data_sub4$cum_prop_surv<-cumprod(data_sub4$prop_surv)
#Median survival time
if(is.na(which(data_sub4$cum_prop_surv<=0.75)[1])){
  Median[4]<-c("No median")
  Quantile[4]<-c("No Quantile")}else
    if(is.na(which(data_sub4$cum_prop_surv<=0.5)[1])){
      Median[4]<-c("No median")
      Quantile[4]<-data_sub4$survival_time[which(data_sub4$cum_prop_surv<=0.75)[1]]
    }else
      if(data_sub4$vital_status[which(data_sub4$cum_prop_surv<=0.5)[1]]==1){
        Median[4]<-data_sub4$survival_time[which(data_sub4$cum_prop_surv<=0.5)[1]]
        Quantile[4]<-c("Don't need")}else 
          if(data_sub4$vital_status[which(data_sub4$cum_prop_surv<=0.75)[1]]==1){
            Quantile[4]<-data_sub4$survival_time[which(data_sub4$cum_prop_surv<=0.75)[1]]}else{
              Median[4]<-c("No median")
              Quantile[4]<-c("No Quantile")}
###########################


#subtype 5(Normal-like)
#set r,d for Kaplan-Meier survival curve
number_alive<-length(data_sub5$sample_id)+1-rev(which(data_sub5$vital_status==1))
deaths<-(c(0,rev(which(data_sub5$vital_status==1)))-c(rev(which(data_sub5$vital_status==1)),0))[-c(1)]

data_sub5$Number_alive<-c(rep(number_alive[1],length(data_sub5$sample_id)-sum(deaths)),rep(number_alive,deaths))
data_sub5$Deaths<-c(rep(0,length(data_sub5$sample_id)))
data_sub5$Deaths[number_alive[]]<-c(1)
#get p=(r-d)/r
data_sub5$prop_surv<-rev((data_sub5$Number_alive-data_sub5$Deaths)/data_sub5$Number_alive)
data_sub5$prop_surv[is.na(data_sub5$prop_surv)]<-1
data_sub5$cum_prop_surv<-cumprod(data_sub5$prop_surv)
#Median survival time
if(is.na(which(data_sub5$cum_prop_surv<=0.75)[1])){
  Median[5]<-c("No median")
  Quantile[5]<-c("No Quantile")}else
    if(is.na(which(data_sub5$cum_prop_surv<=0.5)[1])){
      Median[5]<-c("No median")
      Quantile[5]<-data_sub5$survival_time[which(data_sub5$cum_prop_surv<=0.75)[1]]
    }else
      if(data_sub5$vital_status[which(data_sub5$cum_prop_surv<=0.5)[1]]==1){
        Median[5]<-data_sub5$survival_time[which(data_sub5$cum_prop_surv<=0.5)[1]]
        Quantile[5]<-c("Don't need")}else 
          if(data_sub5$vital_status[which(data_sub5$cum_prop_surv<=0.75)[1]]==1){
            Quantile[5]<-data_sub5$survival_time[which(data_sub5$cum_prop_surv<=0.75)[1]]}else{
              Median[5]<-c("No median")
              Quantile[5]<-c("No Quantile")}
###########################
Median
Quantile

#graph
library(survival)
library(survminer)
ev<-1*(clin$vital_status==1)
fut<-as.numeric(clin$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin$subtype,data=clin)
ggsurvplot(fit,data=clin,pval=TRUE)
fit<-survfit(su~clin$subtype,data=clin)
#find p-value
survdiff(su~clin$subtype,data=clin)
sdf<-survdiff(su~clin$subtype,data=clin)
p_val[1]<-1-pchisq(sdf$chisq, length(sdf$n)-1)

#With all data
ev<-1*(clin_all$vital_status==1)
fut<-as.numeric(clin_all$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin_all$subtype,data=clin_all)
ggsurvplot(fit,data=clin_all,pval=TRUE)
fit<-survfit(su~clin_all$subtype,data=clin_all)
#find p-value
survdiff(su~clin_all$subtype,data=clin_all)
sdf<-survdiff(su~clin_all$subtype,data=clin_all)
p_val[2]<-1-pchisq(sdf$chisq, length(sdf$n)-1)

