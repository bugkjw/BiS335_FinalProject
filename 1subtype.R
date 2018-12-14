clin <- readRDS("C:/Users/Shin Dongho/Desktop/Finalterm-Project/clinical.Rds")
mut <- readRDS("C:/Users/Shin Dongho/Desktop/Finalterm-Project/mutation.Rds")
expre <- readRDS("C:/Users/Shin Dongho/Desktop/Finalterm-Project/expression.Rds")

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

#subtype 1
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
#median survival time
data_sub1$survival_time[which(data_sub1$cum_prop_surv<=0.5&data_sub1$cum_prop_surv!=0)[1]]
###########################

#subtype 2
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
#median survival time
data_sub2$survival_time[which(data_sub2$cum_prop_surv<=0.5&data_sub2$cum_prop_surv!=0)[1]]
###########################

#subtype 3
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
#median survival time
data_sub3$survival_time[which(data_sub3$cum_prop_surv<=0.5&data_sub3$cum_prop_surv!=0)[1]]
c("quantile",data_sub3$survival_time[which(data_sub3$cum_prop_surv<=0.75&data_sub3$cum_prop_surv!=0)[1]])
###########################


#subtype 4
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
#median survival time
data_sub4$survival_time[which(data_sub4$cum_prop_surv<=0.5&data_sub4$cum_prop_surv!=0)[1]]
###########################


#subtype 5
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
#median survival time
data_sub5$survival_time[which(data_sub5$cum_prop_surv<=0.5&data_sub5$cum_prop_surv!=0)[1]]
c("quantile",data_sub5$survival_time[which(data_sub5$cum_prop_surv<=0.75&data_sub5$cum_prop_surv!=0)[1]])
###########################

#graph
library(survminer)
ev<-1*(clin$vital_status==1)
fut<-as.numeric(clin$survival_time)
su<-Surv(fut,ev)
fit<-survfit(su~clin$subtype,data=clin)
ggsurvplot(fit,data=clin,pval=TRUE)
fit<-survfit(su~clin$subtype,data=clin)
#find p-value
survdiff(su~clin$subtype,data=clin)