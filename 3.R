library(ClusterR)
library(cluster)
library(fpc)
library(factoextra)
library(magrittr)
library(ggplot2)
library(tibble)
library(plyr)
library(dplyr)
library(clusteval)
library(NbClust)
#load data
gex<-readRDS("expression.rds")
clin<-readRDS("clinical.rds")

#Demension reduction of gex data
PCA_DATA <- t(gex)
PCA <- prcomp(na.omit(PCA_DATA))
plot(PCA$x[,1]~PCA$x[,2])
plot(PCA$x[,1]~PCA$x[,3])
REDUCED <- PCA$x[,1:3]

#K-means clustering
kmeans_result<-kmeans(REDUCED,3)
table(kmeans_result$cluster)
par(mfrow=c(1,1))
plot(REDUCED[,c(1,2)], col=kmeans_result$cluster, main = "PC1 & PC2")
plot(REDUCED[,c(2,3)], col=kmeans_result$cluster, main = "PC2 & PC3")
plot(REDUCED[,c(1,3)], col=kmeans_result$cluster, main = "PC1 & PC3")


#Hierarchical clustering
hc_complete<-hclust(dist(REDUCED), method="complete")
plot(hc_complete, main="Complete Linkage", xlab="", sub="", cex=.9)
rect.hclust(hc_complete, k=3)

hc_average<-hclust(dist(REDUCED), method="average")
plot(hc_average, main="Average Linkage", xlab="", sub="", cex=.9)
rect.hclust(hc_average, k=3)

hc_single<-hclust(dist(REDUCED), method="single")
plot(hc_single, main="Single Linkage", xlab="", sub="", cex=.9)
rect.hclust(hc_single, k=3)

# Optibal cluster number
fviz_nbclust(REDUCED, kmeans, method="silhouette")
fviz_nbclust(REDUCED, kmeans, method="wss")
fviz_nbclust(REDUCED, kmeans, method="gap_stat")

# Confirm whether dimension reduction to 3 is right
REDUCED_e <- PCA$x[,1:5]
fviz_nbclust(REDUCED_e, kmeans, method="silhouette")

#Scree plot
REDUCED_dist<- dist(REDUCED, method="euclidean")
screeplot(PCA,type="line", ylab="variance", npcs = 20, main = "Result of PCA")

#Gap statistic
gap_stat<- clusGap(REDUCED, FUN=kmeans, nstart=20, K.max = 10)
fviz_gap_stat(gap_stat)


#Make merged data
REDUCED_mer <- data.frame(as.character(rownames(REDUCED)), REDUCED)
colnames(REDUCED_mer)[1] <- "sample_id"
rownames(REDUCED_mer) <- NULL

clin_cluster_sample <- merge(clin[,c("sample_id","survival_index")], REDUCED_mer, by = "sample_id", all = FALSE)
clin_cluster_sample


sample1<- clin_cluster_sample
sample1$survival_index <- NULL
sample1$sample_id <- NULL

data_km <- kmeans(sample1, 4)$cluster
data_hc <- cutree(hclust(dist(sample1), method="complete"), 4)
print(clin_cluster_sample$survival_index)
print(data_km)

print(cluster_similarity(clin_cluster_sample$survival_index, data_km))
print(cluster_similarity(clin_cluster_sample$survival_index, data_hc))



###With best feature set
feature_set1<-c("KDELR1",  "LIPG",    "PER1", "CAPRIN2", "MOBKL2B", "KRT20",  "C1RL", "SORD", "CDKN1C")
PCA_DATA <- t(gex)
number1<-matrix(nrow = 1, ncol = length(feature_set1))
for (i in 1:length(colnames(PCA_DATA))){
  for (j in 1:length(feature_set1)){
    if (colnames(PCA_DATA)[i]==feature_set1[j]){
      number1[j]<-i
    }
  }
}
combined_data1<-rbind(PCA_DATA[,number1[1,]])

PCA_cb1 <- prcomp(na.omit(combined_data1))
REDUCED_cb1 <- PCA_cb$x[,1:3]

REDUCED_mer_cb1 <- data.frame(as.character(rownames(REDUCED_cb1)), REDUCED_cb1)
colnames(REDUCED_mer_cb1)[1] <- "sample_id"
rownames(REDUCED_mer_cb1) <- NULL

clin_cluster_sample_cb1 <- merge(clin[,c("sample_id","survival_index")], REDUCED_mer_cb1, by = "sample_id", all = FALSE)

sample_cb1<- clin_cluster_sample_cb1
sample_cb1$survival_index <- NULL
sample_cb1$sample_id <- NULL

data_km1 <- kmeans(sample_cb1, 4)$cluster
data_hc1 <- cutree(hclust(dist(sample_cb1), method="complete"), 4)

print(cluster_similarity(clin_cluster_sample_cb1$survival_index, data_km1))
print(cluster_similarity(clin_cluster_sample_cb1$survival_index, data_hc1))



feature_set2<-c( "IDE", "ZDHHC20", "DARC", "SNF1LK", "ACOX2", "CLEC12A", "HUWE1",  
                 "BAALC", "ZNF284", "PNCK", "RALGPS2", "ZBTB37")
number2<-matrix(nrow = 1, ncol = length(feature_set2))
for (i in 1:length(colnames(PCA_DATA))){
  for (j in 1:length(feature_set2)){
    if (colnames(PCA_DATA)[i]==feature_set2[j]){
      number2[j]<-i
    }
  }
}
combined_data2<-rbind(PCA_DATA[,number2[1,]])

PCA_cb2 <- prcomp(na.omit(combined_data2))
REDUCED_cb2 <- PCA_cb2$x[,1:3]

REDUCED_mer_cb2 <- data.frame(as.character(rownames(REDUCED_cb2)), REDUCED_cb2)
colnames(REDUCED_mer_cb2)[1] <- "sample_id"
rownames(REDUCED_mer_cb2) <- NULL

clin_cluster_sample_cb2 <- merge(clin[,c("sample_id","survival_index")], REDUCED_mer_cb2, by = "sample_id", all = FALSE)


sample_cb2<- clin_cluster_sample_cb2
sample_cb2$survival_index <- NULL
sample_cb2$sample_id <- NULL

data_km2 <- kmeans(sample_cb2, 4)$cluster
data_hc2 <- cutree(hclust(dist(sample_cb2), method="complete"), 4)

print(cluster_similarity(clin_cluster_sample_cb2$survival_index, data_km2))
print(cluster_similarity(clin_cluster_sample_cb2$survival_index, data_hc2))
