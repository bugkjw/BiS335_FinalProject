#install.packages("e1071")
library(e1071)
library(MASS)
#install.packages("randomForest")
library(randomForest)
#install.packages("gbm")
library(gbm)
#install.packages("caret")
library(caret)

setwd("C:/Users/Chan/Desktop/Finalterm-Project")

# Data import
clin <- readRDS("clinical.rds");
gex <- readRDS("expression.rds");
mut <- readRDS("mutation.rds");
gex_res <- read.csv("gex_anova_result.csv");
mut_res <- read.csv("mut_chisq_result.csv");

# Use the class label found in #1, generate entire labeled dataset
labels <- unique(clin$survival_index);
clin_label <- data.frame(clin$sample_id,clin$survival_index); colnames(clin_label) <- c("sample_id","survival_index")

# Dataset separation: 80% development set, 20% hold-out set (for testing the final model)
shuffle <- clin_label[sample(nrow(clin_label)),]; rownames(shuffle) <- NULL
SELECT <- sample(nrow(shuffle)*0.2)
DEV_DATA <- shuffle[-SELECT,]
HOLDOUT_DATA <- shuffle[SELECT,]

# Filter gene expression and mutation data based on dataset separation
gex_DEV <- gex[,colnames(gex) %in% DEV_DATA$sample_id]
gex_HOLDOUT <- gex[,colnames(gex) %in% HOLDOUT_DATA$sample_id]

# Candidate genomic features to be used as features while fitting
gex_feature_all <- as.character(gex_res$Hugo_Symbol);
gex_feature_all <- gex_feature_all[gex_feature_all %in% rownames(gex_DEV)]; factor(gex_feature_all)
feature_all <- gex_feature_all

# Feature selection by forward stepwise selection
# For each feature set expansion, assess 5-fold CV misclassification error rate
# Expand until no more error rate reduction is possible

# Specifically,
# 1. Try each possible features that are not in the current set
# 2. Get the feature with maximum error (5-fold CV) reduction
# 3. Add to the current feature set
# 4. Repeat until possible error deviation is all >0
gex_feature <- c();
current_feature <- c()
CV_errors_procedure <- c()
last_min_error <- 1;
for(ii in 1:length(feature_all)){
  # Candidate features: Ones not in the current feature set
  candidates <- feature_all[!(feature_all %in% current_feature)]
  candidate_errors <- c()
  
  for (jj in 1:length(candidates)){
    tmp_gex_feature <- c();
    # Temporary expanded feature set
    test_feature <- c(current_feature, candidates[jj])
    # Temporaty expanded gex and mut feature set
    if (candidates[jj] %in% gex_feature_all){
      tmp_gex_feature <- c(gex_feature, candidates[jj])
    }
    # Dataset condtruction
    # With the current feature set (to be tested),
    # make a table containing class labels and feature values of each patient.
    cat(sprintf("Feature %d selection, looking at candidate %d\n", ii, jj))
    
    if (length(tmp_gex_feature) == 0){cat("No gex feature tested. Continue...\n")
      flag <- 1
    }
    
    else{
      flag <- 0
      if (length(tmp_gex_feature) == 1){
        gex_dataset <- data.frame(gex_DEV[tmp_gex_feature,])
        colnames(gex_dataset) <- tmp_gex_feature
      }
      else{
        gex_dataset <- t(gex_DEV[tmp_gex_feature,])
      }
      gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
      clin_dataset <- merge(DEV_DATA, gex_dataset, by = "sample_id", all = FALSE)
    }
    
      if (flag == 1){clin_dataset <- merge(DEV_DATA, mut_dataset, by = "sample_id", all = FALSE)
      }
    else if (flag == 0){clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)}
    }
    # Now we have the feature set and whole dataset.
    # Split the dataset for cross validation
    tmp_folds <- cut(seq(1,nrow(clin_dataset)),breaks=5,labels=FALSE)
    
    # 5-fold CV
    CV_error <- c()
    for (kk in 1:5){
      tmp_testIndexes <- which(tmp_folds==kk,arr.ind=TRUE)
      testData <- clin_dataset[tmp_testIndexes, ]
      trainData <- clin_dataset[-tmp_testIndexes, ]
      
      # Support vector machine
      set.seed(1); cost <- c(0.001, 0.01, 0.1, 1, 5, 10, 100)
      tune.out <- tune(svm,survival_index~.-sample_id, data = trainData, kernel = "linear",
                       ranges = list(cost))
      best_cost <- tune.out$best.parameters
      svm.trainData <- svm(survival_index ~ .-sample_id,
                           data = trainData, kernel = "linear", cost = best_cost[1,1], scale = FALSE)
      # Test
      svm.pred <- predict(svm.trainData, testData, type = "class")
      cMat <- confusionMatrix(svm.pred,testData$survival_index[as.integer(names(svm.pred))-min(as.integer(names(svm.pred)))+1])
      CV_error[kk] <- 1-as.numeric(cMat$overall[1])
    }
    candidate_errors[jj] <- mean(CV_error)
  }
  CV_errors_procedure[ii] <- min(candidate_errors)
  # Check for improvement
  if (min(candidate_errors) > last_min_error){break}
  # Update threshold
  last_min_error <- min(candidate_errors)
  # Greedy feature selection
  new_feature <- candidates[candidate_errors == min(candidate_errors)]
  new_feature <- new_feature[1]
  # Expanded feature set
  current_feature <- c(current_feature, new_feature)
  # Expanded gex and mut feature set
  if (new_feature %in% gex_feature_all){
    gex_feature <- c(gex_feature, new_feature)
  }
}
png(filename="./2/SVM/CV_errors.png")
{par(mfrow = c(1,1))
  plot(1:ii, CV_errors_procedure, xlab = "# of features", ylab = "5-fold CV error", xlim = c(0,ii+1))
  lines(1:ii, CV_errors_procedure)
  text(1:ii, CV_errors_procedure, labels = c(current_feature,"(Terminate)"))
}
dev.off()

# Test error estimate using hold-out set
# Dataset condtruction
# With the feature set, make a table containing class labels and feature values of each patient.
# Full data (for fitting)
  cat("Full dataset construction\n")
  if (length(gex_feature) == 0){cat("No gex feature tested. Continue...\n")
    flag <- 1
  }
  
  else{
    flag <- 0
    if (length(gex_feature) == 1){
      gex_dataset <- data.frame(gex[gex_feature,])
      colnames(gex_dataset) <- gex_feature
    }
    else{
      gex_dataset <- t(gex[gex_feature,])
    }
    gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
    clin_dataset_full <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
  }

  # Hold-out
  cat("Hold-out dataset construction\n")
  if (length(gex_feature) == 0){
    cat("No gex feature tested. Continue...\n")
    flag <- 1
  }
  
  else{
    flag <- 0
    if (length(gex_feature) == 1){
      gex_dataset <- data.frame(gex_HOLDOUT[gex_feature,])
      colnames(gex_dataset) <- gex_feature
    }
    
    else{
      gex_dataset <- t(gex_HOLDOUT[gex_feature,])
    }
    
    gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
    clin_dataset <- merge(HOLDOUT_DATA, gex_dataset, by = "sample_id", all = FALSE)
  }

# ------------------------------------------------------------- Here, we should use the gex excel file
# Dataset ready
library(ClusterR)
library(cluster)
library(fpc)
library(factoextra)
library(magrittr)
library(ggplot2)
library(tibble)
library(plyr)

clin_index <- clin_dataset_full[,2]
  
# Plot basic K-means clustering
clin_dataset_full <- na.omit(clin_dataset_full)
clin_dataset_full$sample_id <- NULL
clin_dataset_full$survival_index <- NULL
km.result <- kmeans(clin_dataset_full, 5, nstart=20)
print(km.result)
plotcluster(clin_dataset_full, km.result$cluster)

# Silhouette plot
fviz_nbclust(clin_dataset_full, kmeans, method="silhouette") # 3이나 5가 좋다

# Scree plot
# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
wss <- (nrow(clin_dataset_full)-1)*sum(apply(clin_dataset_full, 2, var))
for(i in 1:15){
  wss[i] <- sum(kmeans(clin_dataset_full, centers=i)$withinss)
}

plot(1:15, wss, type="b", xlab="Numbers of Clusters", ylab="Within groups sum of squares") # 4~5에서 drop 이 있으므로 5가 좋다?
# We might can conclude in scree plot that we need 5 clusters

# -------------------------------------------------------------
# Basic Hierarchical clustering
hc_complete <- hclust(dist(clin_dataset_full), method="complete")
hc_average <- hclust(dist(clin_dataset_full), method="average")
hc_single <- hclust(dist(clin_dataset_full), method="single")

par(mfrow=c(1,3))
plot(hc_complete, main="Hierarchical ~ Complete linkage", xlab="", sub="", cex=.9)
plot(hc_average, main="Hierarchical ~ Average linkage", xlab="", sub="", cex=.9)
plot(hc_single, main="Hierarchical ~ Single linkage", xlab="", sub="", cex=.9)

# Cutree for data output
print(cutree(hc_complete, 4))
print(cutree(hc_average, 4))
print(cutree(hc_single, 4))

# Silhouette plot
fviz_nbclust(clin_dataset_full, FUN=hcut, method="silhouette") # 3이나 5가 좋ㄷ

# Scree plot
clin_dist <- dist(clin_dataset_full, method="euclidean")
ggplot(hc_complete$height %>% as.tibble() %>% add_column(groups=length(hc_complete$height):1) %>% rename(height=value), aes(x=groups, y=height))+geom_point()+geom_line()


# Hierarchical clustering with scaled features
clin_scaled <- scale(clin_dataset_full)
hc_complete_scaled <- hclust(dist(clin_scaled), method="complete")
hc_average_scaled <- hclust(dist(clin_scaled), method="average")
hc_single_scaled <- hclust(dist(clin_scaled), method="single")

par(mfrow=c(1,3))
plot(hclust(dist(clin_scaled), method="complete"), main="Hierarchical ~ Complete, Scaled", xlab="", sub="", cex=.9)
plot(hclust(dist(clin_scaled), method="average"), main="Hierarchical ~ Average, Scaled", xlab="", sub="", cex=.9)
plot(hclust(dist(clin_scaled), method="single"), main="Hierarchical ~ Single, Scaled", xlab="", sub="", cex=.9)

print(cutree(hc_complete_scaled, 4))
print(cutree(hc_average_scaled, 4))
print(cutree(hc_single_scaled, 4))

par(mfrow=c(1,1))
plot(hc_complete_scaled, hang <- -1, labels=cutree(hc_complete_scaled, k=4), xlab="", sub="", cex=.9)
rect.hclust(hc_complete_scaled, k=4, which=NULL, x=NULL, h=NULL, border=2, cluster=NULL)

# Silhouette plot
fviz_nbclust(clin_scaled, FUN=hcut, method="silhouette") # 3이 좋다

# Scree plot



# ------------------------------------------------------------- a, b complete
# https://www.statmethods.net/advstats/cluster.html
# Similarity test for K-means clustering
print(clin_index)
cluster.stats(clin_dataset_full, km.result$cluster, clin_index)

# Similarity test for Hierarchical clustering
print(cutree())
print(clin_index)
cluster.stats(clin_dataset_full, km.result$cluster, clin_index)

# Similarity test for scaled Hierarchical clustering
print(cutree())
print(clin_index)
cluster.stats(clin_dataset_full, km.result$cluster, clin_index)

# ------------------------------------------------------------- From here, we should use the result of #2
# clin_dataset_full2 <- data 집어넣기

# Plot basic K-means clustering
clin_dataset_full2 <- na.omit(clin_dataset_full2)
clin_dataset_full2$sample_id <- NULL
clin_dataset_full2$survival_index <- NULL
km.result2 <- kmeans(clin_dataset_full2, 5, nstart=20)
print(km.result2)
plotcluster(clin_dataset_full2, km.result2$cluster)

# Silhouette plot
fviz_nbclust(clin_dataset_full2, kmeans, method="silhouette")

# Scree plot
wss <- (nrow(clin_dataset_full2)-1)*sum(apply(clin_dataset_full2, 2, var))
for(i in 1:15){
  wss[i] <- sum(kmeans(clin_dataset_full2, centers=i)$withinss)
}

plot(1:15, wss, type="b", xlab="Numbers of Clusters", ylab="Within groups sum of squares")

# -------------------------------------------------------------
# Basic Hierarchical clustering
hc_complete2 <- hclust(dist(clin_dataset_full2), method="complete")
hc_average2 <- hclust(dist(clin_dataset_full2), method="average")
hc_single2 <- hclust(dist(clin_dataset_full2), method="single")

par(mfrow=c(1,3))
plot(hc_complete2, main="Hierarchical ~ Complete linkage", xlab="", sub="", cex=.9)
plot(hc_average2, main="Hierarchical ~ Average linkage", xlab="", sub="", cex=.9)
plot(hc_single2, main="Hierarchical ~ Single linkage", xlab="", sub="", cex=.9)

# Cutree for data output
print(cutree(hc_complete2, 4))
print(cutree(hc_average2, 4))
print(cutree(hc_single2, 4))

# Silhouette plot
fviz_nbclust(clin_dataset_full2, FUN=hcut, method="silhouette")

# Scree plot


# Hierarchical clustering with scaled features
clin_scaled2 <- scale(clin_dataset_full2)
hc_complete_scaled2 <- hclust(dist(clin_scaled2), method="complete")
hc_average_scaled2 <- hclust(dist(clin_scaled2), method="average")
hc_single_scaled2 <- hclust(dist(clin_scaled2), method="single")

par(mfrow=c(1,3))
plot(hclust(dist(clin_scaled2), method="complete"), main="Hierarchical ~ Complete, Scaled", xlab="", sub="", cex=.9)
plot(hclust(dist(clin_scaled2), method="average"), main="Hierarchical ~ Average, Scaled", xlab="", sub="", cex=.9)
plot(hclust(dist(clin_scaled2), method="single"), main="Hierarchical ~ Single, Scaled", xlab="", sub="", cex=.9)

print(cutree(hc_complete_scaled2, 4))
print(cutree(hc_average_scaled2, 4))
print(cutree(hc_single_scaled2, 4))

par(mfrow=c(1,1))
plot(hc_complete_scaled2, hang <- -1, labels=cutree(hc_complete_scaled2, k=4), xlab="", sub="", cex=.9)
rect.hclust(hc_complete_scaled2, k=4, which=NULL, x=NULL, h=NULL, border=2, cluster=NULL)

# Silhouette plot
fviz_nbclust(clin_scaled2, FUN=hcut, method="silhouette")

# Scree plot

# ------------------------------------------------------------- Similarity test for #2 clustering
# Similarity test for K-means clustering
print(clin_index)
cluster.stats(clin_dataset_full, km.result$cluster, clin_index)

# Similarity test for Hierarchical clustering
print(cutree())
print(clin_index)
cluster.stats(clin_dataset_full, km.result$cluster, clin_index)

# Similarity test for scaled Hierarchical clustering
print(cutree())
print(clin_index)
cluster.stats(clin_dataset_full, km.result$cluster, clin_index)
