#install.packages("tree")
library(tree)
library(MASS)
#install.packages("randomForest")
library(randomForest)
#install.packages("gbm")
library(gbm)
#install.packages("caret")
library(caret)

setwd("C:/Users/VSlab#10/Desktop/JinwooKim/BiS335_FinalProject_Folder")
#setwd("D:/윈도우계정/Desktop/!/3학년3가을학기/BiS335 Biomedical Statistics & Statistical Learning/Final Project/Finalterm-Project")

# Data import
clin <- readRDS("./Data/clinical.rds");
gex <- readRDS("./Data/expression.rds");
gex_res <- read.csv("./Data/gex_anova_result.csv");

# Use the class label found in #1, generate entire labeled dataset
labels <- unique(clin$survival_index);
clin_label <- data.frame(clin$sample_id,clin$survival_index); colnames(clin_label) <- c("sample_id","survival_index")

# Dataset separation: 80% development set, 20% hold-out set (for testing the final model)
shuffle <- clin_label[sample(nrow(clin_label)),]; rownames(shuffle) <- NULL
SELECT <- sample(nrow(shuffle)*0.2)
DEV_DATA <- shuffle[-SELECT,]
HOLDOUT_DATA <- shuffle[SELECT,]

# Filter gene expression data based on dataset separation
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
gex_feature <- c()
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
    # Temporaty expanded gex feature set
    tmp_gex_feature <- c(gex_feature, candidates[jj])
    # Dataset condtruction
    # With the current feature set (to be tested),
    # make a table containing class labels and feature values of each patient.
    cat(sprintf("Feature %d selection, looking at candidate %d\n", ii, jj))
    if (length(tmp_gex_feature) == 1){
      gex_dataset <- data.frame(gex_DEV[tmp_gex_feature,])
      colnames(gex_dataset) <- tmp_gex_feature
    }else{
      gex_dataset <- t(gex_DEV[tmp_gex_feature,])
    }
    gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
    clin_dataset <- merge(DEV_DATA, gex_dataset, by = "sample_id", all = FALSE)
    # Now we have the feature set and whole dataset.
    # Split the dataset for cross validation
    tmp_folds <- cut(seq(1,nrow(clin_dataset)),breaks=5,labels=FALSE)
    
    # 5-fold CV
    CV_error <- c()
    for (kk in 1:5){
      tmp_testIndexes <- which(tmp_folds==kk,arr.ind=TRUE)
      testData <- clin_dataset[tmp_testIndexes, ]
      trainData <- clin_dataset[-tmp_testIndexes, ]
      
      # Linear Discriminant Analysis
      lda.trainData <- lda(survival_index ~ .-sample_id, data = trainData)
      lda.pred <- predict(lda.trainData, testData, type = "class")$class
      cMat <- confusionMatrix(lda.pred,testData$survival_index)
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
  # Expanded gex feature set
  gex_feature <- c(gex_feature, new_feature)
}
png(filename="./Result/LDA/CV_errors.png")
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
if (length(gex_feature) == 1){
  gex_dataset <- data.frame(gex_DEV[gex_feature,])
  colnames(gex_dataset) <- gex_feature
}else{
  gex_dataset <- t(gex_DEV[gex_feature,])
}
gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
clin_dataset_DEV <- merge(DEV_DATA, gex_dataset, by = "sample_id", all = FALSE)
# Hold-out
cat("Hold-out dataset construction\n")
if (length(gex_feature) == 1){
  gex_dataset <- data.frame(gex_HOLDOUT[gex_feature,])
  colnames(gex_dataset) <- gex_feature
}else{
  gex_dataset <- t(gex_HOLDOUT[gex_feature,])
}
gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
clin_dataset_HOLDOUT <- merge(HOLDOUT_DATA, gex_dataset, by = "sample_id", all = FALSE)

# Dataset ready
# Estimate test error using hold-out set
# Train by dev dataset, test by hold-out set!
# LDA model fitting
lda.Final <- lda(survival_index ~ .-sample_id, data = clin_dataset_DEV)
# Training error
lda.pred <- predict(lda.Final, clin_dataset_DEV, type = "class")
cMat <- confusionMatrix(lda.pred$class,clin_dataset_DEV$survival_index)
LDAError_T <- 1-as.numeric(cMat$overall[1])
{
  # Performance
  png(filename="./Result/LDA/LDA_training.png")
  plot(clin_dataset_DEV$survival_index,lda.pred$class,xlab = "Training data", ylab = "Prediction by LDA", col= c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}
# Test
lda.pred <- predict(lda.Final, clin_dataset_HOLDOUT, type = "class")
cMat <- confusionMatrix(lda.pred$class,clin_dataset_HOLDOUT$survival_index)
LDAError_F <- 1-as.numeric(cMat$overall[1])
{
  cat(sprintf("\nLDA performance on test set: %2.3g\n",cMat$overall[1]))
  cat(sprintf("Feature set: "))
  cat(toString(current_feature))
  cat("\n")
}
{
  # Performance
  png(filename="./Result/LDA/LDA_performance.png")
  plot(clin_dataset$survival_index,lda.pred$class,xlab = "Test data", ylab = "Prediction by LDA", col= c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}

# Visualization
#install.packages("klaR")
library(klaR)
{
  par(mfrow = c(1,1))
  # Test data and feature assessment
  lda.Final.values <- predict(lda.Final)
  # Mapping of training data onto a pair of discriminant variables
  # Data label
  png(filename="./Result/LDA/LDA_training_data_12.png")
  plot(lda.Final.values$x[,1],lda.Final.values$x[,2], xlab = "LD1", ylab = "LD2", col=clin_dataset_full$survival_index, pch = 16) # make a scatterplot
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("LDA training data distribution, LD1-LD2")
  dev.off()
  png(filename="./Result/LDA/LDA_training_data_23.png")
  plot(lda.Final.values$x[,2],lda.Final.values$x[,3], xlab = "LD2", ylab = "LD3", col=clin_dataset_full$survival_index, pch = 16) # make a scatterplot
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("LDA training data distribution, LD2-LD3")
  dev.off()
  png(filename="./Result/LDA/LDA_training_data_31.png")
  plot(lda.Final.values$x[,3],lda.Final.values$x[,1], xlab = "LD3", ylab = "LD1", col=clin_dataset_full$survival_index, pch = 16) # make a scatterplot
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("LDA training data distribution, LD3-LD1")
  dev.off()
}
{
  # Fitted model
  png(filename="./Result/LDA/LDA_training_fit_12.png")
  plot(lda.Final.values$x[,1],lda.Final.values$x[,2], xlab = "LD1", ylab = "LD2", col=lda.Final.values$class, pch = 16) # make a scatterplot
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("Fitted LDA, LD1-LD2")
  dev.off()
  png(filename="./Result/LDA/LDA_training_fit_23.png")
  plot(lda.Final.values$x[,2],lda.Final.values$x[,3], xlab = "LD2", ylab = "LD3", col=lda.Final.values$class, pch = 16) # make a scatterplot
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("Fitted LDA, LD2-LD3")
  dev.off()
  png(filename="./Result/LDA/LDA_training_fit_31.png")
  plot(lda.Final.values$x[,3],lda.Final.values$x[,1], xlab = "LD3", ylab = "LD1", col=lda.Final.values$class, pch = 16) # make a scatterplot
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("Fitted LDA, LD3-LD1")
  dev.off()
}
{
  # Test data performance
  # Data label
  png(filename="./Result/LDA/LDA_test_data_12.png")
  plot(lda.pred$x[,1],lda.pred$x[,2], xlab = "LD1", ylab = "LD2", col=clin_dataset$survival_index, pch = 16) # make a scatterplot
  legend("topright", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("LDA test data distribution, LD1-LD2")
  dev.off()
  png(filename="./Result/LDA/LDA_test_data_23.png")
  plot(lda.pred$x[,2],lda.pred$x[,3], xlab = "LD2", ylab = "LD3", col=clin_dataset$survival_index, pch = 16) # make a scatterplot
  legend("topright", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("LDA test data distribution, LD2-LD3")
  dev.off()
  png(filename="./Result/LDA/LDA_test_data_31.png")
  plot(lda.pred$x[,3],lda.pred$x[,1], xlab = "LD3", ylab = "LD1", col=clin_dataset$survival_index, pch = 16) # make a scatterplot
  legend("topright", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("LDA test data distribution, LD3-LD1")
  dev.off()
}
{ 
  # Prediction
  png(filename="./Result/LDA/LDA_test_fit_12.png")
  plot(lda.pred$x[,1],lda.pred$x[,2], xlab = "LD1", ylab = "LD2", col=lda.pred$class, pch = 16) # make a scatterplot
  legend("topright", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("Prediction by LDA, LD1-LD2")
  dev.off()
  png(filename="./Result/LDA/LDA_test_fit_23.png")
  plot(lda.pred$x[,2],lda.pred$x[,3], xlab = "LD1", ylab = "LD2", col=lda.pred$class, pch = 16) # make a scatterplot
  legend("topright", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("Prediction by LDA, LD2-LD3")
  dev.off()
  png(filename="./Result/LDA/LDA_test_fit_31.png")
  plot(lda.pred$x[,3],lda.pred$x[,1], xlab = "LD1", ylab = "LD2", col=lda.pred$class, pch = 16) # make a scatterplot
  legend("topright", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title("Prediction by LDA, LD3-LD1")
  dev.off()
}

save.image("./Result/LDA/v20181204_LDA_data.RData")
