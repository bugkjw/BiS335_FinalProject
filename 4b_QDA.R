library(MASS)
#install.packages("caret")
library(caret)
#install.packages("klaR")
library(klaR)
#install.packages("tree")
library(tree)
#install.packages("randomForest")
library(randomForest)
library(e1071)
library(gbm)
#install.packages("doParallel")
library(doParallel)
registerDoParallel(4)
getDoParWorkers()

#setwd("C:/Users/VSlab#10/Desktop/JinwooKim/BiS335_FinalProject_Folder")
setwd("D:/윈도우계정/Desktop/!/3학년3가을학기/BiS335 Biomedical Statistics & Statistical Learning/Final Project")

# Data import
clin <- readRDS("./Data/clinical.rds");
gex <- readRDS("./Data/expression.rds");
mut <- readRDS("./Data/mutation.rds");
gex_res <- read.csv("./Data/gex_anova_result.csv");
mut_res <- read.csv("./Data/mut_chisq_result.csv");

# Features selected in midterm project
TP1_gex <- c("PSMD3","C12orf28","FBXL20","SOHLH1","CXCL14","FLJ35767","TH","FVT1","ZC3H12D","hCH_2028557","CYP1B1","HS3ST2"
             ,"KCNJ5","RPL36A","PTGES3","FAM84B","RASGEF1A","LYPLA1","C20orf11","BAHCC1","HES6","IYD","LGR4","ZNF544","STAR","MAOB","ZNF498")
TP1_mut <- c("HMCN1","FAM5C","NUMA1","CYP24A1","FOXN1","C20orf112","NLGN4X","ITSN2","DOCK4","NCKAP1L")

# Load model data
load("./Data/LDA_data.Rdata")
LDA.model <- lda.Final
LDA_gex <- gex_feature

load("./Data/QDA_data.Rdata")
QDA.model <- qda.Final
QDA_gex <- gex_feature

load("./Data/Forest_data.Rdata")
Boost.model <- boost.Final
Boost.varimp <- summary(Boost.model)
Boost.varimp[TP1_gex[TP1_gex %in% Boost.varimp$var],]
Boost.varimp[TP1_mut[TP1_mut %in% Boost.varimp$var],]

# Use the class label found in #1, generate entire labeled dataset
labels <- unique(clin$survival_index);
clin_label <- data.frame(clin$sample_id,clin$survival_index); colnames(clin_label) <- c("sample_id","survival_index")

# Dataset separation: 80% development set, 20% hold-out set (for testing the final model)
shuffle <- clin_label[sample(nrow(clin_label)),]; rownames(shuffle) <- NULL
SELECT <- sample(nrow(shuffle)*0.2)
DEV_DATA <- shuffle[-SELECT,]
HOLDOUT_DATA <- shuffle[SELECT,]

# QDA: deficient data in class 2 yields errors. Remove
DEV_DATA <- DEV_DATA[DEV_DATA$survival_index != 2,]
HOLDOUT_DATA <- HOLDOUT_DATA[HOLDOUT_DATA$survival_index != 2,]
DEV_DATA$survival_index <- factor(DEV_DATA$survival_index)
HOLDOUT_DATA$survival_index <- factor(HOLDOUT_DATA$survival_index)

# Filter gene expression and mutation data based on dataset separation
gex_DEV <- gex[,colnames(gex) %in% DEV_DATA$sample_id]
gex_HOLDOUT <- gex[,colnames(gex) %in% HOLDOUT_DATA$sample_id]

# Candidate genomic features to be used as features while fitting
gex_feature_all <- TP1_gex
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
    tmp_gex_feature <- c()
    # Temporary expanded feature set
    test_feature <- c(current_feature, candidates[jj])
    # Temporaty expanded gex and mut feature set
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
      
      # Check for row rank deficiency
      
      # Check for group deficiency
      rows <- length(trainData$survival_index)
      rows_1 <- length(trainData$survival_index[trainData$survival_index==1]);
      rows_3 <- length(trainData$survival_index[trainData$survival_index==3]);
      rows_4 <- length(trainData$survival_index[trainData$survival_index==4]);
      if (rows_1<0.22*rows|rows_3<0.22*rows|rows_4<0.22*rows){
        CV_error[kk] <- NA
        next
      }
      
      # Quadratic Discriminant Analysis
      qda.trainData <- qda(survival_index ~ .-sample_id, data = trainData)
      qda.pred <- predict(qda.trainData, testData)$class
      cMat <- confusionMatrix(qda.pred,testData$survival_index)
      CV_error[kk] <- 1-as.numeric(cMat$overall[1])
    }
    if(length(CV_error[!is.na(CV_error)])==0){
      candidate_errors[jj] <- 1
    }else candidate_errors[jj] <- mean(CV_error, na.rm = T)
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
  gex_feature <- c(gex_feature, new_feature)
}
{
  par(mfrow = c(1,1))
  png(filename="./Result/QDA_4/CV_errors.png")
  plot(1:ii, CV_errors_procedure, xlab = "# of features", ylab = "5-fold CV error", xlim = c(0,ii+1))
  lines(1:ii, CV_errors_procedure)
  text(1:ii, CV_errors_procedure, labels = c(current_feature,"(Terminate)"))
  dev.off()
}

# Test error estimate using hold-out set
# Dataset condtruction
# With the feature set, make a table containing class labels and feature values of each patient.
# Full data (for fitting)
cat("Full dataset construction\n")
gex_dataset <- t(gex_DEV[gex_feature,])
gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
clin_dataset_DEV <- merge(DEV_DATA, gex_dataset, by = "sample_id", all = FALSE)
# Hold-out
cat("Hold-out dataset construction\n")
gex_dataset <- t(gex_HOLDOUT[gex_feature,])
gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
clin_dataset_HOLDOUT <- merge(HOLDOUT_DATA, gex_dataset, by = "sample_id", all = FALSE)

# Dataset ready
# Estimate test error using hold-out set
# Train by full dataset, test by hold-out subset!
# QDA model fitting
qda.Final <- qda(survival_index ~ .-sample_id, data = clin_dataset_DEV)
# Training error
qda.pred <- predict(qda.Final, clin_dataset_DEV, type = "class")
cMat <- confusionMatrix(qda.pred$class,clin_dataset_DEV$survival_index)
QDAError_T <- 1-as.numeric(cMat$overall[1])
{
  par(mfrow = c(1,1))
  # Performance
  png(filename="./Result/QDA_4/QDA_training.png")
  plot(clin_dataset_DEV$survival_index,qda.pred$class,xlab = "Training data", ylab = "Prediction by QDA", col= c(1,3,4))
  legend("topleft", legend = c("Class 1","Class 3","Class 4"), fill = c(1,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}
# Test
qda.pred <- predict(qda.Final, clin_dataset_HOLDOUT, type = "class")
cMat <- confusionMatrix(qda.pred$class,clin_dataset_HOLDOUT$survival_index)
QDAError_F <- 1-as.numeric(cMat$overall[1])
{
  cat(sprintf("\nQDA performance on test set: %2.3g\n",cMat$overall[1]))
  cat(sprintf("Feature set: "))
  cat(toString(current_feature))
  cat("\n")
}
{
  par(mfrow = c(1,1))
  # Performance
  png(filename="./Result/QDA_4/QDA_performance.png")
  plot(clin_dataset_HOLDOUT$survival_index,qda.pred$class,xlab = "Test data", ylab = "Prediction by QDA", col= c(1,3,4))
  legend("topleft", legend = c("Class 1","Class 3","Class 4"), fill = c(1,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}

# Visualization
#install.packages("klaR")
library(klaR)
{
  # Fitted decision boundary in terms of original feature space
  # https://stats.stackexchange.com/questions/143692/plotting-qda-projections-in-r
  # https://stackoverflow.com/questions/23420094/cant-plot-the-result-of-a-quadratic-discriminant-analysis-using-partimat-in-the
  X11(width=20, height=20)
  png(filename="./Result/QDA_4/Feature_3_boundary_.png")
  partimat(survival_index ~ .-sample_id, data = clin_dataset_DEV[,1:5], method = "qda",
           plot.matrix = TRUE, col.correct='blue', col.wrong='red')
  dev.off()
}
