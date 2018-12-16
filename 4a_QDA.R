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

# Use the class label found in #1, generate entire labeled dataset
labels <- unique(clin$survival_index);
clin_label <- data.frame(clin$sample_id,clin$survival_index); colnames(clin_label) <- c("sample_id","survival_index")

# QDA: deficient data in class 2 yields errors. Remove
clin_label <- clin_label[clin_label$survival_index != 2,]
clin_label$survival_index <- factor(clin_label$survival_index)

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

# QDA: deficient data in class 2 yields errors. Remove
clin_label <- clin_label[clin_label$survival_index != 2,]
clin_label$survival_index <- factor(clin_label$survival_index)

# 1. Adding influential gene features in GBM
add_feature <- Boost.varimp[TP1_gex[TP1_gex %in% Boost.varimp$var],]
add_feature <- add_feature[add_feature$rel.inf > 0,]
add_feature <- factor(as.character(add_feature$var))

# Add gene features one by one and assess 5-fold CV error
current_feature <- colnames(QDA.model$means)
added_feature <- c()
unselected_feature <- add_feature
CV_errors_procedure <- c()
CV_train_procedure <- c()

iterlim <- length(add_feature)
for(ii in 0:iterlim){
  if (ii >= 1){
    select_feature <- as.character(sample(unselected_feature,1))
    unselected_feature <- as.character(unselected_feature[unselected_feature != select_feature])
    current_feature <- c(current_feature, select_feature)
    added_feature <- c(added_feature, select_feature)
  }
  # Dataset condtruction
  # With the current feature set (to be tested),
  # make a table containing class labels and feature values of each patient.
  cat(sprintf("Adding feature %d\n", ii))
  gex_dataset <- t(gex[current_feature,])
  gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
  clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
  # Now we have the feature set and whole dataset.
  # Split the dataset for cross validation
  tmp_folds <- cut(seq(1,nrow(clin_dataset)),breaks=5,labels=FALSE)
  # 5-fold CV
  CV_error <- c(); CV_train <- c();
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
    if (rows_1<0.25*rows|rows_3<0.25*rows|rows_4<0.25*rows){
      CV_error[kk] <- NA
      next
    }
    
    # Quadratic Discriminant Analysis
    qda.trainData <- qda(survival_index ~ .-sample_id, data = trainData)
    qda.pred <- predict(qda.trainData, trainData)$class
    cMat <- confusionMatrix(qda.pred,trainData$survival_index)
    CV_train[kk] <- 1-as.numeric(cMat$overall[1])
    # Test
    qda.pred <- predict(qda.trainData, testData, type = "class")$class
    cMat <- confusionMatrix(qda.pred,testData$survival_index)
    CV_error[kk] <- 1-as.numeric(cMat$overall[1])
  }
  CV_errors_procedure[ii+1] <- mean(CV_error, na.rm = T)
  CV_train_procedure[ii+1] <- mean(CV_train, na.rm = T)
}

{
  png(filename="./Result/4/QDA_midterm_addition.png")
  par(mfrow = c(1,1))
  dat <- data.frame(train = CV_train_procedure, test = CV_errors_procedure)
  matplot(dat, type = c("b"),pch = 2,col = c(3,4)
          ,xlab = "# of features added"
          ,ylab = "Misclassification error")
  legend("center",c("Training error","5-fold CV error"),pch = 2,col = c(3,4))
  dev.off()
}

# 2. Adding random gene features
add_feature <- factor(as.character(rownames(gex)))

# Add gene features one by one and assess 5-fold CV error
current_feature <- colnames(LDA.model$means)
added_feature <- c()
unselected_feature <- add_feature
CV_errors_procedure <- c()
CV_train_procedure <- c()

for(ii in 0:iterlim){
  if (ii >= 1){
    select_feature <- as.character(sample(unselected_feature,1))
    unselected_feature <- as.character(unselected_feature[unselected_feature != select_feature])
    current_feature <- c(current_feature, select_feature)
    added_feature <- c(added_feature, select_feature)
  }
  # Dataset condtruction
  # With the current feature set (to be tested),
  # make a table containing class labels and feature values of each patient.
  cat(sprintf("Adding feature %d\n", ii))
  gex_dataset <- t(gex[current_feature,])
  gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
  clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
  # Now we have the feature set and whole dataset.
  # Split the dataset for cross validation
  tmp_folds <- cut(seq(1,nrow(clin_dataset)),breaks=5,labels=FALSE)
  # 5-fold CV
  CV_error <- c(); CV_train <- c();
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
    qda.pred <- predict(qda.trainData, trainData)$class
    cMat <- confusionMatrix(qda.pred,trainData$survival_index)
    CV_train[kk] <- 1-as.numeric(cMat$overall[1])
    # Test
    qda.pred <- predict(qda.trainData, testData, type = "class")$class
    cMat <- confusionMatrix(qda.pred,testData$survival_index)
    CV_error[kk] <- 1-as.numeric(cMat$overall[1])
  }
  CV_errors_procedure[ii+1] <- mean(CV_error, na.rm = T)
  CV_train_procedure[ii+1] <- mean(CV_train, na.rm = T)
}

{
  png(filename="./Result/4/QDA_random_addition.png")
  par(mfrow = c(1,1))
  dat <- data.frame(train = CV_train_procedure, test = CV_errors_procedure)
  matplot(dat, type = c("b"),pch = 2,col = c(3,4)
          ,xlab = "# of features added"
          ,ylab = "Misclassification error")
  legend("center",c("Training error","5-fold CV error"),pch = 2,col = c(3,4))
  dev.off()
}
