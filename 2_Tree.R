#install.packages("tree")
library(tree)
library(MASS)
#install.packages("randomForest")
library(randomForest)
#install.packages("caret")
library(caret)
#install.packages("doParallel")
library(doParallel)
registerDoParallel(4)
getDoParWorkers()


setwd("C:/Users/VSlab#10/Desktop/JinwooKim/BiS335_FinalProject_Folder")

# Data import
clin <- readRDS("./Data/clinical.rds");
gex <- readRDS("./Data/expression.rds");
mut <- readRDS("./Data/mutation.rds");
gex_res <- read.csv("./Data/gex_anova_result.csv");
mut_res <- read.csv("./Data/mut_chisq_result.csv");

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
mut_DEV <- mut[mut$sample_id %in% DEV_DATA$sample_id,]
mut_HOLDOUT <- mut[mut$sample_id %in% HOLDOUT_DATA$sample_id,]

# Candidate genomic features to be used as features while fitting
gex_feature_all <- as.character(gex_res$Hugo_Symbol);
gex_feature_all <- gex_feature_all[gex_feature_all %in% rownames(gex_DEV)]; factor(gex_feature_all)
mut_feature_all <- as.character(mut_res$Hugo_Symbol);
mut_feature_all <- mut_feature_all[mut_feature_all %in% mut_DEV$Hugo_Symbol]; factor(mut_feature_all)
feature_all <- c(gex_feature_all, mut_feature_all)

# Feature selection by forward stepwise selection
# For each feature set expansion, assess 5-fold CV misclassification error rate
# Expand until no more error rate reduction is possible

# Specifically,
# 1. Try each possible features that are not in the current set
# 2. Get the feature with maximum error (5-fold CV) reduction
# 3. Add to the current feature set
# 4. Repeat until possible error deviation is all >0
gex_feature <- c(); mut_feature <- c()
current_feature <- c()
CV_errors_procedure <- c()
last_min_error <- 1;
for(ii in 1:length(feature_all)){
  # Candidate features: Ones not in the current feature set
  candidates <- feature_all[!(feature_all %in% current_feature)]
  candidate_errors <- c()
  
  for (jj in 1:length(candidates)){
    tmp_gex_feature <- c(); tmp_mut_feature <- c()
    # Temporary expanded feature set
    test_feature <- c(current_feature, candidates[jj])
    # Temporaty expanded gex and mut feature set
    if (candidates[jj] %in% gex_feature_all){
      tmp_gex_feature <- c(gex_feature, candidates[jj])
      tmp_mut_feature <- mut_feature
    }else if (candidates[jj] %in% mut_feature_all){
      tmp_gex_feature <- gex_feature
      tmp_mut_feature <- c(mut_feature, candidates[jj])
    }
    # Dataset construction
    # With the current feature set (to be tested),
    # make a table containing class labels and feature values of each patient.
    cat(sprintf("Feature %d selection, looking at candidate %d\n", ii, jj))
    if (length(tmp_gex_feature) == 0){cat("No gex feature tested. Continue...\n")
      flag <- 1
    }else{
      flag <- 0
      if (length(tmp_gex_feature) == 1){
        gex_dataset <- data.frame(gex_DEV[tmp_gex_feature,])
        colnames(gex_dataset) <- tmp_gex_feature
      }else{
        gex_dataset <- t(gex_DEV[tmp_gex_feature,])
      }
      gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
      clin_dataset <- merge(DEV_DATA, gex_dataset, by = "sample_id", all = FALSE)
    }
    if (length(tmp_mut_feature) == 0){cat("No mut feature tested. Continue...\n")
    }else{
      mut_list <- unique(mut_DEV[mut_DEV$Hugo_Symbol %in% tmp_mut_feature, c("sample_id","Hugo_Symbol")])
      mut_ids <- as.character(DEV_DATA$sample_id);
      mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(tmp_mut_feature)));
      rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- tmp_mut_feature
      for (ll in 1:length(mut_ids)){
        for (mm in 1:length(tmp_mut_feature)){
          tmp_pid <- mut_ids[ll]
          tmp_gene <- toString(tmp_mut_feature[mm])
          if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
            mut_dataset[tmp_pid,tmp_gene] <- TRUE
          }else{
            mut_dataset[tmp_pid,tmp_gene] <- FALSE
          }
        }
      }
      mut_dataset <- data.frame(sample_id = rownames(mut_dataset),mut_dataset); rownames(mut_dataset) <- NULL
      if (flag == 1){clin_dataset <- merge(DEV_DATA, mut_dataset, by = "sample_id", all = FALSE)
      }else if (flag == 0){clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)}
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
      
      # Decision tree
      tree.trainData <- tree(survival_index ~ .-sample_id, trainData)
      # Test
      tree.pred <- predict(tree.trainData,testData,type = "class")
      cMat <- confusionMatrix(tree.pred,testData$survival_index)
      CV_error[kk] <- 1-as.numeric(cMat$overall[1])
    }
    candidate_errors[jj] <- mean(CV_error)
  }
  CV_errors_procedure[ii] <- min(candidate_errors)
  # Check for improvement
  if (min(candidate_errors) >= last_min_error){break}
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
  }else if (new_feature %in% mut_feature_all){
    mut_feature <- c(mut_feature, new_feature)
  }
}
{
  png(filename="./Result/Tree/CV_errors.png")
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
if (length(gex_feature) == 0){cat("No gex feature tested. Continue...\n")
  flag <- 1
}else{
  flag <- 0
  if (length(gex_feature) == 1){
    gex_dataset <- data.frame(gex_DEV[gex_feature,])
    colnames(gex_dataset) <- gex_feature
  }else{
    gex_dataset <- t(gex_DEV[gex_feature,])
  }
  gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
  clin_dataset_DEV <- merge(DEV_DATA, gex_dataset, by = "sample_id", all = FALSE)
}
if (length(mut_feature) == 0){cat("No mut feature tested. Continue...\n")
}else{
  mut_list <- unique(mut[mut$Hugo_Symbol %in% mut_feature, c("sample_id","Hugo_Symbol")])
  mut_ids <- as.character(DEV_DATA$sample_id);
  mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mut_feature)));
  rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- mut_feature
  for (ll in 1:length(mut_ids)){
    for (mm in 1:length(mut_feature)){
      tmp_pid <- mut_ids[ll]
      tmp_gene <- toString(mut_feature[mm])
      if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
        mut_dataset[tmp_pid,tmp_gene] <- TRUE
      }else{
        mut_dataset[tmp_pid,tmp_gene] <- FALSE
      }
    }
  }
  mut_dataset <- data.frame(sample_id = rownames(mut_dataset),mut_dataset); rownames(mut_dataset) <- NULL
  if (flag == 1){clin_dataset_DEV <- merge(DEV_DATA, mut_dataset, by = "sample_id", all = FALSE)
  }else if (flag == 0){clin_dataset_DEV <- merge(clin_dataset_DEV, mut_dataset, by = "sample_id", all = FALSE)}
}
# Hold-out
cat("Hold-out dataset construction\n")
if (length(gex_feature) == 0){cat("No gex feature tested. Continue...\n")
  flag <- 1
}else{
  flag <- 0
  if (length(gex_feature) == 1){
    gex_dataset <- data.frame(gex_HOLDOUT[gex_feature,])
    colnames(gex_dataset) <- gex_feature
  }else{
    gex_dataset <- t(gex_HOLDOUT[gex_feature,])
  }
  gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
  clin_dataset_HOLDOUT <- merge(HOLDOUT_DATA, gex_dataset, by = "sample_id", all = FALSE)
}
if (length(mut_feature) == 0){cat("No mut feature tested. Continue...\n")
}else{
  mut_list <- unique(mut[mut$Hugo_Symbol %in% mut_feature, c("sample_id","Hugo_Symbol")])
  mut_ids <- as.character(HOLDOUT_DATA$sample_id);
  mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mut_feature)));
  rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- mut_feature
  for (ll in 1:length(mut_ids)){
    for (mm in 1:length(mut_feature)){
      tmp_pid <- mut_ids[ll]
      tmp_gene <- toString(mut_feature[mm])
      if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
        mut_dataset[tmp_pid,tmp_gene] <- TRUE
      }else{
        mut_dataset[tmp_pid,tmp_gene] <- FALSE
      }
    }
  }
  mut_dataset <- data.frame(sample_id = rownames(mut_dataset),mut_dataset); rownames(mut_dataset) <- NULL
  if (flag == 1){clin_dataset_HOLDOUT <- merge(HOLDOUT_DATA, mut_dataset, by = "sample_id", all = FALSE)
  }else if (flag == 0){clin_dataset_HOLDOUT <- merge(clin_dataset_HOLDOUT, mut_dataset, by = "sample_id", all = FALSE)}
}

# Dataset ready
# Estimate test error using hold-out set
# Train by full dataset, test by hold-out subset!
# Decision tree model fitting
tree.Final <- tree(survival_index ~ .-sample_id, clin_dataset_DEV)
# Training error
tree.pred <- predict(tree.Final,clin_dataset_DEV,type = "class")
cMat_T <- confusionMatrix(tree.pred,clin_dataset_DEV$survival_index)
TreeError_T <- 1-as.numeric(cMat_T$overall[1])
{
  png(filename="./Result/Tree/Simple_training.png")
  plot(clin_dataset_DEV$survival_index,tree.pred, xlab = "Training data", ylab = "Prediction by the simple decision tree", col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat_T$overall[1]))
  dev.off()
}
# Decision tree: Test
tree.pred <- predict(tree.Final,clin_dataset_HOLDOUT,type = "class")
cMat <- confusionMatrix(tree.pred,clin_dataset_HOLDOUT$survival_index)
TreeError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nSimple decision tree performance on test set: %2.3g\n",cMat$overall[1]))
# Plot
{
  par(mfrow = c(1,2))
  png(filename="./Result/Tree/Simple_trained.png")
  plot(tree.Final)
  text(tree.Final,pretty = 0)
  dev.off()
}
{
  png(filename="./Result/Tree/Simple_performance.png")
  plot(clin_dataset_HOLDOUT$survival_index,tree.pred, xlab = "Test data", ylab = "Prediction by the simple decision tree", col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}

# Pruning
pruning <- cv.tree(tree.Final, FUN = prune.misclass, eps = 1e-3, K = 5)
Best <- pruning$size[pruning$dev == min(pruning$dev)]
Ptree.Final <- prune.misclass(tree.Final, best = Best[1])
# Training error
Ptree.pred <- predict(Ptree.Final,clin_dataset_DEV,type = "class")
cMat_T <- confusionMatrix(Ptree.pred,clin_dataset_DEV$survival_index)
TreeError_T <- 1-as.numeric(cMat_T$overall[1])
{
  png(filename="./Result/Tree/Pruned_training.png")
  plot(clin_dataset_DEV$survival_index,Ptree.pred, xlab = "Training data", ylab = "Prediction by the pruned decision tree", col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat_T$overall[1]))
  dev.off()
}
# Pruning: Test
ptree.pred <- predict(tree.Final,clin_dataset_HOLDOUT,type = "class")
cMat <- confusionMatrix(ptree.pred,clin_dataset_HOLDOUT$survival_index)
PTreeError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nPruned decision tree performance on test set: %2.3g\n",cMat$overall[1]))
# Plot
{
  par(mfrow = c(1,2))
  png(filename="./Result/Tree/Pruning.png")
  plot(pruning$size,pruning$dev,type = "b")
  plot(pruning$k,pruning$dev,type = "b")
  dev.off()
}
{
  png(filename="./Result/Tree/Pruned_tree.png")
  plot(Ptree.Final)
  text(Ptree.Final,pretty = 0)
  dev.off()
}
{
  png(filename="./Result/Tree/Pruned_performance.png")
  plot(clin_dataset_HOLDOUT$survival_index,ptree.pred, xlab = "Test data", ylab = "Prediction by the pruned decision tree", col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}
{
  cat(sprintf("Feature set: "))
  cat(toString(current_feature))
  cat(sprintf("\nSimple tree test error estimation:           %2.3f\n",TreeError_F));
  cat(sprintf("Pruned simple tree test error estimation:    %2.3f\n",PTreeError_F));
}

save.image("./Result/Tree/Tree_data.RData")

