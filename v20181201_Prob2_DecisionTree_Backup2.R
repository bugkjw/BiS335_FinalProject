#install.packages("tree")
library(tree)
library(MASS)
#install.packages("randomForest")
library(randomForest)
#install.packages("gbm")
library(gbm)
#install.packages("caret")
library(caret)

setwd("D:/윈도우계정/Desktop/!/3학년3가을학기/BiS335 Biomedical Statistics & Statistical Learning/Final Project/Finalterm-Project")
load("./v20181202_Tree_Final_Model.RData")

# Data import
clin <- readRDS("./Data/clinical.rds");
gex <- readRDS("./Data/expression.rds");
mut <- readRDS("./Data/mutation.rds");
gex_res <- read.csv("./Data/gex_anova_result.csv");
mut_res <- read.csv("./Data/mut_chisq_result.csv");

# Use the class label found in #1
labels <- unique(clin$survival_index);
clin_label <- data.frame(clin$sample_id,clin$survival_index); colnames(clin_label) <- c("sample_id","survival_index")

# Set candidate genomic features to be used as predictors
gex_res_sorted <- gex_res[with(gex_res, order(gex_res$adjust.pval)),]; rownames(gex_res_sorted) <- NULL;
mut_res_sorted <- mut_res[with(mut_res, order(mut_res$p.value)),]; rownames(mut_res_sorted) <- NULL;

gex_feature_all <- as.character(gex_res_sorted$Hugo_Symbol);
gex_feature_all <- gex_feature_all[gex_feature_all %in% rownames(gex)]; factor(gex_feature_all)
mut_feature_all <- as.character(mut_res_sorted$Hugo_Symbol);
mut_feature_all <- mut_feature_all[mut_feature_all %in% mut$Hugo_Symbol]; factor(mut_feature_all)
feature_all <- c(gex_feature_all, mut_feature_all)

# Overall scheme
# 1. Split the data into 5 sets
# 2. For each set (fold), define test (1/5) and training (4/5) data set
# 3. For each fold, perform feature selection and model fitting on the training set
# 4. Estimate the model missclassification error and repeat for all 5 folds
#    At each fold, vote for the selected features
# 5. Estimate the overall accuracy
# 6. Final procedure: Select a particular number of features from the highest votes
#    And train the model on the complete data. Finalize the model

# Split the dataset into 5 sets for cross validation
shuffled_clin <- clin[sample(nrow(clin)),]
shuffled_clin = shuffled_clin[!is.na(shuffled_clin$survival_index),]
rownames(shuffled_clin) <- NULL
sample_num <- length(rownames(shuffled_clin))
folds <- cut(seq(1,nrow(shuffled_clin)),breaks=5,labels=FALSE)

# Iterate over the folds
vote <- data.frame(feature_all, n = 0)
TreeError <- c()
PTreeError <- c()
BagError <- c()
RandError <- c()
for (ff in 1:5){
  testIndexes <- which(folds==ff,arr.ind=TRUE)
  test_clin <- shuffled_clin[testIndexes,]
  train_clin <- shuffled_clin[-testIndexes,]
  
  gex_test <- gex[,colnames(gex) %in% test_clin$sample_id]
  gex_train <- gex[,colnames(gex) %in% train_clin$sample_id]
  mut_test <-mut[mut$sample_id %in% test_clin$sample_id,]
  mut_train <-mut[mut$sample_id %in% train_clin$sample_id,]
  
  # Stepwise feature subset selection: Forward selection
  # 1. Try each possible features that are not in the current set
  # 2. Get the feature with maximum error reduction
  # 3. Add to the current feature set
  # 4. Repeat until possible error deviation is all >0
  gex_feature <- c(); mut_feature <- c()
  current_feature <- c()
  last_min_error <- 1;
  for(ii in 1:length(feature_all)){
    # Candidate features: Ones not in the current feature set
    candidates <- feature_all[!(feature_all %in% current_feature)]
    CV_TreeErrors <- c()
    
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
      # Dataset condtruction: Make a table containing class labels and feature values of each patient.
      cat(sprintf("Feature %d selection, looking at candidate %d\n", ii, jj))
      if (length(tmp_gex_feature) == 0){cat("No gex feature tested. Continue...\n")
        flag <- 1
      }else{
        flag <- 0
        if (length(tmp_gex_feature) == 1){
          gex_dataset <- data.frame(gex_train[tmp_gex_feature,])
          colnames(gex_dataset) <- tmp_gex_feature
        }else{
          gex_dataset <- t(gex_train[tmp_gex_feature,])
        }
        gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
        clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
      }
      
      if (length(tmp_mut_feature) == 0){cat("No mut feature tested. Continue...\n")
      }else{
        mut_list <- unique(mut_train[mut_train$Hugo_Symbol %in% tmp_mut_feature, c("sample_id","Hugo_Symbol")])
        mut_ids <- as.character(train_clin$sample_id);
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
        if (flag == 1){clin_dataset <- merge(clin_label, mut_dataset, by = "sample_id", all = FALSE)
        }else if (flag == 0){clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)}
      }
      # Now we have the feature set and whole dataset.
      # Split the dataset for cross validation
      tmp_shuffled <- clin_dataset[sample(nrow(clin_dataset)),]
      rownames(tmp_shuffled) <- NULL
      tmp_folds <- cut(seq(1,nrow(tmp_shuffled)),breaks=5,labels=FALSE)
      
      # 5-fold CV
      test_TreeError <- c()
      for (kk in 1:5){
        tmp_testIndexes <- which(tmp_folds==kk,arr.ind=TRUE)
        testData <- tmp_shuffled[tmp_testIndexes, ]
        trainData <- tmp_shuffled[-tmp_testIndexes, ]
        
        # Decision tree
        tree.trainData <- tree(survival_index ~ .-sample_id, trainData)
        # Bagging
        #set.seed(1)
        #tree.trainData <- randomForest(survival_index ~ .-sample_id, data = droplevels(trainData),
        #                           mtry = length(colnames(trainData))-2, importance = TRUE, ntree = 200)
        # Test
        tree.pred <- predict(tree.trainData,testData,type = "class")
        cMat <- confusionMatrix(tree.pred,testData$survival_index)
        test_TreeError[kk] <- 1-as.numeric(cMat$overall[1])
      }
      CV_TreeErrors[jj] <- mean(test_TreeError)
    }
    
    # Check for improvement
    if (min(CV_TreeErrors) > last_min_error){break}
    # Update threshold
    last_min_error <- min(CV_TreeErrors)
    # Greedy feature selection
    new_feature <- candidates[CV_TreeErrors == min(CV_TreeErrors)]
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
  
  # Feature selection done
  # Re-cook the training data
  if (length(gex_feature) == 0){cat("No gex feature tested. Continue...\n")
    flag <- 1
  }else{
    flag <- 0
    if (length(gex_feature) == 1){
      gex_dataset <- data.frame(gex_train[gex_feature,])
      colnames(gex_dataset) <- gex_feature
    }else{
      gex_dataset <- t(gex_train[gex_feature,])
    }
    gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
    clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
  }
  if (length(mut_feature) == 0){cat("No mut feature tested. Continue...\n")
  }else{
    mut_list <- unique(mut_train[mut_train$Hugo_Symbol %in% mut_feature, c("sample_id","Hugo_Symbol")])
    mut_ids <- as.character(train_clin$sample_id);
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
    if (flag == 1){clin_dataset <- merge(clin_label, mut_dataset, by = "sample_id", all = FALSE)
    }else if (flag == 0){clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)}
  }
  # Re-cook the test data
  if (length(gex_feature) == 0){cat("No gex feature tested. Continue...\n")
    flag <- 1
  }else{
    flag <- 0
    if (length(gex_feature) == 1){
      gex_dataset_test <- data.frame(gex_test[gex_feature,])
      colnames(gex_dataset_test) <- gex_feature
    }else{
      gex_dataset_test <- t(gex_test[gex_feature,])
    }
    gex_dataset_test <- data.frame(sample_id = rownames(gex_dataset_test),gex_dataset_test); rownames(gex_dataset_test) <- NULL
    clin_dataset_test <- merge(clin_label, gex_dataset_test, by = "sample_id", all = FALSE)
  }
  if (length(mut_feature) == 0){cat("No mut feature tested. Continue...\n")
  }else{
    mut_list <- unique(mut_test[mut_test$Hugo_Symbol %in% mut_feature, c("sample_id","Hugo_Symbol")])
    mut_ids <- as.character(test_clin$sample_id);
    mut_dataset_test <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mut_feature)));
    rownames(mut_dataset_test) <- mut_ids; colnames(mut_dataset_test) <- mut_feature
    for (ll in 1:length(mut_ids)){
      for (mm in 1:length(mut_feature)){
        tmp_pid <- mut_ids[ll]
        tmp_gene <- toString(mut_feature[mm])
        if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
          mut_dataset_test[tmp_pid,tmp_gene] <- TRUE
        }else{
          mut_dataset_test[tmp_pid,tmp_gene] <- FALSE
        }
      }
    }
    mut_dataset_test <- data.frame(sample_id = rownames(mut_dataset_test),mut_dataset_test); rownames(mut_dataset_test) <- NULL
    if (flag == 1){clin_dataset_test <- merge(clin_label, mut_dataset_test, by = "sample_id", all = FALSE)
    }else if (flag == 0){clin_dataset_test <- merge(clin_dataset_test, mut_dataset_test, by = "sample_id", all = FALSE)}
  }
  
  # Decision tree model fitting
  tree.trainData <- tree(survival_index ~ .-sample_id, clin_dataset)
  # Decision tree: Test
  tree.pred <- predict(tree.trainData,clin_dataset_test,type = "class")
  cMat <- confusionMatrix(tree.pred,clin_dataset_test$survival_index)
  TreeError <- c(TreeError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nSimple decision tree performance: %2.3g\n",cMat$overall[1]))
  # Plot
  {par(mfrow = c(1,2))
    plot(tree.trainData)
    text(tree.trainData,pretty = 0)
    plot(clin_dataset_test$survival_index,tree.pred, xlab = "Test data", ylab = "Prediction by the simple decision tree"); abline(0,1)}
  
  # Pruning
  pruning <- cv.tree(tree.trainData, FUN = prune.misclass, eps = 1e-3)
  Best <- pruning$size[pruning$dev == min(pruning$dev)]
  ptree.trainData <- prune.misclass(tree.trainData, best = Best[1])
  # Pruning: Test
  ptree.pred <- predict(ptree.trainData,clin_dataset_test,type = "class")
  cMat <- confusionMatrix(ptree.pred,clin_dataset_test$survival_index)
  PTreeError <- c(PTreeError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nPruned decision tree performance: %2.3g\n",cMat$overall[1]))
  # Plot
  {par(mfrow = c(1,2))
    plot(pruning$size,pruning$dev,type = "b")
    plot(pruning$k,pruning$dev,type = "b")}
  {par(mfrow = c(1,2))
    plot(ptree.trainData)
    text(ptree.trainData,pretty = 0)
    plot(clin_dataset_test$survival_index,ptree.pred, xlab = "Test data", ylab = "Prediction by the pruned decision tree"); abline(0,1)}
  
  # Bagging
  set.seed(1)
  # Use all the features
  bag.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(clin_dataset),
                             mtry = length(colnames(clin_dataset))-2, importance = TRUE)
  # Bagging: Test
  bag.pred <- predict(bag.forest, newdata = clin_dataset_test, type = "class")
  cMat <- confusionMatrix(bag.pred,clin_dataset_test$survival_index)
  BagError <- c(BagError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nRandom forest (no predictor reduction) performance: %2.3g\n",cMat$overall[1]))
  # Variable importance
  {par(mfrow = c(1,1))
    importance(bag.forest)
    varImpPlot(bag.forest)
    plot(clin_dataset_test$survival_index,bag.pred, xlab = "Test data", ylab = "Prediction by the random forest (no predictor reduction)"); abline(0,1)}
  
  # Random forest
  set.seed(1)
  rand.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(clin_dataset),
                              importance = TRUE)
  # Random forest: Test
  rand.pred <- predict(rand.forest,newdata = clin_dataset_test,type = "class")
  cMat <- confusionMatrix(rand.pred,clin_dataset_test$survival_index)
  RandError <- c(RandError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nRandom forest performance: %2.3g\n",cMat$overall[1]))
  # Random forest: Variable importance
  {par(mfrow = c(1,1))
    importance(rand.forest)
    varImpPlot(rand.forest)
    plot(clin_dataset_test$survival_index,rand.pred, xlab = "Test data", ylab = "Prediction by the random forest"); abline(0,1)}

  vote[vote$feature_all %in% current_feature,]$n = vote[vote$feature_all %in% current_feature,]$n + 1
}

# Finalization
# Feature selection
FeatureNum <- length(vote[vote$n!=0,]$n)
Ordered_feature <- vote[order(-vote$n),]
rownames(Ordered_feature) <- NULL
Ordered_feature <- Ordered_feature[1:FeatureNum,]
{
  cat(sprintf("Feature set: "))
  cat(toString(Ordered_feature[1:FeatureNum,]$feature_all))
  cat(sprintf("\nSimple tree CV error estimation:         %2.3f\n",mean(TreeError)));
  cat(sprintf("Pruned simple tree CV error estimation:  %2.3f\n",mean(PTreeError)));
  cat(sprintf("Bagging forest CV error estimation:      %2.3f\n",mean(BagError)));
  cat(sprintf("Random forest CV error estimation:       %2.3f\n",mean(RandError)));
}
# Model fitting
# Final data cooking...
gex_feature <- as.character(Ordered_feature[Ordered_feature$feature_all %in% rownames(gex),]$feature_all); factor(gex_feature)
mut_feature <- as.character(Ordered_feature[!(Ordered_feature$feature_all %in% rownames(gex)),]$feature_all); factor(mut_feature)
if (length(gex_feature) == 0){cat("No gex feature tested. Continue...\n")
  flag <- 1
}else{
  flag <- 0
  if (length(gex_feature) == 1){
    gex_dataset <- data.frame(gex[gex_feature,])
    colnames(gex_dataset) <- gex_feature
  }else{
    gex_dataset <- t(gex[gex_feature,])
  }
  gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
  clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
}
if (length(mut_feature) == 0){cat("No mut feature tested. Continue...\n")
}else{
  mut_list <- unique(mut_test[mut_test$Hugo_Symbol %in% mut_feature, c("sample_id","Hugo_Symbol")])
  mut_ids <- as.character(clin$sample_id);
  mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mut_feature)));
  rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- mut_feature
  for (ll in 1:length(mut_ids)){
    for (mm in 1:length(mut_feature)){
      tmp_pid <- mut_ids[ll]
      tmp_gene <- toString(mut_feature[mm])
      if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
        mut_dataset_test[tmp_pid,tmp_gene] <- TRUE
      }else{
        mut_dataset_test[tmp_pid,tmp_gene] <- FALSE
      }
    }
  }
  mut_dataset <- data.frame(sample_id = rownames(mut_dataset),mut_dataset); rownames(mut_dataset) <- NULL
  if (flag == 1){clin_dataset <- merge(clin_label, mut_dataset, by = "sample_id", all = FALSE)
  }else if (flag == 0){clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)}
}

# Decision tree model fitting
tree.Final <- tree(survival_index ~ .-sample_id, clin_dataset)
# Decision tree: Test
tree.pred <- predict(tree.Final,clin_dataset,type = "class")
cMat <- confusionMatrix(tree.pred,clin_dataset$survival_index)
TreeError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nSimple decision tree performance on training set: %2.3g\n",cMat$overall[1]))
# Plot
{par(mfrow = c(1,2))
  plot(tree.Final)
  text(tree.Final,pretty = 0)
  plot(clin_dataset$survival_index,tree.pred, xlab = "Training data", ylab = "Prediction by the simple decision tree"); abline(0,1)}

# Pruning
pruning <- cv.tree(tree.Final, FUN = prune.misclass, eps = 1e-3)
Best <- pruning$size[pruning$dev == min(pruning$dev)]
tree.Final <- prune.misclass(tree.Final, best = Best[1])
# Pruning: Test
ptree.pred <- predict(tree.Final,clin_dataset,type = "class")
cMat <- confusionMatrix(ptree.pred,clin_dataset$survival_index)
PTreeError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nPruned decision tree performance on training set: %2.3g\n",cMat$overall[1]))
# Plot
{par(mfrow = c(1,2))
  plot(pruning$size,pruning$dev,type = "b")
  plot(pruning$k,pruning$dev,type = "b")}
{par(mfrow = c(1,2))
  plot(tree.Final)
  text(tree.Final,pretty = 0)
  plot(clin_dataset$survival_index,ptree.pred, xlab = "Test data", ylab = "Prediction by the pruned decision tree"); abline(0,1)}

# Bagging
set.seed(1)
# Use all the features
bag.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(clin_dataset),
                           mtry = length(colnames(clin_dataset))-2, importance = TRUE)
# Bagging: Test
bag.pred <- predict(bag.forest, newdata = clin_dataset, type = "class")
cMat <- confusionMatrix(bag.pred,clin_dataset$survival_index)
BagError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nRandom forest (no predictor reduction) performance on training set: %2.3g\n",cMat$overall[1]))
# Variable importance
{par(mfrow = c(1,1))
  importance(bag.forest)
  varImpPlot(bag.forest)
  plot(clin_dataset$survival_index,bag.pred, xlab = "Test data", ylab = "Prediction by the random forest (no predictor reduction)"); abline(0,1)}

# Random forest
set.seed(1)
rand.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(clin_dataset),
                            importance = TRUE)
# Random forest: Test
rand.pred <- predict(rand.forest,newdata = clin_dataset,type = "class")
cMat <- confusionMatrix(rand.pred,clin_dataset$survival_index)
RandError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nRandom forest performance on training set: %2.3g\n",cMat$overall[1]))
# Random forest: Variable importance
{par(mfrow = c(1,1))
  importance(rand.forest)
  varImpPlot(rand.forest)
  plot(clin_dataset$survival_index,rand.pred, xlab = "Test data", ylab = "Prediction by the random forest"); abline(0,1)}

{
  cat(sprintf("Feature set: "))
  cat(toString(Ordered_feature[1:FeatureNum,]$feature_all))
  cat(sprintf("\nSimple tree training set error:           %2.3f\n",TreeError_F));
  cat(sprintf("Pruned simple tree training set error:    %2.3f\n",PTreeError_F));
  cat(sprintf("Bagging forest training set error:        %2.3f\n",BagError_F));
  cat(sprintf("Random forest training set error:         %2.3f\n",RandError_F));
}
