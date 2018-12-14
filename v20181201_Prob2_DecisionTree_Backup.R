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

# Data import
clin <- readRDS("./Data/clinical.rds");
gex <- readRDS("./Data/expression.rds");
mut <- readRDS("./Data/mutation.rds");
gex_res <- read.csv("./Data/gex_anova_result.csv");
mut_res <- read.csv("./Data/mut_chisq_result.csv");

# Use the class label found in #1
labels <- unique(clin$survival_index);
clin_label <- data.frame(clin$sample_id,clin$survival_index); colnames(clin_label) <- c("sample_id","survival_index")

# Select genomic features to be used as predictors
gex_res_sorted <- gex_res[with(gex_res, order(gex_res$adjust.pval)),]; rownames(gex_res_sorted) <- NULL;
mut_res_sorted <- mut_res[with(mut_res, order(mut_res$p.value)),]; rownames(mut_res_sorted) <- NULL;

gex_feature_all <- as.character(gex_res_sorted$Hugo_Symbol);
gex_feature_all <- gex_feature_all[gex_feature_all %in% rownames(gex)]; factor(gex_feature_all)
mut_feature_all <- as.character(mut_res_sorted$Hugo_Symbol);
mut_feature_all <- mut_feature_all[mut_feature_all %in% mut$Hugo_Symbol]; factor(mut_feature_all)

# Stepwise feature subset selection
# Forward selection
# 1. Try each possible features that are not in the current set
# 2. Get the feature with maximum error reduction
# 3. Add to the current feature set
# 4. Repeat until possible error difference is all >0
feature_all <- c(gex_feature_all, mut_feature_all)
gex_feature <- c(); mut_feature <- c()
current_feature <- c()
last_min_error <- 1;
for(ii in 1:length(feature_all)){
  # Candidate features: Ones not in the current feature set
  candidates <- feature_all[!(feature_all %in% current_feature)]
  CV_TreeErrors <- c()
  
  for (jj in 1:length(candidates)){
    test_gex <- c(); test_mut <- c()
    # Temporary expanded feature set
    test_feature <- c(current_feature, candidates[jj])
    # Temporaty expanded gex and mut feature set
    if (candidates[jj] %in% gex_feature_all){
      test_gex <- c(gex_feature, candidates[jj])
      test_mut <- mut_feature
    }else if (candidates[jj] %in% mut_feature_all){
      test_gex <- gex_feature
      test_mut <- c(mut_feature, candidates[jj])
    }
    
    # Test dataset condtruction: Make a table containing class labels and feature values of each patient.
    cat(sprintf("Feature %d selection, looking at candidate %d\n", ii, jj))
    if (length(test_gex) == 0){cat("No gex feature tested. Continue...\n")
      flag <- 1
    }else{
      flag <- 0
      if (length(test_gex) == 1){
        gex_dataset <- data.frame(gex[test_gex,])
        colnames(gex_dataset) <- test_gex
      }else{
        gex_dataset <- t(gex[test_gex,])
      }
      gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
      clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
    }
    
    if (length(test_mut) == 0){cat("No mut feature tested. Continue...\n")
    }else{
      mut_list <- unique(mut[mut$Hugo_Symbol %in% test_mut, c("sample_id","Hugo_Symbol")])
      mut_ids <- as.character(clin$sample_id);
      mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(test_mut)));
      rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- test_mut
      for (ll in 1:length(mut_ids)){
        for (mm in 1:length(test_mut)){
          tmp_pid <- mut_ids[ll]
          tmp_gene <- toString(test_mut[mm])
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
    shuffled_dataset <- clin_dataset[sample(nrow(clin_dataset)),]
    rownames(shuffled_dataset) <- NULL
    sample_num <- length(rownames(shuffled_dataset))
    folds <- cut(seq(1,nrow(shuffled_dataset)),breaks=5,labels=FALSE)
    
    # 5-fold CV
    test_TreeError <- c()
    for (kk in 1:5){
      testIndexes <- which(folds==kk,arr.ind=TRUE)
      testData <- shuffled_dataset[testIndexes, ]
      trainData <- shuffled_dataset[-testIndexes, ]
      
      # Decision tree
      tree.trainData <- tree(survival_index ~ .-sample_id, trainData)
      # Decision tree: Test
      #print("Simple decision tree performance")
      tree.pred <- predict(tree.trainData,testData,type = "class")
      #print(table(tree.pred,testData$survival_index))
      test_TreeError[kk] <- sum(tree.pred != testData$survival_index)/length(testData$survival_index)
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

# Some forward-stepwise-selected feature sets
A <- c("NOX4","IGSF10","ARID5A","CYYR1","PAK3")
B <- c("ARF1","KCNMB3","CD1E","HIST1H3A")
C <- c("NOX4","ARID5A","NTF5","BCDIN3D","LPHN2")
D <- c("RCL1","LILRA1","LPHN2","KIAA0152","FOS","CFP","PCLO","FIGF")
D_gex <- c("RCL1","LILRA1","LPHN2","KIAA0152","FOS","CFP","FIGF"); D_mut <- c("PCLO")
E <- c("ARF1","KCNMB3","IL27RA","SPRED1","DEFB1","PCLO")
E_gex <-c("ARF1","KCNMB3","IL27RA","SPRED1","DEFB1"); E_mut <- c("PCLO")
current_feature <- unique(c(A,B,C,D,E))
gex_feature <- unique(c(A,B,C,D_gex,E_gex))
mux_feature <- unique(c(D_mut,E_mut))
# Various enhancements
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

if (length(mux_feature) == 0){cat("No mut feature tested. Continue...\n")
}else{
  mut_list <- unique(mut[mut$Hugo_Symbol %in% mux_feature, c("sample_id","Hugo_Symbol")])
  mut_ids <- as.character(clin$sample_id);
  mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mux_feature)));
  rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- mux_feature
  for (ll in 1:length(mut_ids)){
    for (mm in 1:length(mux_feature)){
      tmp_pid <- mut_ids[ll]
      tmp_gene <- toString(mux_feature[mm])
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

# Split the dataset for cross validation
shuffled_dataset <- clin_dataset[sample(nrow(clin_dataset)),]
shuffled_dataset = shuffled_dataset[!is.na(shuffled_dataset$survival_index),]
rownames(shuffled_dataset) <- NULL
sample_num <- length(rownames(shuffled_dataset))
folds <- cut(seq(1,nrow(shuffled_dataset)),breaks=5,labels=FALSE)

# 5-fold CV
TreeError <- c()
PTreeError <- c()
BagError <- c()
RandError <- c()
for (ii in 1:5){
  testIndexes <- which(folds==ii,arr.ind=TRUE)
  testData <- shuffled_dataset[testIndexes, ]
  trainData <- shuffled_dataset[-testIndexes, ]
  
  # Decision tree
  tree.trainData <- tree(survival_index ~ .-sample_id, trainData)
  # Decision tree: Test
  tree.pred <- predict(tree.trainData,testData,type = "class")
  cMat <- confusionMatrix(tree.pred,testData$survival_index)
  TreeError <- c(TreeError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nSimple decision tree performance: %2.3g\n",cMat$overall[1]))
  # Plot
  {par(mfrow = c(1,2))
    plot(tree.trainData)
    text(tree.trainData,pretty = 0)
    plot(testData$survival_index,tree.pred, xlab = "Test data", ylab = "Prediction by the simple decision tree"); abline(0,1)}
  
  # Pruning
  pruning <- cv.tree(tree.trainData, FUN = prune.misclass, eps = 1e-3)
  Best <- pruning$size[pruning$dev == min(pruning$dev)]
  ptree.trainData <- prune.misclass(tree.trainData, best = Best[1])
  # Pruning: Test
  ptree.pred <- predict(ptree.trainData,testData,type = "class")
  cMat <- confusionMatrix(ptree.pred,testData$survival_index)
  PTreeError <- c(PTreeError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nPruned decision tree performance: %2.3g\n",cMat$overall[1]))
  # Plot
  {par(mfrow = c(1,2))
    plot(pruning$size,pruning$dev,type = "b")
    plot(pruning$k,pruning$dev,type = "b")}
  {par(mfrow = c(1,2))
    plot(ptree.trainData)
    text(ptree.trainData,pretty = 0)
    plot(testData$survival_index,ptree.pred, xlab = "Test data", ylab = "Prediction by the pruned decision tree"); abline(0,1)}
  
  # Bagging
  set.seed(1)
  train <- sample(1:nrow(trainData), nrow(trainData))
  test <- trainData[-train, "survival_index"]
  # Use all the features
  bag.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(trainData),
                             mtry = length(colnames(shuffled_dataset))-2, importance = TRUE)
  # Bagging: Test
  bag.pred <- predict(bag.forest, newdata = testData, type = "class")
  cMat <- confusionMatrix(bag.pred,testData$survival_index)
  BagError <- c(BagError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nRandom forest (no predictor reduction) performance: %2.3g\n",cMat$overall[1]))
  # Variable importance
  {par(mfrow = c(1,1))
    importance(bag.forest)
    varImpPlot(bag.forest)
    plot(testData$survival_index,bag.pred, xlab = "Test data", ylab = "Prediction by the random forest (no predictor reduction)"); abline(0,1)}
  
  # Random forest
  set.seed(1)
  rand.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(trainData),
                              importance = TRUE)
  # Random forest: Test
  rand.pred <- predict(rand.forest,newdata = testData,type = "class")
  cMat <- confusionMatrix(rand.pred,testData$survival_index)
  RandError <- c(RandError,1-as.numeric(cMat$overall[1]))
  cat(sprintf("\nRandom forest performance: %2.3g\n",cMat$overall[1]))
  # Random forest: Variable importance
  {par(mfrow = c(1,1))
    importance(rand.forest)
    varImpPlot(rand.forest)
    plot(testData$survival_index,rand.pred, xlab = "Test data", ylab = "Prediction by the random forest"); abline(0,1)}
}
{
  cat(sprintf("Simple tree CV error estimation:         %2.3f\n",mean(TreeError)));
  cat(sprintf("Pruned simple tree CV error estimation:  %2.3f\n",mean(PTreeError)));
  cat(sprintf("Bagging forest CV error estimation:      %2.3f\n",mean(BagError)));
  cat(sprintf("Random forest CV error estimation:       %2.3f\n",mean(RandError)));
}

