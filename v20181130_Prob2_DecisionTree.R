#install.packages("tree")
library(tree)
library(MASS)
#install.packages("randomForest")
library(randomForest)
#install.packages("gbm")
library(gbm)

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
gex_feature <- as.character(gex_res_sorted$Hugo_Symbol[1:10]); factor(gex_feature)

mut_res_sorted <- mut_res[with(mut_res, order(mut_res$p.value)),]; rownames(mut_res_sorted) <- NULL;
mut_feature <- as.character(mut_res_sorted$Hugo_Symbol[1:30]); factor(mut_feature)

# Entire dataset condtruction: Make a table containing class labels and predictor values of each patient.
if (length(gex_feature) == 1){
  gex_dataset <- data.frame(gex[gex_feature,])
  colnames(gex_dataset) <- c(gex_feature)
}else{
  gex_dataset <- t(gex[gex_feature,]);
}
gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL

mut_list <- unique(mut[mut$Hugo_Symbol %in% mut_feature, c("sample_id","Hugo_Symbol")])
mut_ids <- as.character(clin$sample_id);
mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mut_feature)));
rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- mut_feature

for (ii in 1:length(mut_ids)){
  for (jj in 1:length(mut_feature)){
    tmp_pid <- mut_ids[ii]
    tmp_gene <- toString(mut_feature[jj])
    if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
      mut_dataset[tmp_pid,tmp_gene] <- TRUE
    }else{
      mut_dataset[tmp_pid,tmp_gene] <- FALSE
    }
  }
}
mut_dataset <- data.frame(sample_id = rownames(mut_dataset),mut_dataset); rownames(mut_dataset) <- NULL

clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)
clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)

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
  {par(mfrow = c(1,1))
    plot(tree.trainData)
    text(tree.trainData,pretty = 0)
    tree.trainData}
  # Decision tree: Test
  print("Simple decision tree performance")
  tree.pred <- predict(tree.trainData,testData,type = "class")
  print(table(tree.pred,testData$survival_index))
  TreeError[ii] <- sum(tree.pred != testData$survival_index)/length(testData$survival_index)
  
  # Pruning
  pruning <- cv.tree(tree.trainData,FUN = prune.misclass, K = 5, eps = 1e-3)
  {par(mfrow = c(1,2))
    plot(pruning$size,pruning$dev,type = "b")
    plot(pruning$k,pruning$dev,type = "b")}
  prune.misclass(tree.trainData)
  # Pruning: Test
  print("Pruned decision tree performance")
  ptree.pred <- predict(tree.trainData,testData,type = "class")
  print(table(ptree.pred,testData$survival_index))
  PTreeError[ii] <- sum(ptree.pred != testData$survival_index)/length(testData$survival_index)
  
  # Bagging
  set.seed(1)
  train <- sample(1:nrow(trainData), nrow(trainData))
  test <- trainData[-train, "survival_index"]
  # Use all the features
  bag.forest <- randomForest(survival_index ~ .-sample_id, data = droplevels(trainData),
                             subset = train, mtry = length(colnames(shuffled_dataset))-2, nTree = 500)
  # Bagging: Test
  print("Bagging forest performance")
  bag.pred <- predict(bag.forest,testData,type = "class")
  {plot(testData$survival_index,bag.pred, xlab = "Test data", ylab = "Prediction by the bagging forest"); abline(0,1)}
  print(table(bag.pred,testData$survival_index))
  BagError[ii] <- sum(bag.pred != testData$survival_index)/length(testData$survival_index)
  # Bagging: Variable importance
  importance(bag.forest)
  varImpPlot(bag.forest)
  
  # Random forest
  set.seed(1)
  rand.forest <- randomForest(survival_index ~ .-sample_id, subset = train, data = droplevels(trainData),
                              nTree = 500)
  # Random forest: Test
  print("Random forest performance")
  rand.pred <- predict(rand.forest,testData,type = "class")
  {plot(testData$survival_index,rand.pred, xlab = "Test data", ylab = "Prediction by the random forest"); abline(0,1)}
  print(table(rand.pred,testData$survival_index))
  RandError[ii] <- sum(rand.pred != testData$survival_index)/length(testData$survival_index)
  # Random forest: Variable importance
  importance(rand.forest)
  varImpPlot(rand.forest)
  
  # Boosting
  set.seed(1)
  boost.model <- gbm(survival_index ~ .-sample_id, data = trainData[train,], distribution = "gaussian",
                    n.trees = 5000)
  summary(boost.model)
  # Boosting: Test
  print("Boosted model performance")
  boost.pred <- predict(boost.model, testData, n.trees = 5000)
  print(table(round(boost.pred),testData$survival_index))
  BoostError[ii] <- sum(round(boost.pred) != droplevels(testData)$survival_index)/length(testData$survival_index)
}
cat(sprintf("Simple tree CV error estimation: %2.3f\n",mean(TreeError)));
cat(sprintf("Pruned simple tree CV error estimation: %2.3f\n",mean(PTreeError)));
cat(sprintf("Bagging forest CV error estimation: %2.3f\n",mean(BagError)));
cat(sprintf("Random forest CV error estimation: %2.3f\n",mean(RandError)));
cat(sprintf("Boosting CV error estimation: %2.3f\n",mean(BoostError)));

cat(sprintf("Simple tree CV error variability: %2.3f\n",sqrt(var(TreeError))));
cat(sprintf("Pruned simple tree CV error variability: %2.3f\n",sqrt(var(PTreeError))));
cat(sprintf("Bagging forest CV error variability: %2.3f\n",sqrt(var(BagError))));
cat(sprintf("Random forest CV error variabiity: %2.3f\n",sqrt(var(RandError))));
cat(sprintf("Boosting CV error variabiity: %2.3f\n",sqrt(var(BoostError))));
