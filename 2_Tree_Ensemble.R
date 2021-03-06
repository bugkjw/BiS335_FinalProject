#install.packages("tree")
library(tree)
library(MASS)
#install.packages("randomForest")
library(randomForest)
#install.packages("caret")
library(caret)
library(e1071)
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

# Candidate genomic features to be used as features while fitting
gex_feature_all <- as.character(gex_res$Hugo_Symbol);
gex_feature_all <- gex_feature_all[gex_feature_all %in% rownames(gex)]; factor(gex_feature_all)
mut_feature_all <- unique(as.character(mut_res$Hugo_Symbol)); factor(mut_feature_all)
feature_all <- c(gex_feature_all, mut_feature_all)

# Whole dataset condtruction
# With the current feature set (to be tested),
# make a table containing class labels and feature values of each patient.
# Training set
gex_dataset <- t(gex[gex_feature_all,])
gex_dataset <- data.frame(sample_id = rownames(gex_dataset),gex_dataset); rownames(gex_dataset) <- NULL
clin_dataset <- merge(clin_label, gex_dataset, by = "sample_id", all = FALSE)

mut_list <- unique(mut[mut$Hugo_Symbol %in% mut_feature_all, c("sample_id","Hugo_Symbol")])
mut_ids <- as.character(clin_label$sample_id);
mut_dataset <- data.frame(matrix(0,nrow = length(mut_ids), ncol = length(mut_feature_all)));
rownames(mut_dataset) <- mut_ids; colnames(mut_dataset) <- mut_feature_all
for (ll in 1:length(mut_ids)){
  for (mm in 1:length(mut_feature_all)){
    tmp_pid <- mut_ids[ll]
    tmp_gene <- toString(mut_feature_all[mm])
    if (tmp_gene %in% mut_list[mut_list$sample_id == tmp_pid,]$Hugo_Symbol){
      mut_dataset[tmp_pid,tmp_gene] <- TRUE
    }else{
      mut_dataset[tmp_pid,tmp_gene] <- FALSE
    }
  }
}
mut_dataset <- data.frame(sample_id = rownames(mut_dataset),mut_dataset); rownames(mut_dataset) <- NULL
clin_dataset <- merge(clin_dataset, mut_dataset, by = "sample_id", all = FALSE)

# Dataset separation: 80% development set, 20% hold-out set (for testing the final model)
shuffle <- clin_dataset[sample(nrow(clin_dataset)),]
shuffle <- na.omit(shuffle); rownames(shuffle) <- NULL
SELECT <- sample(nrow(shuffle)*0.2)
DEV_DATA <- shuffle[-SELECT,]
HOLDOUT_DATA <- shuffle[SELECT,]

# Bagging
set.seed(1)
# Use all the features
bag.Final <- randomForest(survival_index ~ .-sample_id
                           , data = droplevels(DEV_DATA)
                           , mtry = length(colnames(DEV_DATA))-2
                           , importance = TRUE)
# Training error
bag.pred <- predict(bag.Final, newdata = DEV_DATA, type = "class")
cMat_T <- confusionMatrix(bag.pred, DEV_DATA$survival_index)
BagError_T <- 1-as.numeric(cMat_T$overall[1])
{
  png(filename="./Result/Forest/Bagging_training.png")
  plot(DEV_DATA$survival_index
       , bag.pred
       , xlab = "Training data", ylab = "Prediction by the random forest (no predictor reduction)"
       , col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat_T$overall[1]))
  dev.off()
}
# Bagging: Test
bag.pred <- predict(bag.Final, newdata = HOLDOUT_DATA, type = "class")
cMat <- confusionMatrix(bag.pred, HOLDOUT_DATA$survival_index)
BagError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nRandom forest (no predictor reduction) performance on test set: %2.3g\n",cMat$overall[1]))
# Bagging: Visualization
{
  par(mfrow = c(1,1))
  png(filename="./Result/Forest/Bagging_VarImpPlot.png")
  importance(bag.Final)
  varImpPlot(bag.Final)
  dev.off()
}
{
  png(filename="./Result/Forest/Bagging_performance.png")
  plot(HOLDOUT_DATA$survival_index
       , bag.pred
       , xlab = "Test data", ylab = "Prediction by the random forest (no predictor reduction)"
       , col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}

# Random forest
set.seed(2)
fit.control <- tune.control(cross = 5)
tune.out <- tune.randomForest(survival_index ~ .-sample_id
                              , data = droplevels(DEV_DATA)
                              , ntree = 500
                              , mtry = seq(round(sqrt(length(colnames(DEV_DATA))-2))
                                           , length(colnames(DEV_DATA))-2
                                           , 200)
                              , tunecontrol = fit.control)
rf.best <- tune.out$best.parameters
rf.Final <- randomForest(survival_index ~ .-sample_id
                         , data = droplevels(DEV_DATA)
                         , n.trees = 500
                         , mtry = as.numeric(rf.best[1]))
# Training error
rand.pred <- predict(rf.Final, newdata = DEV_DATA,type = "class")
cMat_T <- confusionMatrix(rand.pred, DEV_DATA$survival_index)
RandError_T <- 1-as.numeric(cMat_T$overall[1])
{
  png(filename="./Result/Forest/RandomForest_training.png")
  plot(DEV_DATA$survival_index,rand.pred
       , xlab = "Training data", ylab = "Prediction by the random forest"
       , col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat_T$overall[1]))
  dev.off()
}
{
  # Tune grid
  png(filename = "./Result/Forest/RandomForest_Tunegrid.png")
  plot(tune.out$performances$mtry,tune.out$performances$error, xlab = "mtry", ylab = "5-fold CV error")
  lines(tune.out$performances$mtry,tune.out$performances$error)
  dev.off()
}
{
  par(mfrow = c(1,1))
  png(filename="./Result/Forest/RandomForest_VarImpPlot.png")
  importance(rf.Final)
  varImpPlot(rf.Final)
  dev.off()
}
# Random forest: Test
rand.pred <- predict(rf.Final, newdata = HOLDOUT_DATA,type = "class")
cMat <- confusionMatrix(rand.pred, HOLDOUT_DATA$survival_index)
RandError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nRandom forest performance on test set: %2.3g\n",cMat$overall[1]))
# Random forest: Performance
{
  png(filename="./Result/Forest/RandomForest_performance.png")
  plot(HOLDOUT_DATA$survival_index,rand.pred
       , xlab = "Test data", ylab = "Prediction by the random forest"
       , col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}

# Boosting
set.seed(3)
fit.control <- trainControl(method = "cv", number = 5)
tune.grid <- expand.grid(interaction.depth = c(1,2,3,4)
                         , n.trees = 500
                         , shrinkage = c(0.01,0.001,0.0001)
                         , n.minobsinnode = 10)
DEV_DATA_gbm <- DEV_DATA[,2:length(colnames(DEV_DATA))]
boost.Final <- train(survival_index ~ .
                    , data = droplevels(DEV_DATA_gbm)
                    , method = "gbm"
                    , trControl = fit.control
                    , tuneGrid = tune.grid
                    , verbose = FALSE)
# Training error
boost.pred <- predict(boost.Final
                      , newdata = DEV_DATA
                      , type = "raw"
                      , n.trees = 500)
cMat_T <- confusionMatrix(boost.pred, DEV_DATA$survival_index)
BoostError_T <- 1-as.numeric(cMat_T$overall[1])
{
  png(filename="./Result/Forest/Boosting_training.png")
  plot(DEV_DATA$survival_index
       , boost.pred
       , xlab = "Training data", ylab = "Prediction by the gradient boosted model"
       , col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Training misclassification error: %2.3g",1-cMat_T$overall[1]))
  dev.off()
}
# Boosting: Visualization
png(filename="./Result/Forest/Boosting_Visualization.png")
plot(boost.Final)
dev.off()
{
  par(mfrow = c(1,1))
  png(filename="./Result/Forest/Boosting_VarImpPlot.png")
  vip <- summary(boost.Final)[1:30,]; rownames(vip) <- NULL
  par(mai=c(1,1,1,1))
  barplot(rev(vip$rel.inf), names.arg = rev(vip$var), horiz = TRUE, las = 1, xlab = "Relative influence")
  title("Variable importance plot")
  dev.off()
}
# Boosting: Test
HOLDOUT_DATA_gbm <- HOLDOUT_DATA[,2:length(colnames(HOLDOUT_DATA))]
boost.pred <- predict(boost.Final
                      , newdata = HOLDOUT_DATA_gbm
                      , type = "raw"
                      , n.trees = 500)
cMat <- confusionMatrix(boost.pred, HOLDOUT_DATA_gbm$survival_index)
BoostError_F <- 1-as.numeric(cMat$overall[1])
cat(sprintf("\nBoosting performance on test set: %2.3g\n",cMat$overall[1]))

{
  png(filename="./Result/Forest/Boosting_performance.png")
  plot(HOLDOUT_DATA$survival_index
       , boost.pred
       , xlab = "Test data", ylab = "Prediction by the gradient boosted model"
       , col = c(1,2,3,4))
  legend("topleft", legend = c("Class 1","Class 2","Class 3","Class 4"), fill = c(1,2,3,4))
  title(sprintf("Test misclassification error: %2.3g",1-cMat$overall[1]))
  dev.off()
}
{
  cat(sprintf("Bagging forest training error:         %2.3f\n",BagError_T));
  cat(sprintf("Random forest training error:         %2.3f\n",RandError_T));
  cat(sprintf("Boosted training error:         %2.3f\n",BoostError_T));
  cat(sprintf("Bagging forest test error estimation:        %2.3f\n",BagError_F));
  cat(sprintf("Random forest test error estimation:         %2.3f\n",RandError_F));
  cat(sprintf("Boosted model test error estimation:         %2.3f\n",BoostError_F));
}

save.image("./Result/Forest/Forest_data.RData")
