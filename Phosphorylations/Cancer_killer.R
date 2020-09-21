library(tidyverse)
library(DataExplorer)
library(caret)
library(e1071)

#### Phosphorylations

## load in data 
cancer.train <- read_csv("/Users/graceedriggs/Documents/Stat 495/Phosphorylations/phosphorylation1433/train.csv")
cancer.test <- read_csv("/Users/graceedriggs/Documents/Stat 495/Phosphorylations/phosphorylation1433/test.csv")

cancer <- bind_rows(cancer.test, cancer.train)
str(cancer)

plot_missing(cancer)

## correlation
corrgram::corrgram(cancer[!is.na(cancer$Consensus),])

##linear regression imputation for Consensus
protein.lm <- lm(Consensus ~ SVM, data = cancer)
cancer$Consensus[is.na(cancer$Consensus)] <- predict.lm(protein.lm, newdata = cancer [is.na(cancer$Consensus),])

##linear regression imputation for
protein.pssm.lm <- lm(PSSM ~ SVM + ANN + Consensus, data = cancer)
cancer$PSSM[is.na(cancer$PSSM)] <- predict.lm(protein.pssm.lm, newdata = cancer [is.na(cancer$PSSM),])

plot_missing(cancer)

## plot the density functions for each.


#split
cancer.train <- cancer %>% filter(!is.na(Response))
cancer.test <- cancer %>% filter(is.na(Response))
str(cancer.train)
cancer.train$Response <- as.factor(cancer.train$Response)

#Stochastic Gradient Boosting
#10 folds repeat 3 times
control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=5)
#Metric compare model is Accuracy
metric <- "Accuracy"
set.seed(123)
rf_default <- train(Response ~ ., 
                    data=cancer.train, 
                    method='gbm', 
                    trControl=control)

names(rf_default)
plot(rf_default)
rf_default$bestTune
predict(rf_default, newdata=cancer.test)


#Probability of each being a 1 in Test data set.
probabilities_for_test <- data.frame(Id=cancer.test$SiteNum, Predicted=predict(rf_default, newdata=cancer.test, type = "prob")[,2])
write.csv(probabilities_for_test,"/Users/graceedriggs/Documents/STAT 495/Phosphorylations/GD_GBM_probabilities_for_test.csv", row.names= FALSE)

#Probability of each being a 1 in the Train data set
probabilities_for_train <- data.frame(Id=cancer.train$SiteNum, Predicted=predict(rf_default, newdata=cancer.train, type = "prob")[,2])
write.csv(probabilities_for_train,"/Users/graceedriggs/Documents/STAT 495/Phosphorylations/GD_GBM_probabilities_for_train.csv", row.names= FALSE)

## Overall predictions (using TRUE and FALSE instead of 1 and 0)
raw_predictions <- data.frame(Id=cancer.test$SiteNum, Predicted=predict(rf_default, newdata=cancer.test))
raw_predictions <- data.frame(Id=cancer.test$SiteNum, Predicted=(predict(rf_default, newdata=cancer.test)) %>% as.character %>% as.numeric %>% as.logical)
write.csv(raw_predictions,"/Users/graceedriggs/Documents/STAT 495/Phosphorylations/GD_GBM_CancerPredictions.csv", row.names= FALSE)






## attempt 2
############################3Extreme Gradient Boosting
#10 folds repeat 3 times
control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=2)
#Metric compare model is Accuracy
metric <- "Accuracy"
set.seed(123)
rf_default <- train(Response ~ ., 
                    data=cancer.train, 
                    method='xgbTree', 
                    trControl=control)

names(rf_default)
plot(rf_default)
rf_default$bestTune
predict(rf_default, newdata=cancer.test)

#predict
predictions <- data.frame(Id=cancer.test$SiteNum, Predicted=(predict(rf_default, newdata=cancer.test))%>% as.character %>% as.numeric %>% as.logical)
write.csv(predictions,"/Users/graceedriggs/Documents/STAT 495/Phosphorylations/GD_XGB_CancerPredictions.csv", row.names = FALSE)
