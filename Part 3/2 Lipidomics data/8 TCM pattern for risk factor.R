setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(readr)
library(dplyr)
library(openxlsx)
library(limma)
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(readr)
library(stringi)
library(stringr)
library(ggrepel)
library(ropls)
library(ggforce)
library(ggprism)
library(pheatmap)
library(pROC)
library(Metrics)
library(randomForest)
library(RSNNS)
library(neuralnet)
library(openxlsx)
library(caret)
library(h2o)
library(NeuralNetTools)
library(mice)
library(VIM)
library(magrittr)
library(PerformanceAnalytics)
library(DT)
library(partykit)
library(klaR)
load("D:\\onedrive\\桌面\\Clinical_data_summary\\临床数据与脂质组的处理\\3_三级索引替换成一级索引\\finaly_clin_lipid_list.rdata")
source("E:\\R\\Zl_R_function\\zl_R_Function.r")
for(i in names(finaly_clin_lipid_list)){
    tem_data <- finaly_clin_lipid_list[[i]]
    tem_data$SEQ <- ifelse(is.na(tem_data$SEQ),tem_data$ID,tem_data$SEQ)
    tem_data$ORGAN <- ifelse(is.na(tem_data$ORGAN),strsplit2(tem_data$SEQ,"-")[,2],tem_data$ORGAN)
    
    finaly_clin_lipid_list[[i]] <- tem_data
}

#准备多重插补的数据
tem_data <- finaly_clin_lipid_list[[4]]
tem_data[is.na(tem_data)] <- 0
#删除空白和QC，固定行数
tem_data <- subset(tem_data,ORGAN!="KB")
tem_data <- subset(tem_data,ORGAN!="QC")

#准备做多重插补，考虑到行名的问题，先把行名保存
#统计NA值，然后删除相关的NA值
tem_data[,c(-1:-7)][tem_data[,c(-1:-7)]==0] <- NA
tem_data <- na_delet_data(tem_data,nrow(tem_data)*0.1)
group_data <- tem_data[,c(1:7)]

tem_data1 <- tem_data[,c(-1:-7)]
tem_data1[is.na(tem_data1)] <- 0

tem_data1<- mapply(as.numeric,tem_data1)

#因为脂质的命名有特殊符号，mice的时候会报错，所以必须保存
colname_tem_data <- colnames(tem_data1)
colnames(tem_data1) <- paste0("name",1:ncol(tem_data1))

#数据归一化才能分析
#注意1和2代表行和列，不要做错
tem_data1 <- t(apply(tem_data1, 1, function(x) x/sum(x)))

#0值替换成NA值，有NA值时归一化会全部归成NA
tem_data1[tem_data1==0] <- NA
#mics多重插补
mice_tem_data1 <- mice(tem_data1,m=5,seed = 5,method = "norm")

#提取多重插补的数据
tem_data <- complete(mice_tem_data1,action = 1)

#恢复脂质名称和分组信息
colnames(tem_data) <- colname_tem_data
tem_data <- cbind(group_data,tem_data)

#只提取第0天的数据
tem_data <- subset(tem_data,TIME =="0")
#只提取C-D期的数据
tem_data <- tem_data[tem_data$分组 %in% c("C","D"),]
#剔除Con，只要R
tem_data <- subset(tem_data,ORGAN =="R")

group <- tem_data$中医证型
clini_group <- tem_data[,c(1:7)]
clini_group$分组 |> filter()


#只保留数据部分
tem_data <- tem_data[,c(-1:-7)]

#提取重要的脂质类型
imp_lip <- read.xlsx("clini_coef_data_RBC.xlsx",colNames = F)
colnames(imp_lip) <- imp_lip[1,]
imp_lip <- imp_lip[-1,]
imp_lip <- colnames(imp_lip)[c(-1:-6,-25)]

tcm_lip <- tem_data[,imp_lip]
tcm_lip <- cbind(group=group,tcm_lip)

#构建训练集和测试集
set.seed(123)
index <-  sort(sample(nrow(tcm_lip), nrow(tcm_lip)*0.9))
train <- tcm_lip[index,]
test <-  tcm_lip[-index,]
#train <- tcm_lip

#随机森林
train_data <- as.data.frame(cbind(group=as.factor(train[,1]),train[,-1]))
train_data[,-1] <- sapply(train_data[,-1],as.numeric)
train_data$group <- as.factor(train_data$group)

test_data <- as.data.frame(cbind(group=as.factor(test[,1]),test[,-1]))
test_data[,-1] <- sapply(test_data[,-1],as.numeric)
test_data$group <- as.factor(test_data$group)

testm_data <- as.data.frame(cbind(group=as.factor(test_m[,1]),test_m[,-1]))
testm_data[,-1] <- sapply(testm_data[,-1],as.numeric)
testm_data$group <- as.factor(testm_data$group)

set.seed(2)
rftune <- tuneRF(x=train_data[,-1],y=train_data[,1],stepFactor = 1,ntreeTry = 500,nflod=10)
rftune_min_OBB <- min(which(grepl(min(rftune[,2]),rftune[,2])))
rftune_optimal <- rftune[rftune_min_OBB,1]

RandomForest <- randomForest(x=train_data[,-1], 
                             y=train_data[,1],
                             ntree=500,
                             mtry=rftune_optimal,
                             nflod=10)

rf_train_data <- predict(RandomForest,train_data[,-1])
rf_train_accuracy <- accuracy(as.vector(rf_train_data),train_data[,1])
rf_train_confu <- confusionMatrix(RandomForest$predicted,train_data[,1])

rf_test <- predict(RandomForest,test_data[,-1])
rf_test_accuracy <- accuracy(as.vector(rf_test),test_data$group)
rf_test_confu <- confusionMatrix(rf_test,test_data$group)

rf_testm <- predict(RandomForest,testm_data[,-1])
rf_testm_accuracy <- accuracy(as.vector(rf_testm),testm_data$group)
rf_testm_confu <- confusionMatrix(rf_testm,testm_data$group)

#随机森林的图
pdf(file="2_RandomForest.pdf",width = 12,height = 9)
plot(RandomForest)
varImpPlot(RandomForest,sort=T,main = "The importance of significent genes",color = c("blue"))
dev.off()

#统计随机森林的脂质系数
RF_coef_tem=as.matrix(RandomForest$importance[-1,])
colnames(RF_coef_tem) <- "RF"
RF_coef <- cbind(ID=rownames(RF_coef_tem),RF_coef_tem)

#SVM支持向量机
library(e1071)
#调试支持向量机
set.seed(123)
svm_train = tune.svm(group~.,
                     data=train_data,
                     kernel="radial",
                     cost=1:10,
                     fold=10)

#选择支持向量机的最佳模型
svm_train_again <- predict(svm_train$best.model,train_data[,-1],type = 'response')
svm_train_accuracy <- accuracy(svm_train_again,train_data$group)
svm_train_confu <- confusionMatrix(svm_train_again, train_data$group) 

svm_test <- predict(svm_train$best.model,test_data[,-1],type = 'response')
svm_test_accuracy <- accuracy(svm_test,test_data$group)
svm_test_confu <- confusionMatrix(svm_test, test_data$group) 

svm_testm <- predict(svm_train$best.model,testm_data[,-1],type = 'response')
svm_testm_accuracy <- accuracy(svm_testm,testm_data$group)
svm_test_confu <- confusionMatrix(svm_test, test_data$group) 

#GBM
#梯度提升机
Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk1.8.0_341")
h2o.init(nthreads = 8,max_mem_size = "8G")

#数据转换为h2o格式
gbm_model_data <- as.h2o(train_data)
gbm_model_test <- as.h2o(test_data)

target <- "group"
predictors <- colnames(train_data)[-1]

#构建模型
gbm_grid_cla <- h2o.gbm(x=predictors,
                        y=target,
                        distribution="AUTO",
                        training_frame = gbm_model_data,
                        ntree=100,
                        learn_rate=0.01,
                        sample_rate = 0.8,
                        col_sample_rate = 0.6,
                        seed=1234)

aml <- h2o.automl(y = target,
                  training_frame = gbm_model_data,
                  max_models = 6,
                  seed = 1)

#GBM作图
pdf(file="4_GBM_Variable_importance.PDF",width = 10,height = 6)
h2o.varimp_heatmap(aml)
h2o.varimp_plot(gbm_grid_cla,round(ncol(train_data)/2,0))
dev.off()

#判断测试集和验证机的准确性
#GBM还需要重新整理
gbm_train_data <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_data[,-1]))
gbm_train_accuracy <- accuracy(as.vector(gbm_model_data$group),gbm_train_data$predict)
gbm_train_confu <- confusionMatrix(gbm_train_data$predict,train_data$group)

gbm_test <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_test[,-1]))
gbm_test_accuracy <- accuracy(as.vector(gbm_model_test$group),gbm_test$predict)
gbm_test_confu <- confusionMatrix(gbm_test$predict,test_data$group)

#寻找GBM中脂质的重要性
GBM_result <- summary(gbm_grid_cla)
gbm_coef <- GBM_result[,c(1,3)]
rownames(gbm_coef) <- gbm_coef[,1]
colnames(gbm_coef) <- c("ID","GBM")

#BP神经网络
set.seed(234)
colnames_train_data <- colnames(train_data)
colnames(train_data) <- paste0("name",1:ncol(train_data))
colnames(test_data) <- paste0("name",1:ncol(train_data))

n <- colnames(train_data)

form <- as.formula(paste("group~",paste(n,collapse = "+")))
train_data[,1:ncol(train_data)] <- lapply(train_data[,1:ncol(train_data)], as.numeric)
test_data[,1:ncol(test_data)] <- lapply(test_data[,1:ncol(test_data)],as.numeric)

#BP设置独热变量
train_data1 <- train_data
test_data1 <- test_data

train_data1$one <- train_data[,1] ==1
train_data1$two <- train_data[,1] ==2
train_data1$three <- train_data[,1] ==3

bp_train <- neuralnet(one+two+three~.,
                      data=train_data1[,-1], 
                      threshold=0.01,
                      stepmax=1000000000,
                      err.fct="ce",
                      linear.output=F,
                      hidden=5)

#dev.new()
#plot(BP_train)

#模型结果
bp_train_pred_tem <- as.data.frame(bp_train$net.result)
bp_train_pred =c("1","2","3")[apply(bp_train_pred_tem, 1, which.max)] %>%  as.factor()

bp_train_confu <- confusionMatrix(as.factor(bp_train_pred),as.factor(train_data[,1]))
bp_train_accuracy <- accuracy(bp_train_pred,train_data[,1])

#测试验证集的结果
test_data1$one <- test_data1[,1] ==1
test_data1$two <- test_data1[,1] ==2
test_data1$three <- test_data1[,1] ==3
bp_test_result <- compute(bp_train,test_data1[,-1])

#查看混淆矩阵

bp_test_pred =c("1","2","3")[apply(bp_test_result$net.result, 1, which.max)] %>%  as.factor()
predict.table = table(test_data1[,1],bp_test_pred)
bp_test_confu <- confusionMatrix(bp_test_pred,as.factor(test_data1[,1]))
bp_test_accuracy <- accuracy(bp_test_pred,test_data1[,1])

bp_test_roc <- multiclass.roc(test_data1[,1],as.numeric(bp_test_pred),direction = "<")
bp_train_roc <- multiclass.roc(train_data[,1],as.numeric(bp_train_pred),direction = "<")

#nn的图
pdf(file="6_Neuralnet.pdf",width=7,height = 20)
par(cex = 1)
plotnet(bp_train,max_sp=1,pad_x=0.7,prune_lty=0.5,circle_cex=3,pos_col = "red", neg_col = "grey")
dev.off()

#将NN中train和test的列名转换回来
colnames(train_data) <- colnames_train_data
colnames(test_data) <- colnames_train_data

#决策树模型
library(rpart)
library(rpart.plot)
train_data <- as.data.frame(train_data)
rpart_model <- rpart(group~.,data=train_data, method = "class",cp=0.000001)

#作图
pdf(file="5_Decision trees_no_filter.pdf",width = 15,height = 15)
rpart.plot(rpart_model,type=3,extra = "auto",
           under=T,fallen.leaves=F,cex=0.7,main="Decision trees")
plotcp(rpart_model)
dev.off()

#剪枝
bestcp <- rpart_model$cptable[which.min(rpart_model$cptable[,"xerror"]),"CP"]
rpart_model.pruned <- prune(rpart_model,cp=bestcp)
par(family="STKaiti")
#剪枝后的图
pdf(file="5_Decision trees_filter.pdf",width = 15,height = 15)
rpart.plot(rpart_model.pruned,type=2,extra = "auto",
           under=T,fallen.leaves=F,cex=2,main="Optimal Decision trees")
dev.off()

#训练集和验证集的正确率
dt_train_data <- predict(rpart_model,train_data[,-1])
dt_train_pre <- ifelse(dt_train_data[,1]>dt_train_data[,2]&dt_train_data[,1]>dt_train_data[,3],1,ifelse(dt_train_data[,2]>dt_train_data[,3],2,3))
dt_train_accuracy <- accuracy(as.vector(dt_train_pre),train_data[,1])
dt_train_confu <- confusionMatrix(as.factor(dt_train_pre),train_data[,1])

dt_test_data <- predict(rpart_model,test_data[,-1])
dt_test_pre <- ifelse(dt_test_data[,1]>dt_test_data[,2]&dt_test_data[,1]>dt_test_data[,3],1,ifelse(dt_test_data[,2]>dt_test_data[,3],2,3))
dt_test_accuracy <- accuracy(as.vector(dt_test_pre),test_data[,1])
dt_test_confu <- confusionMatrix(as.factor(dt_test_pre),test_data[,1])

#提取DT的系数
dt_coef_tem <- as.matrix(rpart_model.pruned$variable.importance)
dt_coef <- cbind(rownames(dt_coef_tem),dt_coef_tem)
colnames(dt_coef) <- c("ID","DT")

#总结所有的coef
rpart_roc_test <- multiclass.roc(test_data[,1],as.numeric(dt_test_pre),direction = "<") 
rf_roc_test <- multiclass.roc(test_data[,1], as.numeric(rf_test) ,direction = "<") 
gbm_roc_test <- multiclass.roc(test_data[,1], as.numeric(gbm_test$predict),direction = "<") 
svm_roc_test <- multiclass.roc(test_data[,1], as.numeric(svm_test),direction = "<")
bp_roc_test <- multiclass.roc(test_data[,1],as.numeric(bp_test_pred),direction = "<")

rpart_roc_train <- multiclass.roc(train_data[,1],as.numeric(dt_train_pre),direction = "<") 
rf_roc_train <- multiclass.roc(train_data[,1], as.numeric(rf_train_data) ,direction = "<") 
gbm_roc_train <- multiclass.roc(train_data[,1], as.numeric(gbm_train_data$predict),direction = "<") 
svm_roc_train <- multiclass.roc(train_data[,1], as.numeric(svm_train_again),direction = "<")
bp_roc_train <- multiclass.roc(train_data[,1],as.numeric(bp_train_pred),direction = "<")


train_test_accuracy <- data.frame(train="train",test="test",AUCTrain="AUCTrain",AUCTest="AUCTest")
train_test_accuracy <- rbind(train_test_accuracy,c(rf_train_accuracy,rf_test_accuracy,rf_roc_train$auc,rf_roc_test$auc))
train_test_accuracy <- rbind(train_test_accuracy,c(svm_train_accuracy,svm_test_accuracy,svm_roc_train$auc,svm_roc_test$auc))[-1,]
train_test_accuracy <- rbind(train_test_accuracy,c(bp_train_accuracy,bp_test_accuracy,bp_roc_train$auc,bp_roc_test$auc))
train_test_accuracy <- rbind(train_test_accuracy,c(dt_train_accuracy,dt_test_accuracy,rpart_roc_train$auc,rpart_roc_test$auc))
train_test_accuracy <- rbind(train_test_accuracy,c(gbm_train_accuracy,gbm_test_accuracy,gbm_roc_train$auc,gbm_roc_test$auc))

rownames(train_test_accuracy) <- c("RF","SVM","BP","DT","GBM")

write.xlsx(cbind(rownames(train_test_accuracy),train_test_accuracy),file = "train_test_accuracy.xlsx")
