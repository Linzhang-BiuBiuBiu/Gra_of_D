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
load("finaly_clin_lipid_list.rdata")
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
save(mice_tem_data1,file = "mice_tem_data1.rdata")

#提取多重插补的数据
tem_data <- complete(mice_tem_data1,action = 1)

#恢复脂质名称和分组信息
colnames(tem_data) <- colname_tem_data
tem_data <- cbind(group_data,tem_data)

#只提取第0天的数据
tem_data <- subset(tem_data,TIME =="0")
#A-B,C-D统一分期
tem_data$分组 <- ifelse(tem_data$分组=="Con","Con",
                      ifelse(tem_data$分组=="A","AB",
                             ifelse(tem_data$分组=="B","AB",
                                    ifelse(tem_data$分组=="C","CD","CD"))))

#剔除Con，只要R
#tem_data <- subset(tem_data,分组 =="B"|分组 =="D")
tem_data <- subset(tem_data,分组 !="Con")
tem_data <- subset(tem_data,ORGAN =="R")

group <- tem_data$分组
clini_group <- tem_data[,c(1:7)]
#只保留数据部分
tem_data <- tem_data[,c(-1:-7)]

#opls-da分析
pdf("A_B-C_D_oplsda_diff.pdf")                   
df1_oplsda <- opls(tem_data, group, predI = 1,orthoI = NA)                 
dev.off() 

#提取VIP值
data_VIP <- as.data.frame(df1_oplsda@vipVn)
colnames(data_VIP) <- "VIP"

#差异分析
i <- unique(group)[2]
j <- unique(group)[1]
design <- model.matrix(~0+factor(group))
colnames(design) <- c(i,j)
fit <- lmFit(t(tem_data),design)
cont.matrix<-makeContrasts(i-j,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

diff <- topTable(fit2,adjust='fdr',number=200000)[rownames(data_VIP),]

iave <- as.data.frame(apply(tem_data[which(group==i),], 2, mean))
jave <- as.data.frame(apply(tem_data[which(group==j),], 2, mean))

FC_LOG <- data.frame("FC"=iave/jave,"log2"=log2(iave)-log2(jave))
colnames(FC_LOG) <- c("FC","Log2FC")

allDiff1 <- cbind(diff,FC_LOG)
allDiff1 <- cbind(allDiff1,data_VIP)

#差异的筛选条件,首先找出差异脂质，做出方便后续统计的表格
allDiff2 <- allDiff1[as.factor(allDiff1$FC) != "Inf",]
allDiff2 <- allDiff2[as.factor(allDiff2$FC) != "0",]
allDiff2 <- allDiff2[!is.nan(allDiff2$Log2FC),]

#设置筛选条件
sig_diff <- allDiff2[allDiff2$adj.P.Val <0.05,]
sig_diff <- sig_diff[order(as.numeric(as.vector(abs(sig_diff$VIP))),decreasing = T),]

#筛选出机器学习的差异脂质
dif_li <- tem_data[,colnames(tem_data) %in% rownames(sig_diff)]
dif_lip <- cbind(group=group,dif_li)

#构建训练集和测试集
#统计这些脂质的ROC
#load("train_test.rdata")
dif_lip[,-1] <- mapply(as.numeric, dif_lip[,-1])
roc_all <- apply(dif_lip[,-1], 2, function(x) roc(as.factor(dif_lip[,1]),as.numeric(x)))
roc_all1 <- sapply(roc_all, function(x) x$auc)
roc_all2 <- cbind(ID=rownames(roc_all1),roc_all1)
roc_all2 <- cbind(ID=rownames(roc_all2),roc_all2)
roc_all2 <- roc_all2[order(roc_all2[,2],decreasing = T),]

#根据排名，画出前5个脂质的AUC
for(i in 1:5){
    j <- rownames(roc_all2)[i]
    roc1 <- roc_all[[j]]
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=round(as.numeric(ci1),3)
    
    pdf(file=paste0("ROC_top_",i,".pdf"), width=5, height=5)
    plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=j)
    text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"−",sprintf("%.03f",ciVec[3])), col="red")
    dev.off()
}

#广义线性模型寻找差异脂质
dif_lip[,-1] <- apply(dif_lip[,-1], 2, as.numeric)

dif_lip[,1] <- as.factor(dif_lip[,1])
model_1<-glm(group~.,data = dif_lip, family = binomial(link ="logit"))

#统计模型中的差异脂质，用差异脂质再做一次glm
sumar_glm <- summary(model_1)
glm_coef <- sumar_glm$coefficients[-1,]
glm_sig_lip <- glm_coef[glm_coef[,4] < 0.05,]
glm_sig_lip <- glm_sig_lip[order(glm_sig_lip[,4],decreasing = F),]
glm_sig_lip <- gsub("\\`","",rownames(glm_sig_lip))

#广义线性模型again
model_2<-glm(group~.,data = dif_lip[,c("group",glm_sig_lip)], family = binomial(link ="logit"))

#tbl做出相关的表，但是由于变量太多，可能会影响到运行速度
sumar_glm <- summary(model_2)
glm_coef <- sumar_glm$coefficients[-1,]
glm_sig_lip <- glm_coef[glm_coef[,4] < 0.05,]
glm_sig_lip <- glm_sig_lip[order(glm_sig_lip[,4],decreasing = F),]
glm_sig_lip <- gsub("\\`","",rownames(glm_sig_lip))

#广义线性模型again2
model_3<-glm(group~.,data = dif_lip[,c("group",glm_sig_lip)], family = binomial(link ="logit"))

#tbl做出相关的表，但是由于变量太多，可能会影响到运行速度
sumar_glm <- summary(model_3)
glm_coef <- sumar_glm$coefficients[-1,]
glm_sig_lip <- glm_coef[glm_coef[,4] < 0.00000001,]
glm_sig_lip <- glm_sig_lip[order(glm_sig_lip[,4],decreasing = F),]
glm_sig_lip <- gsub("\\`","",rownames(glm_sig_lip))
glm_sig_lip <- c("DCER(24:1)","HexCer d18:1/12:0")

#glm的again4,用预测出的值，来做ROC
model_4<-glm(group~.,data = dif_lip[,c("group",glm_sig_lip)], family = binomial(link ="logit"))
sumar_glm <- summary(model_4)
fitted.prob<-predict(model_4, newdata = dif_lip, type = "response")
roc_multivar_3<-roc(dif_lip$group,model_3$fitted.values)

#做出最终的glm合并的图
pdf("glm_dif_ROC.pdf",width = 5,height = 5)
plot.roc(roc_multivar_3,print.auc=TRUE, col="blue",legacy.axes=T, main=paste(glm_sig_lip,collapse = "+")) 
ci1=ci.auc(roc_multivar_3, method="bootstrap")
ciVec=round(as.numeric(ci1),3)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"−",sprintf("%.03f",ciVec[3])), col="red")
dev.off()

#分组整理成factor
dif_lip$group <- ifelse(dif_lip$group=="AB",0,1)
index <-  sort(sample(nrow(dif_lip), nrow(dif_lip)*0.9))
train <- dif_lip[index,]
test <-  dif_lip[-index,]
#train <- dif_lip

#LASSO筛选重要脂质
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)

x=as.matrix(train[,-1])
y=train[,1]

cvfit=cv.glmnet(x, y, family="binomial", 
                alpha=1,
                type.measure='deviance',
                nfolds = 10)

#lasso筛选的重要变量
pdf(file="1_lasso_cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#统计lasso筛选出来的系数和重要脂质
lasso_coef=as.matrix(coef(cvfit, s = cvfit$lambda.1se)[-1,])
colnames(lasso_coef) <- "LASSO"
lasso_coef <- cbind(ID=rownames(lasso_coef),lasso_coef)

index=which(lasso_coef != 0)
lassoLipid=row.names(lasso_coef)[index]
lassoLipid=lassoLipid[-1]

#随机森林
train_data <- as.data.frame(cbind(group=as.factor(train[,1]),train[,-1]))
train_data[,-1] <- sapply(train_data[,-1],as.numeric)
train_data$group <- as.factor(train_data$group)

test_data <- as.data.frame(cbind(group=as.factor(test[,1]),test[,-1]))
test_data[,-1] <- sapply(test_data[,-1],as.numeric)
test_data$group <- as.factor(test_data$group)

x <- train_data[,-1]
y <- train_data[,1]

set.seed(2)
rftune <- tuneRF(x=x,y=y,stepFactor = 1,ntreeTry = 500,nflod=10)
rftune_min_OBB <- min(which(grepl(min(rftune[,2]),rftune[,2])))
rftune_optimal <- rftune[rftune_min_OBB,1]

RandomForest <- randomForest(x=x, 
                             y=y,
                             ntree=500,
                             mtry=rftune_optimal,
                             nflod=10)

rf_train_data <- predict(RandomForest,train_data)
rf_train_accuracy <- accuracy(as.vector(rf_train_data),y)

rf_test <- predict(RandomForest,test_data)
rf_test_accuracy <- accuracy(as.vector(rf_test),test_data$group)

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
                     kernel="polynomial",
                     cost=1:10,
                     fold=10)

#选择支持向量机的最佳模型
svm_train_again <- as.character(predict(svm_train$best.model,train_data))
svm_train_accuracy <- accuracy(svm_train_again,train_data$group)

svm_test <- as.character(predict(svm_train$best.model,test_data))
svm_test_accuracy <- accuracy(svm_test,test_data$group)

#寻找SVM的脂质重要性系数

Profile <- rfe(x=train_data[,-1],
            y=train_data[,1],
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")

#????ͼ??
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#??ע??????֤??????С?ĵ?
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

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
gbm_train_data <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_data))
gbm_train_accuracy <- accuracy(as.vector(gbm_model_data$group),gbm_train_data$predict)

gbm_test <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_test))
gbm_test_accuracy <- accuracy(as.vector(gbm_model_test$group),gbm_test$predict)

#寻找GBM中脂质的重要性
GBM_result <- summary(gbm_grid_cla)
gbm_coef <- GBM_result[,c(1,3)]
rownames(gbm_coef) <- gbm_coef[,1]
colnames(gbm_coef) <- c("ID","GBM")

#神经网络
set.seed(234)
colnames_train_data <- colnames(train_data)
colnames(train_data) <- paste0("name",1:ncol(train_data))
colnames(test_data) <- paste0("name",1:ncol(train_data))

n <- colnames(train_data)

form <- as.formula(paste("group~",paste(n,collapse = "+")))

train_data[,1:ncol(train_data)] <- lapply(train_data[,1:ncol(train_data)], as.numeric)

test_data[,1:ncol(test_data)] <- lapply(test_data[,1:ncol(test_data)],
                                          as.numeric)

nn_train <- neuralnet(name1~.,data=train_data,
                      hidden = c(1,2,3),
                      err.fct = "sse",
                      linear.output = T)

nn_train_data <- predict(nn_train,train_data)
nn_train_data <- ifelse(nn_train_data>1.5,2,1)
nn_train_accuracy <- accuracy(as.vector(nn_train_data),train_data[,1])

nn_test <- predict(nn_train,test_data)
nn_test <- ifelse(nn_test>1.5,2,1)
nn_test_accuracy <- accuracy(as.vector(nn_test),test_data[,1])

#nn的图
pdf(file="6_Neuralnet.pdf",width=7,height = 20)
par(cex = 1)
plotnet(nn_train,max_sp=1,pad_x=0.7,prune_lty=0.5,circle_cex=3,
        pos_col = "red", neg_col = "grey")
dev.off()

#nn的系数
nn_coef <- nn_train$result.matrix
nn_coef <- data.frame(nn_coef[grep("^name",rownames(nn_coef)),])
rownames(nn_coef) <- strsplit2(rownames(nn_coef),"\\.")[,1]
colnames(nn_coef) <- "NN"
nn_name_index <- data.frame(colnames_train_data[-1],colnames(train_data)[-1])
rownames(nn_name_index) <- nn_name_index[,2]
nn_coef <- cbind(nn_name_index,nn_coef)[,-2]
rownames(nn_coef) <- nn_coef[,1]
colnames(nn_coef)[1] <- "ID"

#将NN中train和test的列名转换回来
colnames(train_data) <- colnames_train_data
colnames(test_data) <- colnames_train_data

#决策树模型
library(rpart)
library(rpart.plot)
train_data <- as.data.frame(train_data)
rpart_model <- rpart(group~.,
                     data=train_data,
                     method = "class",
                     cp=0.000001)

#作图
pdf(file="5_Decision trees_no_filter.pdf",width = 15,height = 15)
rpart.plot(rpart_model,type=2,extra = "auto",
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
dt_train_data <- predict(rpart_model,train_data)
dt_train_pre <- ifelse(dt_train_data[,1]>=dt_train_data[,2],1,2)
dt_train_accuracy <- accuracy(as.vector(dt_train_pre),train_data[,1])

dt_test_data <- predict(rpart_model,test_data)
dt_test_pre <- ifelse(dt_test_data[,1]>=dt_test_data[,2],1,2)
dt_test_accuracy <- accuracy(as.vector(dt_test_pre),test_data[,1])

#提取DT的系数
dt_coef_tem <- as.matrix(rpart_model.pruned$variable.importance)
dt_coef <- cbind(rownames(dt_coef_tem),dt_coef_tem)
colnames(dt_coef) <- c("ID","DT")

#总结所有的coef
train_test_accuracy <- data.frame(train="train",test="test")
train_test_accuracy <- rbind(train_test_accuracy,c(rf_train_accuracy,rf_test_accuracy))
train_test_accuracy <- rbind(train_test_accuracy,c(svm_train_accuracy,svm_test_accuracy))[-1,]
train_test_accuracy <- rbind(train_test_accuracy,c(nn_train_accuracy,nn_test_accuracy))
train_test_accuracy <- rbind(train_test_accuracy,c(dt_train_accuracy,dt_test_accuracy))
train_test_accuracy <- rbind(train_test_accuracy,c(gbm_train_accuracy,gbm_test_accuracy))

rownames(train_test_accuracy) <- c("RF","SVM","NN","DT","GBM")

#opls-da的coef
oplsda_coef <- as.matrix(cbind(ID=rownames(sig_diff),sig_diff[,c("VIP","Log2FC")]))
colnames(oplsda_coef)[2] <- "OPLSDA"

all_ML_coef <- merge(lasso_coef,nn_coef,all=T)
all_ML_coef <- merge(all_ML_coef,gbm_coef,all=T)
all_ML_coef <- merge(all_ML_coef,RF_coef,all=T)
all_ML_coef <- merge(all_ML_coef,dt_coef,all = T)
all_ML_coef <- merge(all_ML_coef,oplsda_coef,all=T)
all_ML_coef[is.na(all_ML_coef)] <- 0

all_ML_coef[,-1] <- mapply(as.numeric, all_ML_coef[,-1])
all_ML_coef[,c(-1,-8)] <- apply(all_ML_coef[,c(-1,-8)], 2, function(x) abs(x)/max(abs(x)))
all_ML_coef[all_ML_coef=="NaN"] <- 0
all_ML_coef$all_wei <- apply(all_ML_coef[,c(-1,-8)],1,function(x) sum(abs(x)))
all_ML_coef$sum_coef <- ifelse(all_ML_coef$Log2FC>0,-all_ML_coef$all_wei,all_ML_coef$all_wei)
all_ML_coef <- all_ML_coef[order(all_ML_coef$all_wei,decreasing = T),]

all_ML_coef <- all_ML_coef[order(all_ML_coef$OPLSDA,decreasing = T),]
 #all_ML_coef_list <- list()
#all_ML_coef_list[[5]] <- all_ML_coef

#train_test_accuracy_list <- list()
#train_test_accuracy_list[[5]] <- train_test_accuracy

#write.xlsx(all_ML_coef_list,file = "all_ML_coef_list.xlsx")
#write.xlsx(train_test_accuracy_list,file = "train_test_accuracy_list.xlsx")

#根据ROC筛选具有诊断意义的脂质
roc_lipid <- c("DCER(24:1)", 
               "DCER(22:0)", "HexCer.d18:1/12:0", "DCER(24:0)", "PE(O-18:0/20:4)", 
               "PE(P-18:0/20:4)", "HCER(16:0)", "PE(O-18:0/22:5)", "PE(P-18:0/22:4)", 
               "PE(O-18:0/22:4)", "HCER(24:1)", "PE(16:0/18:2)", "PE(P-18:0/22:5)", 
               "PE(16:0/20:2)")

#计算新的数据分级
#根据总权重计算coef,只留下权重大于1的数据
all_ML_coef1 <- all_ML_coef[all_ML_coef$ID %in% roc_lipid,]
all_ML_coef1$sum_coef <- ifelse(all_ML_coef1$Log2FC < 0,all_ML_coef1$all_wei,-all_ML_coef1$all_wei)

coef_lipi <- dif_lip[,all_ML_coef1$ID]
coef_lipi$coef <- apply(coef_lipi,1,function(x) sum(x*all_ML_coef1$sum_coef))
clini_coef_data <- cbind(group,coef_lipi$coef)

#重新载入AB分期的数据
tem_data <- complete(mice_tem_data1,action = 1)
#恢复脂质名称和分组信息
colnames(tem_data) <- colname_tem_data
tem_data <- cbind(group_data,tem_data)

#只提取第0天的数据
tem_data <- subset(tem_data,TIME =="0")
tem_data <- subset(tem_data,ORGAN=="R")

coef_lipi <- tem_data[,all_ML_coef1$ID]
coef_lipi$coef <- apply(coef_lipi,1,function(x) sum(x*all_ML_coef1$all_wei))

coef_lipi[,-13] <- apply(coef_lipi[,-13], 2, function(x) x+(2*abs(min(x))))
coef_lipi <- apply(coef_lipi, 2, function(x) x/min(x))

clini_coef_data <- cbind(tem_data[,1:6],coef_lipi)

write.xlsx(clini_coef_data,file = "clini_coef_data.xlsx")

#根据筛选出来的重要脂质，重新做分类模型
#分组整理成factor
dif_lip2 <- dif_lip[,c("group",roc_lipid)]

#opls-da分析
pdf("A_B-C_D_oplsda_v2.pdf")                   
df1_oplsda <- opls(dif_lip2, group, predI = 1,orthoI = NA)                 
dev.off() 

index <-  sort(sample(1:nrow(dif_lip2), nrow(dif_lip2)*0.9))
train <- dif_lip2[index,]
test <-  dif_lip2[-index,]

#随机森林
train_data <- as.data.frame(cbind(group=as.factor(train[,1]),train[,-1]))
train_data[,-1] <- sapply(train_data[,-1],as.numeric)
train_data$group <- as.factor(train_data$group)

test_data <- as.data.frame(cbind(group=as.factor(test[,1]),test[,-1]))
test_data[,-1] <- sapply(test_data[,-1],as.numeric)
test_data$group <- as.factor(test_data$group)

x <- train_data[,-1]
y <- train_data[,1]

set.seed(2)
rftune <- tuneRF(x=x,y=y,stepFactor = 1,ntreeTry = 500,nflod=10)
rftune_min_OBB <- min(which(grepl(min(rftune[,2]),rftune[,2])))
rftune_optimal <- rftune[rftune_min_OBB,1]

RandomForest <- randomForest(x=x, 
                             y=y,
                             ntree=500,
                             mtry=rftune_optimal,
                             nflod=10)

rf_train_data <- predict(RandomForest,train_data)
rf_train_accuracy <- accuracy(as.vector(rf_train_data),y)

rf_test <- predict(RandomForest,test_data)
rf_test_accuracy <- accuracy(as.vector(rf_test),test_data$group)

#SVM支持向量机
library(e1071)
#调试支持向量机
set.seed(123)
svm_train = tune.svm(group~.,
                     data=train_data,
                     kernel="polynomial",
                     cost=1:10,
                     fold=10)

#选择支持向量机的最佳模型
svm_train_again <- as.character(predict(svm_train$best.model,train_data))
svm_train_accuracy <- accuracy(svm_train_again,train_data$group)

svm_test <- as.character(predict(svm_train$best.model,test_data))
svm_test_accuracy <- accuracy(svm_test,test_data$group)

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
gbm_train_data <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_data))
gbm_train_accuracy <- accuracy(as.vector(gbm_model_data$group),gbm_train_data$predict)

gbm_test <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_test))
gbm_test_accuracy <- accuracy(as.vector(gbm_model_test$group),gbm_test$predict)

#寻找GBM中脂质的重要性
GBM_result <- summary(gbm_grid_cla)
gbm_coef <- GBM_result[,c(1,3)]
rownames(gbm_coef) <- gbm_coef[,1]
colnames(gbm_coef) <- c("ID","GBM")

#神经网络
set.seed(234)
colnames_train_data <- colnames(train_data)
colnames(train_data) <- paste0("name",1:ncol(train_data))
colnames(test_data) <- paste0("name",1:ncol(train_data))

n <- colnames(train_data)

form <- as.formula(paste("group~",paste(n,collapse = "+")))

train_data[,1:ncol(train_data)] <- lapply(train_data[,1:ncol(train_data)], as.numeric)

test_data[,1:ncol(test_data)] <- lapply(test_data[,1:ncol(test_data)],
                                        as.numeric)

nn_train <- neuralnet(name1~.,data=train_data,
                      hidden = c(1,2,3),
                      err.fct = "sse",
                      linear.output = T)

nn_train_data <- predict(nn_train,train_data)
nn_train_data <- ifelse(nn_train_data>1.5,2,1)
nn_train_accuracy <- accuracy(as.vector(nn_train_data),train_data[,1])

nn_test <- predict(nn_train,test_data)
nn_test <- ifelse(nn_test>1.5,2,1)
nn_test_accuracy <- accuracy(as.vector(nn_test),test_data[,1])

#将NN中train和test的列名转换回来
colnames(train_data) <- colnames_train_data
colnames(test_data) <- colnames_train_data

#决策树模型
library(rpart)
library(rpart.plot)
train_data <- as.data.frame(train_data)
rpart_model <- rpart(group~.,
                     data=train_data,
                     method = "class",
                     cp=0.000001)

#剪枝
bestcp <- rpart_model$cptable[which.min(rpart_model$cptable[,"xerror"]),"CP"]
rpart_model.pruned <- prune(rpart_model,cp=bestcp)
par(family="STKaiti")

#训练集和验证集的正确率
dt_train_data <- predict(rpart_model,train_data)
dt_train_pre <- ifelse(dt_train_data[,1]>=dt_train_data[,2],1,2)
dt_train_accuracy <- accuracy(as.vector(dt_train_pre),train_data[,1])

dt_test_data <- predict(rpart_model,test_data)
dt_test_pre <- ifelse(dt_test_data[,1]>=dt_test_data[,2],1,2)
dt_test_accuracy <- accuracy(as.vector(dt_test_pre),test_data[,1])
save.image(file = "tem_all_result.RData")

#用opls-da做出的重要性评分
sig_diff$coef <- ifelse(sig_diff$Log2FC <0, sig_diff$VIP,-sig_diff$VIP)
sig_diff2 <- sig_diff[sig_diff$VIP >1.5,]

#重新载入AB分期的数据
tem_data <- complete(mice_tem_data1,action = 1)
#恢复脂质名称和分组信息
colnames(tem_data) <- colname_tem_data
tem_data <- cbind(group_data,tem_data)

#只提取第0天的数据
tem_data <- subset(tem_data,TIME =="0")
tem_data <- subset(tem_data,ORGAN=="R")

coef_lipi <- tem_data[,rownames(sig_diff2)]
coef_lipi$coef <- apply(coef_lipi,1,function(x) sum(x*sig_diff2$coef))

#血浆
coef_lipi[,c(-7,-14)] <- apply(coef_lipi[,c(-7,-14)], 2, function(x) x+(2*abs(min(x))))
#红细胞
coef_lipi[,c(-17)] <- apply(coef_lipi[,c(-17)], 2, function(x) x+(2*abs(min(x))))

coef_lipi <- apply(coef_lipi, 2, function(x) x/min(x))

clini_coef_data <- cbind(tem_data[,1:6],coef_lipi)

write.xlsx(clini_coef_data,file = "clini_coef_data.xlsx")

