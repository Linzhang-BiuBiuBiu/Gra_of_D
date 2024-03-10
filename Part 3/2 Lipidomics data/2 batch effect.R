setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(readr)
library(dplyr)
library(openxlsx)
library(limma)
library(ggplot2)
library(ggpubr)
library(sva)

load("HF_lipid.rdata")
source("E:\\R\\Zl_R_function\\zl_R_Function.r")
#根据QC来校正数据

pos_qc_adj_list <- lapply(pos_lipid_list, qc_normalize)
neg_qc_adj_list <- lapply(neg_lipid_list, qc_normalize)

#将所有原始数据数据合并，如果需要合并QC校正后的，直接用pos_lipid_list1
#合并正负离子原始数据和QC校正后的数据

pos_neg_primry_qcAdj_list <- list(pos_lipid_list,
                             neg_lipid_list,
                             pos_qc_adj_list,
                             neg_qc_adj_list)

pos_neg_primry_qcAdj_list_data <- lapply(pos_neg_primry_qcAdj_list, 
                                         function(x) 
                                             lipid_data_merge(lipid_list=x))

names(pos_neg_primry_qcAdj_list_data) <- c("primary_pos","primary_neg",
                                           "qc_adj_pos","qc_adj_neg")



i <- names(pos_neg_primry_qcAdj_list_data)[2]
#对列表中的数据进行循环过滤
lipid_names <- c(rownames(pos_lipid_list[["pos-5"]]),rownames(neg_lipid_list[["neg-5"]]))

#先分析S1P的变化系数
pos_data_primary <- pos_neg_primry_qcAdj_list_data[[1]]
s1p_data <- pos_data_primary[,c("ID","S-1-P","S-1-P_2","S-1-P_3","S-1-P_4")]
s1p_data[s1p_data==0] <- NA

#计算DCER24：1与SIP的相关性系数
pos_data_primary <- pos_neg_primry_qcAdj_list_data[[1]]
s1p_der_data <- pos_data_primary[,c("ID","S-1-P","S-1-P_2","S-1-P_3","S-1-P_4","DCER(24:1)","DCER(24:1)_2")]
s1p_der_data[,-1] <- mapply(as.numeric, s1p_der_data[,-1])
s1p_der_data[s1p_der_data==0] <- NA
s1p_der_data <- na.omit(s1p_der_data)

cor(s1p_der_data$`DCER(24:1)`,s1p_der_data$`DCER(24:1)_2`)
cor(s1p_der_data$`S-1-P`,s1p_der_data$`DCER(24:1)_2`)

sip_der_coef <- s1p_der_data$`S-1-P`/s1p_der_data$`DCER(24:1)_2`
    
s1p_data_qc_index <- which(strsplit2(s1p_data$ID,"-")[,3]=="QC")
s1p_data_QC <- s1p_data[s1p_data_qc_index,]
s1p_data_QC1 <- s1p_data_QC[!is.na(s1p_data_QC$`S-1-P`)|!is.na(s1p_data_QC$`S-1-P_4`),]

s1p_data_QC_04_both <- s1p_data_QC1[!is.na(s1p_data_QC1$`S-1-P`)&!is.na(s1p_data_QC1$`S-1-P_4`),]
summary(s1p_data_QC_04_both[,-1])
s1p_1_4_coef <- s1p_data_QC_04_both$`S-1-P`/s1p_data_QC_04_both$`S-1-P_4`
summary(s1p_1_4_coef)
mean(s1p_data_QC_04_both$`S-1-P_4`)
mean(s1p_data_QC_04_both$`S-1-P`)

s1p_data_QC_04_na <- s1p_data_QC1[is.na(s1p_data_QC1$`S-1-P`)|is.na(s1p_data_QC1$`S-1-P_4`),]
summary(s1p_data_QC_04_na)
mean(na.omit(s1p_data_QC_04_na$`S-1-P_4`))
mean(na.omit(s1p_data_QC_04_na$`S-1-P`))

s1p_data_QC_02_both <- s1p_data_QC1[!is.na(s1p_data_QC1$`S-1-P`)&!is.na(s1p_data_QC1$`S-1-P_2`),]
summary(s1p_data_QC_02_both[,-1])
mean(s1p_data_QC_02_both$`S-1-P_2`)
mean(s1p_data_QC_02_both$`S-1-P`)

s1p_data_QC_03_both <- s1p_data_QC1[!is.na(s1p_data_QC1$`S-1-P`)&!is.na(s1p_data_QC1$`S-1-P_3`),]
summary(s1p_data_QC_03_both[,-1])
mean(s1p_data_QC_03_both$`S-1-P_3`)
mean(s1p_data_QC_03_both$`S-1-P`)

s1p_data_QC_03_na <- s1p_data_QC1[is.na(s1p_data_QC1$`S-1-P`)|is.na(s1p_data_QC1$`S-1-P_3`),]
summary(s1p_data_QC_02_na)
mean(na.omit(s1p_data_QC_03_na$`S-1-P_3`))
mean(na.omit(s1p_data_QC_03_na$`S-1-P`))
s1p_data[,-1] <- mapply(as.numeric, s1p_data[,-1])
s1p_data$`S-1-P` <- ifelse(is.na(s1p_data$`S-1-P`),s1p_data$`S-1-P_4`/0.7881304,s1p_data$`S-1-P`)
s1p_data$`S-1-P` <- ifelse(is.na(s1p_data$`S-1-P`),s1p_data$`S-1-P_2`/0.5423002,s1p_data$`S-1-P`)

0.5423002
s1p_data_no_0 <- na.omit(s1p_data)

s1p_data_no_0[,-1] <- mapply(as.numeric,s1p_data_no_0[,-1])
cor(s1p_data_no_0[,-1])


s1p_data_no_0_qc <- s1p_data_no_0[which(strsplit2(s1p_data_no_0[,1],"-")[,3]=="QC"),]
cor(s1p_data_no_0_qc[,-1])

s1p_data_QC[is.na(s1p_data_QC)] <- 0
s1p_data_QC[,-1] <- mapply(as.numeric, s1p_data_QC[,-1])


s1p_cof <- data.frame(ID=s1p_data_QC[,1],
                         y_vs_2=s1p_data_QC[,2]/s1p_data_QC[,3],
                         y_vs_3=s1p_data_QC[,2]/s1p_data_QC[,4],
                         y_vs_4=s1p_data_QC[,2]/s1p_data_QC[,5],
                         e_vs_3=s1p_data_QC[,3]/s1p_data_QC[,4],
                         e_vs_4=s1p_data_QC[,3]/s1p_data_QC[,5],
                         s_vs_4=s1p_data_QC[,4]/s1p_data_QC[,5])

s1p_cof[s1p_cof=="Inf"] <- 0
rownames(s1p_cof) <- s1p_cof[,1]
s1p_cof <- s1p_cof[,-1]
s1p_cof <- s1p_cof[rowMeans(s1p_cof)>0,]
count_0 <- apply(s1p_cof, 2, function(x) sum(x==0))
summary(s1p_cof)

lipid_data <- list()
for (i in names(pos_neg_primry_qcAdj_list_data)) {
    tem_data <- pos_neg_primry_qcAdj_list_data[[i]]
    tem_rownames <- tem_data[,1]
    tem_data[,-1] <- apply(tem_data[,-1], 2, as.numeric)
    rownames(tem_data) <- tem_rownames
    tem_data <- tem_data[,-1]
    
    #过滤S-1-P
    if(i=="primary_pos"){
        tem_data[is.na(tem_data)] <- 0
        tem_data$`S-1-P` <- ifelse(tem_data$`S-1-P`!=0,tem_data$`S-1-P`,
                                   ifelse(tem_data$`S-1-P_4`!=0,tem_data$`S-1-P_4`/1.276725354231723,
                                          ifelse(tem_data$`S-1-P_2`!=0,tem_data$`S-1-P_2`*18.3061219945177,tem_data$`S-1-P_3`*0.204336963221697)))
        }
        
    if(i=="qc_adj_pos"){
        tem_data[,"S-1-P"] <- ifelse(tem_data[,"S-1-P"]!=0,ifelse(tem_data[,"S-1-P_2"]>0,
                                                                    tem_data[,"S-1-P_2"],tem_data[,"S-1-P_3"]),
                                     tem_data[,"S-1-P"])}
    
    #提取原始脂质
    tem_data1 <- tem_data[,intersect(colnames(tem_data),lipid_names)]
    lipid_data[[i]] <- tem_data1
    
    #校正批次效应
    #qx=as.numeric(quantile(tem_data1, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    #LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    #if(LogC){
    #   tem_data1[tem_data1<0]=0
    #   tem_data1=log2(tem_data1+2)}
    #tem_data1=normalizeBetweenArrays(tem_data1)
    
    combat_data <- ComBat(t(tem_data1), strsplit2(rownames(tem_data1),"-")[,2],par.prior = TRUE)
    lipid_data[[paste(i,"_comb")]] <- t(combat_data)
    
    #去除第11批
    no_11_data <- tem_data1[-grep("-11-",rownames(tem_data1)),]
    
    lipid_data[[paste(i,"no_11")]] <- no_11_data
}


#画出所有不归一化的PCA图
for(i in names(lipid_data)){
    
    PCA_lipid(data=lipid_data[[i]],
              nor=0,file = paste0("pca_",i,".pdf"),
              main =paste0("pca_",i) )
}

#画出所有归一化的数据
for(i in names(lipid_data)){
    
    PCA_lipid(data=lipid_data[[i]],
              nor=1,file = paste0("nor-PCA_",i,".pdf"),
              main =paste0("nor-PCA_",i) )
}
save(lipid_data,file = "lipid_data.rdata")
load("lipid_data.rdata")

i <- names(lipid_data)[3]

#强制删除QC和KB
lipid_data_noQCKB <- list()
for(i in names(lipid_data)){
    tem_data <- lipid_data[[i]]
    type <- strsplit2(rownames(tem_data),"-")[,3]
    tem_data <- tem_data[-which(type=="QC"|type=="KB"),]
    lipid_data_noQCKB[[i]] <- tem_data
}

#画出删除QC和KB后不归一化的数据
for(i in names(lipid_data_noQCKB)){
    
    PCA_lipid(data=lipid_data_noQCKB[[i]],
              nor=0,file = paste0("qcKb_no_pca_",i,".pdf"),
              main =paste0("qcKb_no_pca_",i) )
}

#画出所有没有QC和KB归一化的数据
for(i in names(lipid_data_noQCKB)){
    
    PCA_lipid(data=lipid_data_noQCKB[[i]],
              nor=1,file = paste0("qcKb_no_pca_",i,"-nor.pdf"),
              main =paste0("qcKb_no_pca_",i,"-nor") )
}


#最后决定用combat,原始数据,QC校正的数据进行分析
#对数据融合后进行分析

for (i in names(lipid_data)) {
    rownames(lipid_data[[i]]) <- gsub("POS-|NEG-","",rownames(lipid_data[[i]]),)
    
}

#合并正负离子模式原始数据
same_primary <- intersect(rownames(lipid_data[["primary_pos"]]),
                           rownames(lipid_data[["primary_neg"]]))

lipid_all_primary <- cbind(lipid_data[["primary_pos"]][same_primary,],
                           lipid_data[["primary_neg"]][same_primary,])

#合并正负离子模式Combat数据
same_combat <- intersect(rownames(lipid_data[["primary_pos _comb"]]),
                          rownames(lipid_data[["primary_neg _comb"]]))

lipid_all_combat <- cbind(lipid_data[["primary_pos _comb"]][same_combat,],
                           lipid_data[["primary_neg _comb"]][same_combat,])

#合并正负离子模式QC校正数据
same_qcadj <- intersect(rownames(lipid_data[["qc_adj_pos"]]),
                         rownames(lipid_data[["qc_adj_neg"]]))

lipid_all_qcadj <- cbind(lipid_data[["qc_adj_pos"]][same_qcadj,],
                          lipid_data[["qc_adj_neg"]][same_qcadj,])



finaly_lipid_list <- list(primary_lipid=lipid_all_primary,
                          combat_lipid=lipid_all_combat,qcadj_lipid=lipid_all_qcadj)


save(finaly_lipid_list,file = "finaly_lipid_list.rdata")

