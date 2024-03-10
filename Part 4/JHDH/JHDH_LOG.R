library(readxl)
library(openxlsx)
library(stringr)
library(dplyr)
library(nhanesR)
library(limma)
library(pheatmap)
library(tidyverse)
library(MatchIt)#倾向性评分控制基线
library(mice)
library(gtsummary)          
library(flextable)
source("E:\\R\\Zl_R_function\\zl_R_Function_v2.r")

#设置工作目录为R脚本所在的路径
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(digits = 15)
load("analysis_datav3_RMR.RData")
#analysis_data <- read.xlsx("analysis_datav3.xlsx",1)

#analysis_data <- analysis_data[!(analysis_data$group==0&analysis_data$B型钠尿肽>100),]
table(analysis_data$group)

#批量做差异并统计成表格
# Factor(analysis_data$group)
analysis_data$group <- factor(analysis_data$group, levels = c("0","1"))
analysis_data$sex <- factor(analysis_data$sex.y, levels = c("1","0"))
analysis_data <- drop_col(analysis_data,"sex.y")
#数据的分组统计,要求第一列是group，group的值只能是0或者1
analysis_data$group <- factor(analysis_data$group, levels = c("0","1"))
analysis_data$平均血红蛋白浓度 <- analysis_data$平均血红蛋白浓度/100
analysis_data$EFI <- round(analysis_data$RDW_CV/analysis_data$平均血红蛋白浓度,2)
analysis_data$group <- ifelse(analysis_data$B型钠尿肽>=300,1,0)
analysis_data <- analysis_data[analysis_data$B型钠尿肽<100|analysis_data$B型钠尿肽>300,]
analysis_data <- analysis_data[!is.na(analysis_data$EFI),]
analysis_data$EFIQ <- quant(analysis_data$EFI, n = 4,Q = TRUE,round=2)
analysis_data$age60 <- factor(ifelse(analysis_data$age>=60,1,0))

data <- analysis_data[,drop_vector(colnames(analysis_data),c("name","ID","血小板压积","RMR"))]
data <- na_median_sub(data,100)

diff_name_list <- name_list(c("group","sex","EFIQ","age60"))

lapply(diff_name_list, function(x) tbl_summary(data,by = x,
                                               statistic = list(all_continuous()~"{mean} ±{sd}",
                                                                all_categorical() ~ "{n}({p}%)"),
                                               digits = all_categorical() ~ 2)|> 
           add_p()|> 
           as_flex_table() |> 
           save_as_docx(path = paste0(x,"_diff.docx")) )


#logistic回归

#校正1，用年龄和性别
ad_mol1_index <- c("sex","age")

#校正2，用较正1+病史
ad_mol2_index <- c(ad_mol1_index,"X.低密度脂蛋白", "X.高密度脂蛋白", 
                   "X.胆固醇", "X.甘油三酯", "X.葡萄糖", "总胆汁酸", 
                   "间接胆红素", "直接胆红素", "总胆红素", "白球比值", 
                   "球蛋白", "X.白蛋白.溴甲酚绿法.", "X.总蛋白", "a.L.岩藻糖苷酶", 
                   "前白蛋白", "X.Y.谷氨酰转移酶", "碱性磷酸酶", 
                   "AST.ALT", "X.谷草转氨酶", "X.谷丙转氨酶")

#校正3，用校正2+部分临床化验
ad_mol3_index <- c(ad_mol2_index, "中性粒细胞比率", 
                   "大型血小板比率", "平均血小板体积", "血小板分布宽度", 
                   "嗜碱性粒细胞数", "嗜酸性粒细胞数", "中性细胞数", "单核细胞数", 
                   "淋巴细胞数", "嗜酸性粒细胞比率", "中性细胞比率", 
                   "单核细胞比率.1", "淋巴细胞比率.1", "X.血小板", 
                  "红细胞平均体积", "X.红细胞压积", "X.血红蛋白", "X.红细胞", "X.白细胞")

#关注的变量
adj_var <- name_list(c("EFI","红细胞分布宽度CV","平均血红蛋白浓度","EFIQ"))

#批量单因素
adj_glm_model_0 <- adj_log_glm(data=data,
                               signle_list=adj_var,
                               adj_immo=NULL,
                               group_index = "group",
                               outfile="ad_model_0.docx")

#adj模型1
adj_glm_model_1 <- adj_log_glm(data=data,
                               signle_list=adj_var,
                               adj_immo=ad_mol1_index,
                               group_index = "group",
                               outfile="ad_model_1.docx")

#adj模型2
adj_glm_model_2 <- adj_log_glm(data=data,
                               signle_list=adj_var,
                               adj_immo=ad_mol2_index,
                               group_index = "group",
                               outfile="ad_model_2.docx")

#adj模型3
adj_glm_model_3 <- adj_log_glm(data=data,
                               signle_list=adj_var,
                               adj_immo=ad_mol3_index,
                               group_index = "group",
                               outfile="ad_model_3.docx")

#根据性别和年龄，再重新构建一下模型

#age>=60 无校正
ad_mol_list <- list(NULL,
                    drop_vector(ad_mol1_index,"age"),
                    drop_vector(ad_mol2_index,"age"),
                    drop_vector(ad_mol3_index,"age"))

lapply(ad_mol_list, function(x) adj_log_glm_split(data=data,
                                                  signle_list=adj_var,
                                                  adj_immo=x,
                                                  data_split="age60",
                                                  group_index="group"))

#sex 无校正
ad_mol_list <- list(NULL,
                    drop_vector(ad_mol1_index,"sex"),
                    drop_vector(ad_mol2_index,"sex"),
                    drop_vector(ad_mol3_index,"sex"))

lapply(ad_mol_list, function(x) adj_log_glm_split(data=data,
                                                  signle_list=adj_var,
                                                  adj_immo=x,
                                                  data_split="sex",
                                                  group_index="group"))
