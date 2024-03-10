library(readxl)
library(openxlsx)
library(dplyr)

#设置工作目录为R脚本所在的路径
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(digits = 15)

load("cly_lab_result.rdata")

source("E:\\R\\Zl_R_function\\zl_R_Function.r")

#提取数据，病人的门诊号或者住院号
tem_index2 <- read.xlsx("CLY病历门诊号231108-时间节点-v2.xlsx",1)

data <- unite_lab_resu[unite_lab_resu$病历号 %in% tem_index2$病历号,]

cly_1_10_lab <- lapply(list(1,2,3,4,5,6,7,8,9,10), 
                       function(x) 
                           extrade_n_lab(data=data,n=x,
                                                     human_ID = "病历号",
                                                     lab_col="项目名称",
                                                     lab_result="结果"))
#添加表头
i <-1
for(i in 1:length(cly_1_10_lab)){
    tem_data <- cly_1_10_lab[[i]]
    clini_data <- right_join(tem_index2[,c(2:5)],tem_data)
    cly_1_10_lab[[i]] <- clini_data
}

write.xlsx(cly_1_10_lab,file = "cly_1_2_lab.xlsx")
#删除超过800个NA的化验
na_count(first_lab,2)
first_lab1 <- na_delet_data(first_lab,67)
colnames(first_lab1)[1] <- "病历号"
first_lab2 <- left_join(tem_index2,first_lab1)
write.xlsx(first_lab2,file = "cly_lab_first_function_filter.xlsx")
