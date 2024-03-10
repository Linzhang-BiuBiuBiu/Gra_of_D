library(xlsx)
library(readxl)
library(openxlsx)
library(stringr)
library(dplyr)
library(nhanesR)
library(DataEditR)
library(reshape2)
library(data.table)
library(mice)
library(survival) 
library(survminer)
library(survey)
library(rms)
library(haven)
library(pacman)
library(pROC)
library(gtsummary)          
library(tidyverse)
library(survival)
library(flextable)

source("E:\\R\\Zl_R_function\\zl_R_Function_v2.r")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rx <- read.csv("Mimic_iii_clean_data.csv")

icd_merge <- read.csv("icd_merge.csv")
rx1 <- left_join(rx,icd_merge)

group_patients <- read.csv("all_patiens.csv")

primary_data  <- rx1[which(rx1$HADM_ID %in% group_patients$HADM_ID),]
primary_data <- merge(group_patients,primary_data[,-3])

analysis_data <- primary_data[,c(-1,-3)]
analysis_data <- analysis_data[analysis_data$age>18,]

# Factor(analysis_data$Group)
analysis_data$Group <- factor(analysis_data$Group, levels = c("0","1"))
analysis_data$Group <- as.factor(analysis_data$Group)
# Factor(analysis_data$GENDER)
analysis_data$GENDER <- ifelse(analysis_data$GENDER=="F",0,1)
# Factor(analysis_data$GENDER)
analysis_data$GENDER <- factor(analysis_data$GENDER, levels = c(0,1))

#删除既没有RDW，也没有MCH的数据
analysis_data1 <- analysis_data[!is.na(analysis_data$RDW_Bl_H),]
analysis_data1 <- analysis_data1[!is.na(analysis_data1$MCH_Bl_H.1),]
analysis_data1 <- analysis_data1[!(analysis_data1$MCH_Bl_H.1 %in% 0),]
analysis_data1 <- analysis_data1[!(analysis_data1$RDW_Bl_H %in% 0),]

#删除缺失数据超过0.05的化验,删除糖尿病的数据
for(i in colnames(analysis_data1)[-1:-9]){
    na_count1 <- sum(is.na(analysis_data1[,i]))
    if (na_count1 >= nrow(analysis_data1)*0.05) {
        analysis_data1 <- analysis_data1[,!(colnames(analysis_data1) %in% i)]
    }else{
        na_where <- is.na(analysis_data1[,i])
        mean_na <- mean(as.numeric(na.omit(analysis_data1[,i])))
        analysis_data1[na_where,i] <- mean_na
    }
    
}

#读取Item的简称
Iten_id <- read.csv("item_more.csv")
#将对应的化验替换成简称

analysis_data2 <- data.frame(colnames(analysis_data1),1:36)
colnames(analysis_data2) <- c("simit","kk")
analysis_data3 <- merge(analysis_data2,Iten_id[,c("simit","simply")],all = T)
analysis_data3 <- analysis_data3[order(analysis_data3$kk),]

analysis_data3$simply <- ifelse(is.na(analysis_data3$simply),
                                analysis_data3$simit,
                                analysis_data3$simply)

colnames(analysis_data1) <- analysis_data3$simply

#批量做差异并统计成表格
#数据的分组统计,要求第一列是group，group的值只能是0或者1
analysis_data1 <- drop_col(analysis_data1,"RAR")
analysis_data1 <- apply(analysis_data1, 2, as.numeric)
analysis_data1 <- as.data.frame(analysis_data1)
analysis_data1$MCHC <- analysis_data1$MCHC/10
analysis_data1$RMR <- round(analysis_data1$RDW/analysis_data1$MCHC,2)

analysis_data1$RMRQ <- quant(analysis_data1$RMR, n = 4,Q = TRUE,round=2)
analysis_data1$RMRQ.median <- quant.median(analysis_data1$RMR, n = 4,round=2)
analysis_data1$age60 <- factor(ifelse(analysis_data1$age >=60,1,0))


diff_name_list <- name_list(c("Group","GENDER","RMRQ","age60"))

#做差异
lapply(diff_name_list, function(x) tbl_summary(analysis_data1,by = x,
                                               statistic = list(all_continuous()~"{mean} ±{sd}",
                                                                all_categorical() ~ "{n}({p}%)"),
                                               digits = all_categorical() ~ 2)|> 
           add_p()|> 
           as_flex_table() |> 
           save_as_docx(path = paste0(x,"_diff.docx")) )


#logistic回归
#校正1，用年龄和性别
ad_mol1_index <- c("GENDER","age")

#校正2，用较正1+病史
ad_mol2_index <- c(ad_mol1_index, "CAD", "CKD", "CRD", "DIABETES", 
                   "HYPERTENSION")

#校正3，用校正2+部分临床化验
ad_mol3_index <- c(ad_mol2_index, "ALB", "Pho", "Mg", "Ca", "WBC", "PTT", "PT", 
                   "Pla", "INR", "HCT", "BUN", "Na", "K", "Glu", "SCR", "Cl", "HCO3", "AG")

#关注的变量
adj_var <- name_list(c("RMR","RDW","MCHC","RMRQ"))

adj_name_list <- list(NULL,ad_mol1_index,
                      ad_mol2_index,ad_mol3_index)

lapply(adj_name_list,function(x) adj_log_glm(data=analysis_data1,
                                            signle_list=adj_var,
                                            adj_immo=x,
                                            group_index = "Group") )

#敏感性分析
#age<=60 无校正
ad_mol_list <- list(NULL,
                    drop_vector(ad_mol1_index,"age"),
                    drop_vector(ad_mol2_index,"age"),
                    drop_vector(ad_mol3_index,"age"))

lapply(ad_mol_list, function(x) adj_log_glm_split(data=analysis_data1,
                                                  signle_list=adj_var,
                                                  adj_immo=x,
                                                  data_split="age60",
                                                  group_index="Group"))
#sex 
ad_mol_list <- list(NULL,
                    drop_vector(ad_mol1_index,"GENDER"),
                    drop_vector(ad_mol2_index,"GENDER"),
                    drop_vector(ad_mol3_index,"GENDER"))

lapply(ad_mol_list, function(x) adj_log_glm_split(data=analysis_data1,
                                                  signle_list=adj_var,
                                                  adj_immo=x,
                                                  data_split="GENDER",
                                                  group_index="Group"))

#筛选随访数据
HF_initial_data <- analysis_data1[analysis_data1$Group==1,][,c(-1,-8)]#将NtProBNP删除
HF_initial_data <- HF_initial_data[HF_initial_data$EXPIRE_FLAG==1,]
HF_initial_data$EXPIRE_FLAG <- ifelse(HF_initial_data$liv_time<1087,1,0)

#分析HF患者的两组特征

diff_name_list <- name_list(c("EXPIRE_FLAG","GENDER","RMRQ","age60"))
lapply(diff_name_list, function(x) tbl_summary(HF_initial_data,by = x,
                                               statistic = list(all_continuous()~"{mean} ±{sd}",
                                                                all_categorical() ~ "{n}({p}%)"),
                                               digits = all_categorical() ~ 2)|> 
           add_p()|> 
           as_flex_table() |> 
           save_as_docx(path = paste0(x,"_follow_diff.docx")) )


#Cox回归
lapply(adj_name_list,function(x)adj_cox_reg (data=HF_initial_data, signle_list=adj_var,
                                            statue="EXPIRE_FLAG",adj_immo=x,live_time="liv_time"))




H_chi_t_data <- differenc_test(HF_initial_data,"EXPIRE_FLAG")

tbl_summary(HF_initial_data,by = "EXPIRE_FLAG",
            statistic = all_continuous()~"{mean}±{sd}") |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "HSW_group_follow.docx")

tbl_summary(HF_initial_data,by = "RMRQ",
            statistic = all_continuous()~"{mean}±{sd}") |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "HSW_group_RMRQ_follow.docx")

#单因素Cox
uni_cox_data <- uni_cox_test(HF_initial_data,"EXPIRE_FLAG","liv_time")

#多因素cox
d_v <- c("RDW","MCHC","NTproBNP","HYPERTENSION","RMR","RMRQ.dedian")
adj_cox_data_RMR <- adj_cox_test(HF_initial_data,"EXPIRE_FLAG","liv_time",d_v)

d_v <- c("RMR","NTproBNP","HYPERTENSION")
adj_cox_data_RDW <- adj_cox_test(HF_initial_data,"EXPIRE_FLAG","liv_time",d_v)

adj_log_RMR <- cbind(rownames(adj_log_RMR),adj_log_RMR) 
adj_log_rdw_mch <- cbind(rownames(adj_log_rdw_mch),adj_log_rdw_mch)
mimic_log_cox_result <- list(H_chi_t_data=H_chi_t_data,
                             uni_log=uni_glm_result,adj_log1=adj_log_RMR,
                             adj_log2=adj_log_rdw_mch,uni_cox=uni_cox_data,
                             adj_cox1=adj_cox_data_RMR,adj_cox2=adj_cox_data_RDW)

write.xlsx(mimic_log_cox_result,file = "mimic_log_cox_result1.xlsx")
write.xlsx(uni_cox_data,file = "TEM_UNI_cox_result1.xlsx")

#寻找RAR的最佳RCS
ddist <- datadist(HF_initial_data)
options(datadist="ddist")
serch_index <- "MCHC"
adju_index <-  "INR"

S <- Surv(HF_initial_data$liv_time,HF_initial_data$EXPIRE_FLAG==1)
for (knot in 3:10) {
    cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",knot,")+",adju_index))
    
    fit <- cph(cph_formula,data=HF_initial_data,x= TRUE, y= TRUE, surv = TRUE)
    tmp <- extractAIC(fit)
    if(knot==3){AIC=tmp[2];nk=3}
    if(tmp[2]<AIC){AIC=tmp[2];nk=knot}
}


#---------基于cox绘制 单因素LnAl 的 RCS
pacman::p_load(rms,survminer,ggplot2,ggsci)
#构建cph 函数获取rcs, 单指标LnAl；也可校正其他因素 S ~ rcs(LnAl,4)+ BMI 等
cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",nk,")+",adju_index))

fit.RAR <- cph(cph_formula, x=TRUE, y=TRUE,data=HF_initial_data)
# PH 检验,P>0.05 符合PH假设
cox.zph(fit.RAR, "rank")    
#残差图的横轴是时间，纵轴是残差，残差均匀分布则表示残差与时间相互独立，满足ph假设
ggcoxzph(cox.zph(fit.RAR, "rank"))

pdf_name <- paste0("RCS_",serch_index,"-",adju_index,".pdf")
pdf(pdf_name,width = 5,height = 5)
ggcoxzph(cox.zph(fit.RAR, "rank")) 
dev.off()

#非线性检验p-non-linear，0.9289
anova(fit.RAR)                      
p <-round(anova(fit.RAR)[,3],3)
## HR计算，fun是转化函数
Pre_HR.RAR <-rms::Predict(fit.RAR,MCHC,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
ggplot(Pre_HR.RAR)
# y-hat=1.00000, 对应LnAl= 5.006
#View(Pre_HR.RAR)

# 全人群的COX模型 RCS绘图
rcs_all <-  ggplot()+
    geom_line(data=Pre_HR.RAR,
              aes(MCHC,yhat),# x y轴距
              linetype="solid",#曲线加粗
              size=1,
              alpha=0.7,
              colour="green")+
    scale_fill_nejm()+ ##采用ggsci包中英格兰调色，也可以其他
    geom_ribbon(data=Pre_HR.RAR,
               aes(MCHC, ymin=lower,ymax=upper,fill="zx"),alpha=0.1)+
    theme_classic()+
    scale_fill_nejm()+
    geom_hline(yintercept=1,linetype=2,size=0.75) #y=1水平线+
labs(title ="风险随LnAl变化曲RCS",
     x="RAR", 
     y="HR (95%CI)"
)
pdf_name1 <- paste0("rcs_",serch_index,"_",adju_index,"_.pdf")
pdf(pdf_name1,width = 5,height = 4)
rcs_all
dev.off()


 
#KM曲线
HF_initial_data$RMR_4.7 <- ifelse(HF_initial_data$RMR>4.7,"high","low")
HF_initial_data$ALB_3.2 <- ifelse(HF_initial_data$ALB>3.2,"high","low")
HF_initial_data$RDW_15.6 <- ifelse(HF_initial_data$RDW>15.6,"high","low")
HF_initial_data$Age_68 <- ifelse(HF_initial_data$age>68,"high","low")
HF_initial_data$MCHC_3.3 <- ifelse(HF_initial_data$MCHC>3.3,"high","low")

#KM分析的批量做图

Var_KM <- grep("_KM",colnames(HF_initial_data),value = T)
Var_KM <- c("RMR_4.7","ALB_3.2","RDW_15.6","Age_68","MCHC_3.3")
km_formula <- sapply(Var_KM,function(x) as.formula(paste0('Surv(liv_time, EXPIRE_FLAG)~',x)))
km_analys <- lapply(km_formula,function(x) survfit(x,data=HF_initial_data))

for(i in Var_KM){
    
    fit_formula <- as.formula(paste0('Surv(liv_time, EXPIRE_FLAG)~',i))
    fit <- survfit(fit_formula,data=HF_initial_data)
    fit$call$formula <- fit_formula

    fig <- ggsurvplot(fit,
               pval = TRUE, conf.int = TRUE,
               risk.table = TRUE, 
               risk.table.col = "strata", 
               linetype = "strata", 
               surv.median.line = "hv", # 同时显示垂直和水平参考线
               ggtheme = theme_bw(), 
               palette = c("#E7B800", "#2E9FDF"))
    
    pdf(file = paste0(i,"_mimic.pdf"),width=8, height=6)
    print(fig,newpage = FALSE)
    dev.off()
}

#构建列线图
#将数据打包好
HF_initial_data$RARQ <- quant(HF_initial_data$RAR, n = 4,Q = TRUE,round=5)
HF_initial_data$RARQ.median <- quant.median(HF_initial_data$RAR, n = 4,round=2)

ddist <- datadist(HF_initial_data)
options(datadist='ddist')

#构建多因素的Cox回归模型
cox <- cph(Surv(liv_time, EXPIRE_FLAG) ~ age+RARQ.median+ALB+RDW,
    data = HF_initial_data,x=T,y=T,surv = T)

surv <- Survival(cox)
sur_1_year<-function(x)surv(365,lp=x)#1年生存
sur_2_year<-function(x)surv(730,lp=x)#2年生存
sur_3_year<-function(x)surv(1095,lp=x)#3年生存
#sur_4_year<-function(x)surv(1460,lp=x)#4年生存
#sur_5_year<-function(x)surv(1825,lp=x)#5年生存

# 做列线图
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_2_year,sur_3_year),
                    lp= F,
                    funlabel=c('1-Year Survival','2-Year survival','3-Year survival'),maxscale=100,fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))

pdf("nomogram.pdf",width=8, height=6)
plot(nom_sur)
dev.off()

#将RAR多次分段，比较病死率
HF_initial_data$RARQ <- quant(HF_initial_data$RMR, n = 20,Q = TRUE,round=2)
HF_initial_data$RARQ.median <- quant.median(HF_initial_data$RMR, n = 20,round=2)

deth_rar_P <- c()
for (i in unique(HF_initial_data$RARQ.median)) {
    mordi_data <-HF_initial_data[HF_initial_data$RARQ.median==i,]
    dedth_indi <- mordi_data[mordi_data$EXPIRE_FLAG==1,]
    deth_perc <- nrow(dedth_indi)/nrow(mordi_data)
    q_deth <- c(Q=i,number=nrow(mordi_data),dead_number=nrow(dedth_indi),mor=deth_perc)
    deth_rar_P <- rbind(q_deth,deth_rar_P)
}




deadt_all <- list(RMR_death_rate,RDW_death_rate,MCHC_death_rate)

write.xlsx(deadt_all,file = "deadt_RMR_RDW_MCHC.xlsx")
