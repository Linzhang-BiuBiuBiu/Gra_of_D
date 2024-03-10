library(nhanesR)
library(reshape2)
library(dplyr)
library(survey)
library(openxlsx)
library(plyr)
library(dplyr)
library(stringr)
library(rms)
library(pROC)
library(tidyverse)
library(gtsummary)
source("E:\\R\\Zl_R_function\\zl_R_Function_v2.r")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("svrvival_data.Rdata")

#增加RAR和PAR
svrvival_data$MCHC <- svrvival_data$MCHC/10
#svrvival_data <- svrvival_data[!is.na(svrvival_data$MCHC),]
#svrvival_data <- svrvival_data[!is.na(svrvival_data$RDW),]

svrvival_data$HSW <- round((svrvival_data$RDW/svrvival_data$MCHC),2)
#删除0权重数据，删除CHF为NA的数据
svrvival_data <- svrvival_data[!is.na(svrvival_data$CHF),]
svrvival_data <- svrvival_data[(svrvival_data$nhs_wt!=0),]
# Recode(svrvival_data$sex)
svrvival_data$sex <- Recode(svrvival_data$sex,
	"Male::1", 
	"Female::0",
	to.numeric = T)

#识别是否有HF的危险因素，提取变量名称
svrvival_data$pre_HF <- apply(svrvival_data[,c("ASCVD", "angina", "CVD", "heart.attack", "Hypertension")],1,
                                function(x) str_detect(paste(x,collapse = "_"),"1"))

#选择有CHF危险因素的人口
#svrvival_data <- svrvival_data[svrvival_data$pre_HF,]
cla_va <- c("Sex", "Ethnicity", "ASCVD", "angina", "CVD", "heart.attack", "Hypertension")

gre_va <- c("Age", "HSW","MCH","MCHC","ALB","RDW", "WBC", "LymP", "Mon", "SegneP", "EoP", "BaP", 
            "Lym", "MonP", "SeneP", "Eo", "Ba", "RBC", "Hg", "Hem", 
            "Plt", "MPV", "ALT", "AST", "BUN", "Ca", "TC", "HCO3", 
            "GGT", "Glu", "Fe", "TP", "TG", "UA", "Na", "Cl", "GLB", "ntProBNP")

#比较差异
#制作图表
#删除确实超过50000的变量
for(i in colnames(svrvival_data)){
    na_cou_tem <- sum(is.na(svrvival_data[,i]))
    if(na_cou_tem>50000){
        svrvival_data <- drop_col(svrvival_data,i)
    }
}
svrvival_data <- svrvival_data[!is.na(svrvival_data$HSW),]
svrvival_data$HSWQ <- quant(svrvival_data$HSW, n = 4,Q = TRUE,round=2)
svrvival_data$HSWQ.median <- quant.median(svrvival_data$HSW, n = 4,round=2)
svrvival_data$age_60 <- ifelse(svrvival_data$Age>=60,1,0)

nhs <- svy_design(svrvival_data)
differ_variable <- svy_tableone(design = nhs,
                                c_meanPMse   =T,
                                cv=gre_va,
                                gv=c(cla_va,"HSWQ","HSWQ.median"),
                                by = "CHF",
                                xlsx = "nhs_con_HF_difference_IQS.xlsx")

differ_with_RMR_q <- svy_tableone(design = nhs,
                               c_meanPMse=T,
                               cv=gre_va,
                               gv=c(cla_va,"CHF"),
                               by = "HSWQ",
                               xlsx = "nhs_HSWQ_difference.xlsx")

#提取差异变量
differ_variable1 <- differ_variable[!(differ_variable$Pvalue==""),]
diff_varia <- differ_variable1[differ_variable1$Pvalue<0.05,"variable"]
formula_glm <- formula(paste0("CHF~",paste(diff_varia,collapse = " +"))) 

#计算的结果
#将CHF换成数值，否则会报错
svrvival_data$CHF <- Recode(svrvival_data$CHF,
	"0::0", 
	"1::1",
	to.numeric = T)

#批量单因素logistic
nhs <- svy_design(svrvival_data)

single_logi <- svy_uv.logit(design = nhs,y="CHF",
             x=diff_varia,xlsx = "log_con_HF_difference_IQS1.xlsx")

#整理单因素的结果，做多因素logistic
single_logi1 <- single_logi[!(single_logi$`Pr(>|t|)`==""),]
single_logi1 <- single_logi1[1:32,]
single_logi1 <- single_logi1[single_logi1$`Pr(>|t|)`<0.05,"character"]
single_logi1 <- na.omit(single_logi1)

#多个校正模型
#做多因素logistic
#做HSW的多因素logistic
#模型1，用age和Sex校正
single_logi1 <- c("Age", "MCHC", "ALB", "RDW", "LymP", "Mon", "SegneP", 
                  "EoP", "Lym", "MonP", "SeneP", "Eo", "Ba", "RBC", "Hg", "Hem", 
                  "Plt", "MPV", "BUN", "Ca", "TC", "HCO3", "GGT", "Glu", "Fe", 
                  "TP", "TG", "UA", "Cl", "GLB")

variab <- c("MCHC","RDW","HSW","HSWQ")

adj_log_result <- c()
for(target in variab){
    m1 <- c("Sex","Age")
    m2 <- unique(c(m1,"heart.attack","ASCVD","angina","Hypertension","Ethnicity"))
    m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
    
    log_adi_m1 <- formula(paste0("CHF~",paste(c(target,m1),collapse = "+")))
    log_adi_m2 <- formula(paste0("CHF~",paste(c(target,m2),collapse = "+")))
    log_adi_m3 <- formula(paste0("CHF~",paste(c(target,m3),collapse = "+")))
    log_adj_list <- list(log_adi_m1,log_adi_m2,log_adi_m3)
    
    adj_log_tem_result <- lapply(log_adj_list, 
                                 function(x) 
                                     svyglm(x, design=nhs,
                                            family =quasibinomial() ) 
                                 |> reg_table())
    names_index <- paste(paste0(target),c(1:3))
    
    names(adj_log_tem_result) <- names_index
    
    adj_log_result <- c(adj_log_tem_result,adj_log_result)}
write.xlsx(adj_log_result,file = "adj_log_result.xlsx")
write.xlsx(single_logi,file = "uni_log_result.xlsx")

#单独将age(60)，sex，DM分层，观察HSW与CVD患病率关系
single_logi1 <- c("Age", "HSW","MCHC", "ALB", "RDW", "LymP", "Mon", "SegneP", 
                  "EoP", "Lym", "MonP", "SeneP", "Eo", "Ba", "RBC", "Hg", "Hem", 
                  "Plt", "MPV", "BUN", "Ca", "TC", "HCO3", "GGT", "Glu", "Fe", 
                  "TP", "TG", "UA", "Cl", "GLB","HSWQ")

svrvival_data$age_60 <- ifelse(svrvival_data$Age>60,">60","<=60")
sprci_va <- c(">60","<=60","Female","Male")

adj_log_spe_result <- c()
uni_log_spe_result <- list()
for(i in 1:length(sprci_va)){
    
    index <- ifelse(i<=2,"age_60","Sex")
    
    ana_data <- svrvival_data[which(svrvival_data[,index]==sprci_va[i]),]
    ana_data <- drop_col(ana_data,index)
    
    nhs <- svy_design(ana_data)
    
    single_spe_logi <- svy_uv.logit(design = nhs,y="CVD",
                                    x=c("HSW","HSWQ"))
    
    single_spe_log_name <- paste(index,"_",sprci_va[i])
    uni_log_spe_result[[single_spe_log_name]] <- single_spe_logi
    
    variab <- c("RDW","MCHC",index)
    
    variab1 <- c("HSW","HSWQ")
    
    for(target in variab1){
        if (index=="age_60") {
            m1 <- c("Sex",target)  
        }else{
            m1 <- drop_vector(c("Age","Sex",target),index)
        }
        
        m2 <- drop_vector(unique(c(m1,"heart.attack","ASCVD","angina","Hypertension","Ethnicity")),index)
        m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
        
        log_adi_m1 <- formula(paste0("CHF~",paste(c(target,m1),collapse = "+")))
        log_adi_m2 <- formula(paste0("CHF~",paste(c(target,m2),collapse = "+")))
        log_adi_m3 <- formula(paste0("CHF~",paste(c(target,m3),collapse = "+")))
        log_adj_list <- list(log_adi_m1,log_adi_m2,log_adi_m3)
        
        adj_log_tem_result <- lapply(log_adj_list, 
                                     function(x) 
                                         svyglm(x, design=nhs,
                                                family =quasibinomial() ) 
                                     |> reg_table())
        names_index <- paste(paste0(index,"_",target,"_",sprci_va[i]),c(1:3))
        
        names(adj_log_tem_result) <- names_index
        
        adj_log_spe_result <- c(adj_log_tem_result,adj_log_spe_result)}
    
}
write.xlsx(uni_log_spe_result,file="uni_log_spe_result.xlsx")
write.xlsx(adj_log_spe_result,file = "adj_log_result_2.xlsx")

#整理单因素结果，准备做cox
svrvival_data1 <- svrvival_data[!svrvival_data$mortstat %in% NA,]
svrvival_data1 <- svrvival_data1[svrvival_data1$CHF==1,]
svrvival_data1 <- svrvival_data1[svrvival_data1$eligstat=="Eligible",]

# Recode(svrvival_data1$mortstat)
svrvival_data1$mortstat <- Recode(svrvival_data1$mortstat,
	"Assumed deceased::1", 
	"Assumed alive::0",
	to.numeric = T)

#多因素cox回归
var_count <- apply(svrvival_data1, 2, 
                   function(x) 
                       ifelse(length(unique(x))>5,"conti","cla"))

gre_va <- names(var_count[var_count=="conti"])
gre_va <- drop_vector(gre_va,c("Year","ucod_leading","sdmvstra","seqn",
                               "wtmec2yr","wtmec4yr","ntProBNP",
                               "nhs_wt","permth_exm","permth_int","MCH","Hg"))

cla_va <- names(var_count[var_count=="cla"])
cla_va <- drop_vector(cla_va,c("CHF","diabetes","CVD","sdmvpsu","hyperten","mortstat",
                               "eligstat","pre_HF","HSWQ.median","CKD_prognosis","HSWQ"))

#做表
hs <- svy_design(svrvival_data1)
differ_variable <- svy_tableone(design = nhs,
                                c_meanPMse   =T,
                                cv=c(gre_va,"permth_int"),
                                gv=c(cla_va,"HSWQ","HSWQ.median"),
                                by = "mortstat")

differ_with_RMR_q <- svy_tableone(design = nhs,
                                  c_meanPMse=T,
                                  cv=c(gre_va,"permth_int"),
                                  gv=c(cla_va,"mortstat"),
                                  by="HSWQ")

write.xlsx(list(differ_variable,differ_with_RMR_q),
           file = "NHANES_diff_follow.xlsx")

#单因素cox
uni_cox <- svy_uv.cox(nhs,time = "permth_int",
                      status = "mortstat",
                      x=c(gre_va,cla_va,"HSWQ"))

#多因素cox回归
#多因素cox回归
variab <- c("HSW","RDW","MCHC","HSWQ")
adj_cox_result <- c()

for(target in variab){
    m1 <- c("Age","Sex")
    m2 <- unique(c(m1,"heart.attack","ASCVD","angina","Hypertension","Ethnicity"))
    m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
    
    cox_adi_m1 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m1),collapse = "+")))
    cox_adi_m2 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m2),collapse = "+")))
    cox_adi_m3 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m3),collapse = "+")))
    
    cox_adi_list <- list(cox_adi_m1,cox_adi_m2,cox_adi_m3)
    
    cox_log_tem_result <- lapply(cox_adi_list, 
                                 function(x) 
                                     svycoxph(x, design=nhs) 
                                 |> reg_table())
    
    adj_cox_result <- c(cox_log_tem_result,adj_cox_result)}
write.xlsx(adj_cox_result,file = "adj_cox_result.xlsx")
write.xlsx(uni_cox,file = "uni_cox_result.xlsx")

#单独将age(60)，sex，DM分层，观察HSW与CVD患病率关系

sprci_va <- c(">60","<=60","Female","Male")

adj_cox_spe_result <- c()
uni_cox_spe_result <- c()
for(i in 1:length(sprci_va)){
    index <- ifelse(i<=2,"age_60","Sex")
    
    ana_data <- svrvival_data1[which(svrvival_data1[,index]==sprci_va[i]),]
    ana_data <- drop_col(ana_data,index)
    
    nhs <- svy_design(ana_data)
    
    uni_spe_cox <- svy_uv.cox(nhs,time = "permth_int",
                              status = "mortstat",
                              x=c("HSW","HSWQ"))
    
    single_spe_cox_name <- paste(index,"_",sprci_va[i])
    
    uni_cox_spe_result[[single_spe_cox_name]] <- uni_spe_cox
    
    variab <- c("RDW","MCHC",index)
    
    variab1 <- c("HSW","HSWQ")
    
    for(target in variab1){
        if (index=="age_60") {
            m1 <- c("Sex",target)  
        }else{
            m1 <- drop_vector(c("Age","Sex",target),index)
        }
        
        m2 <- drop_vector(unique(c(m1,"heart.attack","ASCVD","angina","Hypertension","Ethnicity")),index)
        m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
        
        cox_adi_m1 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m1),collapse = "+")))
        cox_adi_m2 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m2),collapse = "+")))
        cox_adi_m3 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m3),collapse = "+")))
        
        cox_adi_list <- list(cox_adi_m1,cox_adi_m2,cox_adi_m3)
        
        cox_tem_result <- lapply(cox_adi_list, 
                                 function(x) 
                                     svycoxph(x, design=nhs) 
                                 |> reg_table())
        
        names_index <- paste(paste0(index,"_",target,"_",sprci_va[i]),c(1:3))
        
        names(cox_tem_result) <- names_index
        
        adj_cox_spe_result <- c(adj_cox_spe_result,cox_tem_result)}
    
}
write.xlsx(adj_cox_spe_result,file = "adj_cox_result_3.xlsx")
write.xlsx(uni_cox_spe_result,file = "uni_cox_spe_result.xlsx")
#生存分析
CHF_survi <-subset(svrvival_data1,CHF==1) 
CHF_survi <- CHF_survi[CHF_survi$nhs_wt!=0,]
CHF_survi <- CHF_survi[!is.na(CHF_survi$nhs_wt),]

#RCS寻找最佳的界值
ddist <- datadist(CHF_survi)
options(datadist='ddist')

serch_index <- "RDW"
adju_index <-  "Age"

S <- Surv(CHF_survi$permth_int,CHF_survi$mortstat==1)
for (knot in 3:10) {
    cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",knot,")+",adju_index))
    
    fit <- cph(cph_formula,data=CHF_survi,x= TRUE, y= TRUE, surv = TRUE)
    tmp <- extractAIC(fit)
    if(knot==3){AIC=tmp[2];nk=3}
    if(tmp[2]<AIC){AIC=tmp[2];nk=knot}
}

#---------基于cox绘制 单因素LnAl 的 RCS
pacman::p_load(rms,survminer,ggplot2,ggsci)
#构建cph 函数获取rcs, 单指标LnAl；也可校正其他因素 S ~ rcs(LnAl,4)+ BMI 等
cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",5,")+",adju_index))

fit.RAR <- cph(cph_formula, x=TRUE, y=TRUE,data=CHF_survi)
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
Pre_HR.RAR <-rms::Predict(fit.RAR,RDW,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
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
              colour="purple")+
    scale_color_nejm()+ ##采用ggsci包中英格兰调色，也可以其他
    geom_ribbon(data=Pre_HR.RAR,
                aes(MCHC, ymin=lower,ymax=upper,fill="purple"),alpha=0.1)+
    theme_classic()+
    scale_fill_nejm()+
    geom_hline(yintercept=1,linetype=2,size=0.75) #y=1水平线
labs(title ="风险随LnAl变化曲RCS",
     x="RAR", 
     y="HR (95%CI)")

pdf("rcs_RDW.pdf",width = 6,height = 4)
rcs_all
dev.off()

CHF_survi$RMR_4.09 <- ifelse(CHF_survi$RMR>4.09,"High","Low")
CHF_survi$Age_70 <- ifelse(CHF_survi$Age>70,"High","Low")
CHF_survi$RDW_13.70 <- ifelse(CHF_survi$RDW>13.70,"High","Low")
CHF_survi$MCHC_3.36 <- ifelse(CHF_survi$MCHC > 3.36,"High","Low")

#根据分层统计数据
nhs <- svy_design(CHF_survi)
svy_population(nhs)

pdf(file = "survival1.pdf",width = 8,height = 6)
RMR <- svykm(Surv(permth_int,mortstat)~RMR_4.09,design = nhs)
svy_kmplot(RMR,ci = T)

RDW <- svykm(Surv(permth_int,mortstat)~RDW_13.70,design = nhs)
svy_kmplot(RDW,ci = T)

Age <- svykm(Surv(permth_int,mortstat)~Age_70,design = nhs)
svy_kmplot(Age,ci = T)

MCHC <- svykm(Surv(permth_int,mortstat)~MCHC_3.358,design = nhs)
svy_kmplot(MCHC,ci = T)
dev.off()
#分段生存率
#将RAR多次分段，比较病死率
CHF_patients <- svrvival_data1[svrvival_data1$CHF==1,]
CHF_patients <- CHF_patients[!is.na(CHF_patients$RMR),]
CHF_patients$RMRQ <- quant(CHF_patients$RMR, n = 10,Q = TRUE,round=3)
CHF_patients$RMRQ.median <- quant.median(CHF_patients$RMR, n = 10,round=3)

deth_rar_q <- c()
for (i in unique(CHF_patients$RMRQ)) {
    mordi_data <-CHF_patients[CHF_patients$RMRQ==i,]
    dedth_indi <- mordi_data[mordi_data$mortstat==1,]
    deth_perc <- nrow(dedth_indi)/nrow(mordi_data)
    q_deth <- c(Q=i,mor=deth_perc)
    deth_rar_q <- rbind(q_deth,deth_rar_q)
}

plot(deth_rar_q[,1],deth_rar_q[,2])

write.xlsx(list(deth_rar_q,CHF_patients),file="NHANES_CH_PATIENTS.xlsx")
RMR_15 <- quantile_mortality(CHF_patients,"mortstat","RMR",15)
RDW_15 <- quantile_mortality(CHF_patients,"mortstat","RDW",15)
MCHC_15 <- quantile_mortality(CHF_patients,"mortstat","MCHC",15)

all_mor <- merge(RMR_15,RDW_15,by="Q")
all_mor <- merge(all_mor,MCHC_15,by="Q")
colnames(all_mor) <- c("Q","RMR_15","RDW_15","MCHC_15")

write.xlsx(all_mor,file = "all_mor.xlsx")

ck <- c( "CHF", "ASCVD", "angina", "CVD", "heart.attack", "Hypertension", 
         "Age", "Sex", "Ethnicity", "WBC", "LymP", "MonP", "SegneP", "EoP", "BaP", "Lym", "Mon", 
         "SeneP", "Eo", "Ba", "RBC", "Hg", "Hem", "MCV", "MCH", "MCHC", 
         "RDW", "Plt", "MPV", "ALB", "ALT", "AST", "BUN", "Ca", "TC", 
         "HCO3", "GGT", "Glu", "Fe", "TP", "TG", "UA", "Na", "Cl", "GLB", 
         "ntProBNP","RMR", "RMR_4.09","Age_70","RDW_13.70","MCHC_3.358","straty")
colnames(CHF_survi)
d <- CHF_survi |> select(ck) |> tbl_summary(by=RMR_4.09,
                                            missing = "no",
                                            statistic = list(all_continuous() ~ "{mean}±{sd}",
                                                             all_categorical() ~ "{n}/{N}({p}%)"),
                                            digits = all_categorical() ~ 2,
                                            missing_text = "(Missing)",
                                            ) |> add_p(test = list(all_continuous() ~ "t.test",            
                                                                   all_categorical() ~ "fisher.test")) %>% 
    separate_p_footnotes() %>%            
    modify_caption("**Table 1. Patient Characteristics**") %>%            
    bold_labels() 

write.xlsx(d,file = "CHF_level_diff.xlsx")

