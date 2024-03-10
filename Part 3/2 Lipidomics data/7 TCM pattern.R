setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(gtsummary)          
library(flextable)
library(openxlsx)
library(mice)
source("E:\\R\\Zl_R_function\\zl_R_Function.r")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#读取脂质数据以及临床数据
lipid_coef <- read.xlsx("clini_coef_data_RBC.xlsx")

clini_data <- read.xlsx("clin_data_v3.xlsx",3)

clini_data$NUM <- as.character(clini_data$NUM)
#数据合并，删除临床数据中的分组
lipid_clin_coef <- left_join(lipid_coef,clini_data[,-11:-13])
lipid_clin_coef$HSW <- (as.numeric(lipid_clin_coef$红细胞分布宽度CV)/as.numeric(lipid_clin_coef$平均血红蛋白浓度))*100

#选择需要的列
colnames(lipid_clin_coef) |> dput()
clin_index <- c("分组", "中医证型", "性别",  "冠心病", 
  "房颤", "肾功能不全", "高血压", "糖尿病", "高脂血症","年龄", "DCER(24:1)", 
  "DCER(22:0)", "HexCer.d18:1/12:0", "DCER(24:0)", "PE(P-18:0/20:4)", 
  "HCER(24:1)", "PE(O-18:0/20:4)", "DCER(16:0)", "PE(16:0/18:2)", 
  "PE(O-18:0/22:5)", "PE(P-18:0/22:5)", "PE(P-18:0/22:4)", "PE(O-18:0/22:4)", 
  "PE(P-18:2/18:2)", "PE(18:2/18:2)", "HCER(16:0)", "PE(16:0/20:2)", 
  "PE(P-18:0/22:6)",  "coef", "HSW",
   "BNP", "EF", "左房", "左室舒末内径", 
  "右房", "右室", "A/E", "肺动脉压mmHg", "左室缩末内径", 
  "室间隔厚度", "左室后壁厚度", "肺动脉:内径", 
  "肺动脉流速m/s", "白细胞", "上皮细胞", 
  "肌钙蛋白T", "葡萄糖", "★白细胞", "★红细胞", 
  "★血红蛋白", "★红细胞压积", "红细胞平均体积", 
  "平均血红蛋白量", "平均血红蛋白浓度", "★血小板", 
  "淋巴细胞比率", "单核细胞比率", "中性细胞比率", 
  "嗜酸性粒细胞比率", "嗜碱性粒细胞比率", "淋巴细胞数", 
  "单核细胞数", "中性细胞数", "嗜酸性粒细胞数", 
  "嗜碱性粒细胞数", "红细胞分布宽度CV", "红细胞分布宽度SD", 
  "血小板分布宽度", "平均血小板体积", "大型血小板比率", 
  "血小板压积", "中性粒细胞比率", "★谷草转氨酶", 
  "肌酸激酶同工酶", "★超敏促甲状腺素", "★凝血酶原时间", 
  "★凝血酶原国际标准化比值", "活化部分凝血活酶时间", 
  "C反应蛋白", "★游离三碘甲状腺原氨酸", "★游离甲状腺素", 
  "钾", "钠", "氯", "★葡萄糖", "★谷丙转氨酶", "AST/ALT", 
  "★白蛋白溴甲酚绿法", "★总蛋白", "球蛋白", "白球比值", 
  "总胆红素", "直接胆红素", "间接胆红素", "总胆汁酸", 
  "★甘油三酯", "★胆固醇", "★高密度脂蛋白", "★低密度脂蛋白", 
  "碱性磷酸酶", "★Y谷氨酰转移酶", "前白蛋白", 
  "aL岩藻糖苷酶", "同型半胱氨酸", "比重", "尿蛋白", 
  "尿胆原", "酮体", "PH", "VC", "微白蛋白", "隐血", 
  "胆红素", "红细胞", "结晶", "管型", "白细胞高倍视野", 
  "红细胞高倍视野", "上皮细胞高倍视野", "管型低倍视野", 
  "小圆细胞", "类酵母菌", "未溶红细胞绝对值", "未溶红细胞比率", 
  "病理管型", "电导率", "总粒子数", "★乳酸脱氢酶", 
  "★肌酸激酶", "a羟丁酸脱氢酶", "铁", "镁", "★钙", 
  "★磷", "锌", "潜血免疫法", "粪便颜色", "粪便性状", 
  "★糖化血红蛋白", "肌酐", "尿素", "尿酸", "纤维蛋白原", 
  "★钾", "★钠", "★氯", "嗜酸性细胞比率", "D一二聚体", 
  "纤维蛋白原降解产物FDP", "★肌酐酶法", "★尿素", 
  "★尿酸", "二氧化碳结合力", "ABO血", "RH血型", "★乙型肝炎表面抗原", 
  "★乙型肝炎表面抗体", "乙型肝炎e抗原", "乙型肝炎e抗体", 
  "乙型肝炎核心抗体", "乙肝病毒Pre.S1抗原", "★丙型肝炎病毒抗体", 
  "谷草转氨酶", "肌酸激酶", "乳酸脱氢酶")


lipid_clin_coef <- lipid_clin_coef[,clin_index]

#删除NA值太多的列，超过30%的数据
lipid_clin_coef1 <- na_delet_data(lipid_clin_coef,nrow(lipid_clin_coef)*0.3)

#化验的特殊符号处理
ckkk <- gsub("★","P",colnames(lipid_clin_coef1))
ckkk <- gsub("\\(","",ckkk)
ckkk <- gsub("\\)","",ckkk)
ckkk <- gsub("\\/","_",ckkk)
ckkk <- gsub("\\.","",ckkk)
ckkk <- gsub("\\:","_",ckkk)
ckkk <- gsub("-","_",ckkk)
colnames(lipid_clin_coef1) <- ckkk

#写出数据,转换成numeric
lipid_clin_coef1[,c(-1:-9)] <- apply(lipid_clin_coef1[,c(-1:-9)],2,as.numeric)

lipid_clin_coef1[,c(1:9)] <- apply(lipid_clin_coef1[,c(1:9)],2,as.factor)

#mice多重插补
mice_lipid_clin_coef1 <- mice(lipid_clin_coef1,m=1,seed = 5)

#读取多重插补信息,并替换列名
lipid_clin_coef2 <- complete(mice_lipid_clin_coef1,action = 1)
lipid_clin_coef2[,"HSW"] <- (lipid_clin_coef2$红细胞分布宽度CV/lipid_clin_coef2$平均血红蛋白浓度)*100
lipid_clin_coef2$coef <- lipid_clin_coef2$coef^2

#先筛选出CD期的患者
tcm_cd <- lipid_clin_coef2[lipid_clin_coef2$分组 %in% c("C","D"),]

#统计数据
tbl_summary(tcm_cd[,-1],by = "中医证型",
            statistic = list(all_continuous()~"{median} ({p25}, {p75})",
                             all_categorical() ~ "{n}({p}%)"),
            digits = all_continuous() ~ 2) |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "HF-TCM_period.docx")

#分类汇总
#替换A,B,Con都为0，C和D都是1
lipid_clin_coef2$分组 <- ifelse(lipid_clin_coef2$分组=="C",1,
                              ifelse(lipid_clin_coef2$分组=="D",1,0))

#统计出0组con-A-B和C-D组的差异
tbl_summary(lipid_clin_coef2,by = "分组",
            statistic = list(all_continuous()~"{median} ({p25}, {p75})",
                             all_categorical() ~ "{n}({p}%)"),
            digits = all_continuous() ~ 2) |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "noHF-HF_p.docx")

lipid_clin_coef2$分组 <- as.factor(lipid_clin_coef2$分组)

#批量单因素分析

uni_glm_result <- uni_log_test(lipid_clin_coef2,"分组")
cor_lipid <- cor(lipid_clin_coef2[,c(-1:-9)])
write.xlsx(lipid_clin_coef2,file = "lipid_clin_coef2.xlsx")

cor_lipd_data <- read.xlsx("lipid_clin_coef3.xlsx",1,header = F)
colnames(cor_lipd_data) <- cor_lipd_data[1,]
cor_lipd_data <- cor_lipd_data[-1,]
cor_lipd_data <- apply(cor_lipd_data,2,as.numeric)
cor_lipid_result <- cor(cor_lipd_data)
pmtcars <- cor_pmat(cor_lipid_result)

ggcorrplot(cor_lipid)
ggcorrplot(cor_lipid,method = "circle",lab=T)


pdf("cor_result_1.pdf",width = 10,height = 10)
ggcorrplot(cor_lipid_result,hc.order = T,  method = "circle",#分等级聚类重排矩阵
           ggtheme = ggplot2::theme_void(base_size = 10), #主题修改
           colors = c("CornflowerBlue","white","Salmon"), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
           lab = T,lab_size = 3,    #相关系数文本字体大小
           tl.cex = 15,             #坐标轴字体大小
           p.mat = pmtcars,         #添加显著性信息
           sig.level = 0.01,        #显著性水平
           pch = 5,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
           pch.cex = 8)            #不显著标记的大小，使用insig = "blank"将不显著的空白处理
dev.off()

colnames(lipid_clin_coef2) |> dput()
#需要校正的脂质或者评分
colnames(lipid_clin_coef2) |> dput()

adj_lip <- c("HSW","DCER24_1","DCER22_0", "HexCerd18_1_12_0", "DCER24_0", "PEP_18_0_20_4", 
             "HCER24_1", "PEO_18_0_20_4", "DCER16_0", "PE16_0_18_2", "PEO_18_0_22_5", 
             "PEP_18_0_22_5", "PEP_18_0_22_4", "PEO_18_0_22_4", "PEP_18_2_18_2", 
             "PE18_2_18_2", "HCER16_0", "PE16_0_20_2", "PEP_18_0_22_6", "coef")

adj_lip_list <- name_list(adj_lip)

#校正1，用年龄和性别
ad_mol1_index <- c("年龄","性别")

#校正2，用较正1+病史
ad_mol2_index <- c(ad_mol1_index,"冠心病","房颤","肾功能不全","高血压","糖尿病","高脂血症")

#校正3，用校正2+部分临床化验
ad_mol3_index <- c(ad_mol2_index, "P谷草转氨酶", 
                   "肌酸激酶同工酶", "P葡萄糖", "P谷丙转氨酶", "AST_ALT", 
                   "P白蛋白溴甲酚绿法", "P总蛋白", "球蛋白", "白球比值", 
                   "总胆红素", "直接胆红素", "间接胆红素", "总胆汁酸", 
                   "P甘油三酯", "P胆固醇", "P高密度脂蛋白", "P低密度脂蛋白", 
                   "碱性磷酸酶", "PY谷氨酰转移酶", "前白蛋白", "aL岩藻糖苷酶")


#批量单因素
uni_glm_result <- uni_log_test(lipid_clin_coef2,"分组")

#adj模型1
adj_glm_model_1 <- adj_log_glm(data=lipid_clin_coef2,
                               signle_list=adj_lip_list,
                               adj_immo=ad_mol1_index,
                               group_index = "分组",
                               outfile="ad_model_1_p.docx")

#adj模型2
adj_glm_model_2 <- adj_log_glm(data=lipid_clin_coef2,
                               signle_list=adj_lip_list,
                               adj_immo=ad_mol2_index,
                               group_index = "分组",
                               outfile="ad_model_2_p.docx")

#adj模型3
adj_glm_model_3 <- adj_log_glm(data=lipid_clin_coef2,
                               signle_list=adj_lip_list,
                               adj_immo=ad_mol3_index,
                               group_index = "分组",
                               outfile="ad_model_3_p.docx")

write.xlsx(uni_glm_result,file = "uni_glm_result_p.xlsx")


#统计C和D期的数据
tbl_summary(clini_data[,-1:-2],by = "中医证型",
            statistic = list(all_continuous()~"{median} ({p25}, {p75})",
                             all_categorical() ~ "{n}({p}%)"),
            digits = all_continuous() ~ 2) |> 
    add_p()|> 
    as_flex_table() |> 
    save_as_docx(path = "HF-TCM_period_v2.docx")
