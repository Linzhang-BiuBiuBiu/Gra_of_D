setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(gtsummary)          
library(flextable)
library(openxlsx)
library(mice)
source("E:\\R\\Zl_R_function\\zl_R_Function_v2.r")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#读取脂质数据以及临床数据
lipid_coef <- read.xlsx("clini_coef_data.xlsx")

clini_data <- read.xlsx("clin_data_v2.xlsx",2)

#数据合并，删除临床数据中的分组
lipid_clin_coef <- left_join(lipid_coef,clini_data[,-21])
lipid_clin_coef$HSW <- (as.numeric(lipid_clin_coef$红细胞分布宽度CV)/as.numeric(lipid_clin_coef$平均血红蛋白浓度))*100

#选择需要的列
colnames(lipid_clin_coef) |> dput()
clin_index <- c("分组",  "性别",  "冠心病", 
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
lipid_clin_coef1[,c(-1:-8)] <- apply(lipid_clin_coef1[,c(-1:-8)],2,as.numeric)

lipid_clin_coef1[,c(1:8)] <- apply(lipid_clin_coef1[,c(1:8)],2,as.factor)

#mice多重插补
mice_lipid_clin_coef1 <- mice(lipid_clin_coef1,m=1,seed = 5)

#读取多重插补信息,并替换列名
lipid_clin_coef2 <- complete(mice_lipid_clin_coef1,action = 1)
lipid_clin_coef2[,"HSW"] <- (lipid_clin_coef2$红细胞分布宽度CV/lipid_clin_coef2$平均血红蛋白浓度)*100
lipid_clin_coef2$coef <- lipid_clin_coef2$coef^2

#作图
library(ggpubr) # 继承ggplot语法
library(patchwork) # 拼图包
library(ggsci) #配色包

lipid_clin_coef2$分组 <- factor(lipid_clin_coef2$分组,levels = c("Con","A","B","C","D"))

sum_col <- c("分组", "性别", "冠心病", "房颤", "肾功能不全", 
             "高血压", "糖尿病", "高脂血症", "年龄","HSW", "BNP", "EF")

data <- lipid_clin_coef2

colnames(data)[colnames(data) %in% "HSW"] <- "EFI"

vio_box_jit_mulGrou(data=data,
                    group = "分组",
                    target = "年龄",
                    step=10,
                    size=20)

count(data$分组)
fivenum(data[,c("分组", "性别", "冠心病", "房颤", "肾功能不全", 
              "高血压", "糖尿病", "高脂血症")])

#统计证型数据，并逐步作图

lipid_clin_coef2$分组 <- factor(lipid_clin_coef2$分组,levels = c("Con","A","B","C","D"))

CHF_CD <- lipid_clin_coef2 |> filter(分组 %in% c("D","C"))

table(CHF_CD$中医证型)




