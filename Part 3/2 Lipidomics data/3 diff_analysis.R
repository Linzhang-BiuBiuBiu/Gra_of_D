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

load("finaly_clin_lipid_list.rdata")
source("E:\\R\\Zl_R_function\\zl_R_Function.r")

#将QC和KB进行编组
for(i in names(finaly_clin_lipid_list)){
    tem_data <- finaly_clin_lipid_list[[i]]
    tem_data$SEQ <- ifelse(is.na(tem_data$SEQ),tem_data$ID,tem_data$SEQ)
    tem_data$ORGAN <- ifelse(is.na(tem_data$ORGAN),strsplit2(tem_data$SEQ,"-")[,2],tem_data$ORGAN)
    
    finaly_clin_lipid_list[[i]] <- tem_data
}

#首先做R，P,KB,QC的PCA
tem_data <- finaly_clin_lipid_list[[2]]
tem_data[is.na(tem_data)] <- 0

#删除空白KB
tem_data <- subset(tem_data,ORGAN!="KB")

group <- tem_data$ORGAN
tem_data <- tem_data[,c(-1:-5)]
tem_data <- apply(tem_data, 2, as.numeric)
tem_data <- t(apply(tem_data, 1, function(x) x/sum(x)))

data.pca=prcomp(tem_data)

pcaPredict=predict(data.pca)

PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=group)

main <- "PCA_RP_HF_combat"
pdf(file=paste0(main,".pdf"), width=5.5, height=4.25)
p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape=19, 
             ellipse=F, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
             size=2, main=main, legend="right")+
    theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
print(p1)
dev.off()


#提取分期的数据
tem_data <- finaly_clin_lipid_list[[4]]
tem_data[is.na(tem_data)] <- 0

#删除空白KB
tem_data <- subset(tem_data,ORGAN!="KB")
tem_data <- subset(tem_data,ORGAN!="QC")
tem_data <- subset(tem_data,ORGAN =="R")

#A-B,C-D统一分期
tem_data$分组 <- ifelse(tem_data$分组=="Con","Con",
                      ifelse(tem_data$分组=="A","AB",
                             ifelse(tem_data$分组=="B","AB",
                                                           ifelse(tem_data$分组=="C","CD","CD"))))
group <- tem_data$分组
tem_data <- tem_data[,c(-1:-6)]
tem_data <- apply(tem_data, 2, as.numeric)
tem_data <- t(apply(tem_data, 1, function(x) x/sum(x)))

data.pca=prcomp(tem_data)

pcaPredict=predict(data.pca)

PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=group)

main <- "PCA_HF_A-D_R_combat"
pdf(file=paste0(main,".pdf"), width=5.5, height=4.25)
p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape=19, 
             ellipse=F, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
             size=2, main=main, legend="right")+
    theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
print(p1)
dev.off()

#单独分析A-B和C-D组的数据，进行差异分析
#做OPLS-DA，只做A-B和C-D组的数据
tem_data <- finaly_clin_lipid_list[[4]]
tem_data[is.na(tem_data)] <- 0
#删除空白KB
tem_data <- subset(tem_data,ORGAN!="KB")
tem_data <- subset(tem_data,ORGAN!="QC")
tem_data <- subset(tem_data,ORGAN =="R")

#只提取第0天的数据
tem_data <- subset(tem_data,TIME =="0")
#A-B,C-D统一分期
tem_data$分组 <- ifelse(tem_data$分组=="Con","Con",
                      ifelse(tem_data$分组=="A","AB",
                             ifelse(tem_data$分组=="B","AB",
                                    ifelse(tem_data$分组=="C","CD","CD"))))

tem_data <- subset(tem_data,分组 !="Con")
group <- tem_data$分组

tem_data <- tem_data[,c(-1:-6)]
tem_data <- apply(tem_data, 2, as.numeric)
tem_data <- t(apply(tem_data, 1, function(x) x/sum(x)))
tem_data <- tem_data[,colMeans(tem_data)>0]
tem_data <- tem_data[,colSums(tem_data)>0]

pdf("A_B-C_D.oplsda.pdf")                   
df1_oplsda <- opls(tem_data, group, predI = 1,orthoI = NA)                 
dev.off() 

data <- as.data.frame(df1_oplsda@scoreMN)
o1 <- df1_oplsda@orthoScoreMN[,1]
data$o1 <- o1
data$Type = group
data$samples = rownames(data)

#提取解释度
x_lab <- df1_oplsda@modelDF[1, "R2X"] * 100
col=c("#1597A5","#FFC24B")
p1 <- ggplot(data,aes(x=p1,y=o1,color=Type))+#指定数据、X轴、Y轴，颜色
    theme_bw()+#主题设置
    geom_point(size=3)+#绘制点图并设定大小
    theme(panel.grid = element_blank())+
    geom_vline(xintercept = 0,lty="dashed",color="red")+
    geom_hline(yintercept = 0,lty="dashed",color="red")+#图中虚线
    # guides(color=guide_legend(title=NULL))+#去除图例标题
    labs(x=paste0("P1 (",x_lab,"%)"),
         y=paste0("to1"))+#将x、y轴标题改为贡献度
    stat_ellipse(data=data,
                 geom = "polygon",level = 0.95,
                 linetype = 2,size=0.5,
                 aes(fill=Type),
                 alpha=0.2,
                 show.legend = T)+
    scale_color_manual(values = col) +#点的颜色设置
    
    scale_fill_manual(values = c("#1597A5","#FFC24B"))+
    theme(axis.title.x=element_text(size = 12),#修改X轴标题文本
          axis.title.y=element_text(size =12,angle=90),#修改y轴标题文本
          axis.text.y=element_text(size =10),#修改x轴刻度标签文本
          axis.text.x=element_text(size =10),#修改y轴刻度标签文本
          panel.grid=element_blank())#隐藏网格线

pdf(file = paste0("A_B-C_D_nor.oplsda_v2.pdf"),width = 8,height = 6)
print(p1)
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
allDiff1 <- cbind(diff,FC_LOG[rownames(data_VIP),])
allDiff1 <- cbind(allDiff1,data_VIP[rownames(allDiff1),])
colnames(allDiff1)[length(allDiff1)] <- "VIP"

write.xlsx(cbind(rownames(allDiff1),allDiff1),file = "AB-CD_combat_diff.xlsx")

##聚类热图
load("all_data.RData")
heat_map_data <-cbind(as.data.frame(group),tem_data) 
heat_map_data2 <- heat_map_data[sample(1:1367,50),]
heat_map_data2 <- heat_map_data2[order(heat_map_data2$group),]

diffSig <- allDiff1 |> filter(P.Value<0.05&VIP>1)

group <- heat_map_data2$group

heat_map_data2 <- t(heat_map_data2[,rownames(diffSig)])
heat_map_data2[is.infinite(heat_map_data2)] <- 0

names(group)=colnames(heat_map_data2)
group=as.data.frame(group)

pdf(file="AB-CD_heatmap.pdf", width=10, height=100)

pheatmap(heat_map_data2, 
         annotation=group, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 16,
         fontsize_row=16,
         fontsize_col=16)

dev.off()

library(dplyr)
library(ggplot2)
library(ggrepel)

rt = read.table(inputFile, header=T, sep="\t", check.names=F)
#??????????
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#???ƻ?ɽͼ
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "black","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 40, hjust = 0.5, face = "bold"))
#???ڲ????????Ļ??򣬱?ע??????????
p1=p+geom_label_repel(data=filter(rt, ((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter))),
                      box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                      size=2.5, aes(label=id)) + theme_bw()
#??????ɽͼ
pdf(file="vol.pdf", width=6, height=6)
print(p1)
dev.off()
