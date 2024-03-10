setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(readr)
library(limma)
#读取离子数据
neg_file <- grep("^neg",dir(),value = T)
pos_file <- grep("^pos",dir(),value = T)

load("HF_lipid.rdata")

neg_file_number <- grep("^neg-33-|^neg-33.",dir(),value = T)
#neg_file_number <- grep("^neg-10-",dir(),value = T)

neg_file_number
for (i in 1:length(neg_file_number)) {
    tem_neg <- read.table(neg_file_number[i],sep = "\t",header = F)
    index_name <- gsub(".txt","",neg_file_number[i])
    
    tem_neg <- as.data.frame(tem_neg[,-2])
    
    tem_neg[tem_neg=="N/A"] <- 0
    
    #将英文统一成大写，_替换成-
    tem_neg[,1] <- toupper(gsub("(.*)\\:(.*)", "\\2", tem_neg[,1]))
    
    tem_neg[,1] <- sub("_","-",tem_neg[,1])
    tem_neg[,1] <- sub("_","-",tem_neg[,1])
    tem_neg[,1] <- sub("_","-",tem_neg[,1])
    
    #整理进样序列，分析进样批次与时间
    type <- strsplit2(tem_neg[,1],"-")
    
    if(i == 1){ 
    tem_neg[which(type[,3]=="QC"),1] <- paste(paste0(type[2,1],"-",type[2,2]),1:length(which(type[,3]=="QC")),sep="-QC-")
    tem_neg[which(type[,3]=="KB"),1] <- paste(paste0(type[2,1],"-",type[2,2]),1:length(which(type[,3]=="KB")),sep="-KB-")
    qc_max <- length(which(type[,3]=="QC"))
    kb_max <- length(which(type[,3]=="KB"))
    
    }else{
        qc_min_loc <- as.numeric(type[which(type[,3]=="QC")[1],4])
        kb_min_loc <- as.numeric(type[which(type[,3]=="KB")[1],4])
        
        tem_neg[which(type[,3]=="QC"),1] <- paste(paste0(type[2,1],"-",type[2,2]),qc_min_loc+qc_max:(length(which(type[,3]=="QC"))+qc_max-1),sep="-QC-")
        tem_neg[which(type[,3]=="KB"),1] <- paste(paste0(type[2,1],"-",type[2,2]),kb_min_loc+kb_max:(length(which(type[,3]=="KB"))+kb_max-1),sep="-KB-")
        
        qc_max <- length(which(type[,3]=="QC"))+qc_max
        kb_max <- length(which(type[,3]=="KB"))+kb_max
    }
   
    tem_neg <- t(tem_neg)
    rownames(tem_neg) <- tem_neg[,1]
    colnames(tem_neg) <- tem_neg[1,]
    tem_neg <- tem_neg[-1,-1]
    tem_rowname <- rownames(tem_neg)
    
    tem_neg <- apply(tem_neg, 2, as.numeric)
    rownames(tem_neg) <- tem_rowname
    neg_lipid_list[[index_name]] <- tem_neg
}

pos_lipid_list <- list()

pos_file <- grep("^pos",dir(),value = T)

#pos_file_number <- grep("^pos-14-|^pos-14.",dir(),value = T)
pos_file_number <- grep("^pos-1-",dir(),value = T)

pos_file_number
for (i in 1:length(pos_file_number)) {
    tem_pos <- read.table(pos_file_number[i],sep = "\t",header = F)
    index_name <- gsub(".txt","",pos_file_number[i])
    
    tem_pos <- as.data.frame(tem_pos[,-2])
    
    tem_pos[tem_pos=="N/A"] <- 0
    
    #将英文统一成大写，_替换成-
    tem_pos[,1] <- toupper(gsub("(.*)\\:(.*)", "\\2", tem_pos[,1]))
    
    tem_pos[,1] <- sub("_","-",tem_pos[,1])
    tem_pos[,1] <- sub("_","-",tem_pos[,1])
    tem_pos[,1] <- sub("_","-",tem_pos[,1])
    
    #整理进样序列，分析进样批次与时间
    type <- strsplit2(tem_pos[,1],"-")
    
    if(i == 1){ 
        tem_pos[which(type[,3]=="QC"),1] <- paste(paste0(type[2,1],"-",type[2,2]),1:length(which(type[,3]=="QC")),sep="-QC-")
        tem_pos[which(type[,3]=="KB"),1] <- paste(paste0(type[2,1],"-",type[2,2]),1:length(which(type[,3]=="KB")),sep="-KB-")
        qc_max <- length(which(type[,3]=="QC"))
        kb_max <- length(which(type[,3]=="KB"))
        
    }else{
        qc_min_loc <- as.numeric(type[which(type[,3]=="QC")[1],4])
        kb_min_loc <- as.numeric(type[which(type[,3]=="KB")[1],4])
        
        tem_pos[which(type[,3]=="QC"),1] <- paste(paste0(type[2,1],"-",type[2,2]),qc_min_loc+qc_max:(length(which(type[,3]=="QC"))+qc_max-1),sep="-QC-")
        tem_pos[which(type[,3]=="KB"),1] <- paste(paste0(type[2,1],"-",type[2,2]),kb_min_loc+kb_max:(length(which(type[,3]=="KB"))+kb_max-1),sep="-KB-")
        
        qc_max <- length(which(type[,3]=="QC"))+qc_max
        kb_max <- length(which(type[,3]=="KB"))+kb_max
    }
    
    tem_pos <- t(tem_pos)
    rownames(tem_pos) <- tem_pos[,1]
    colnames(tem_pos) <- tem_pos[1,]
    tem_pos <- tem_pos[-1,-1]
    tem_rowname <- rownames(tem_pos)
    
    tem_pos <- apply(tem_pos, 2, as.numeric)
    rownames(tem_pos) <- tem_rowname
    pos_lipid_list[[index_name]] <- tem_pos
}


save(pos_lipid_list,neg_lipid_list,file = "HF_lipid.rdata")