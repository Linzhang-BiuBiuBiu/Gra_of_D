library(openxlsx)##读写excel
library(stringr)##加载包
library(officer)##读写word
library(magrittr)##增加管道操作符
library(eoffice)##编辑word

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

file_name <- grep("^\\d.*(xlsx)$",dir(),value = T)

#导出出院
for (k in file_name){
  temp <- tryCatch(
    {   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
      rx <- read.xlsx(k,2,detectDates=T)
      
      dur_data <- character()
      for(i in 1:nrow(rx)){
        for (j in 1:ncol(rx)) {
          dur_data <- paste0(dur_data,rx[i,j]) 
        }}
    
      dur_data <- gsub("NA|\\s","",dur_data)
      my_doc<-read_docx()
      my_doc<-my_doc %>% body_add_par(dur_data, style = "Normal")
      
      #步骤五：输出文档到本地
      setwd("D:\\onedrive\\桌面\\Clinical_data_summary\\入院+出院+手术-v2\\出院")
      filnames <- paste0(gsub("入院出院手术.xlsx","出院",k),".docx")
      print(my_doc, target = filnames )
    },
    warning = function(w) { return(NA) },
    error = function(e) { file_name[k] },
    finally = { file_name[k] }
  )
}
write.csv(file_name,file="excel_names.csv")
