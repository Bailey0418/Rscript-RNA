#####E-MTAB-1980 data#####
setwd("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC预后模型检验/E-MTAB-1980")
library(tidyverse)
path <- "D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC预后模型检验/E-MTAB-1980/data"
files <- list.files(path, pattern = "ccRCC-.*-tumor_expression\\.txt$", full.names = TRUE)
read_agilent <- function(file) {
  # 找到 "FEATURES" 所在行
  lines <- readLines(file)
  start_line <- grep("^FEATURES", lines)
  
  # 没找到就跳过
  if (length(start_line) == 0) return(NULL)
  
  # 从FEATURES那行开始读表
  dat <- read.delim(file, skip = start_line - 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 只保留探针名和表达值
  dat <- dat %>%
    select(ProbeName, gProcessedSignal) %>%
    rename(Expression = gProcessedSignal)
  
  # 添加样本名
  dat$Sample <- basename(file)
  
  return(dat)
}

expr_list <- lapply(files, read_agilent)
expr_list <- expr_list[!sapply(expr_list, is.null)] 
expr_mat <- do.call(cbind, lapply(expr_list, function(df) df$Expression))
rownames(expr_mat) <- expr_list[[1]]$ProbeName
colnames(expr_mat) <- basename(files)