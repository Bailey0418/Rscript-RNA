options(stringsAsFactors = F) 
library(tidyverse)

a1 <- read.csv("/mnt/raid5/User/bailin/project/241209Constipation_lncRNA/result/01Data/mRNA_count.csv",row.names = 1)
### counts矩阵的构建
counts <- a1[,7:ncol(a1)] #截取样本基因表达量的counts部分作为counts 
rownames(counts) <- a1$Geneid #将基因名作为行名

### 从featurecounts 原始输出文件counts.txt中提取Geneid、Length(转录本长度)，
geneid_efflen <- subset(a1,select = c("Geneid","Length"))
colnames(geneid_efflen) <- c("geneid","efflen")  
geneid_efflen_fc <- geneid_efflen #用于之后比较

### 取出counts中geneid的对应的efflen
dim(geneid_efflen)
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]

### 计算 TPM
#TPM (Transcripts Per Kilobase Million)  每千个碱基的转录每百万映射读取的Transcripts
counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)       #每千碱基reads (“per million” scaling factor) 长度标准化
  PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
  RPK/PMSC_rpk                    
}  
tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm)

### 计算FPKM
#FPKM/RPKM (Fragments/Reads Per Kilobase Million )  每千个碱基的转录每百万映射读取的Fragments/reads
#RPKM与FPKM分别针对单端与双端测序而言，计算公式是一样的
counts2FPKM <- function(count=count, efflength=efflen){ 
  PMSC_counts <- sum(count)/1e6   #counts的每百万缩放因子 (“per million” scaling factor) 深度标准化
  FPM <- count/PMSC_counts        #每百万reads/Fragments (Reads/Fragments Per Million) 长度标准化
  FPM/(efflength/1000)                                    
}

#FPKM与TPM的转化
FPKM2TPM <- function(fpkm){
  fpkm/sum(fpkm)*1e6
}

fpkm <- as.data.frame(apply(counts,2,counts2FPKM))
colSums(fpkm)


#FPKM与TPM的转化
FPKM2TPM <- function(fpkm){
  fpkm/sum(fpkm)*1e6
}

tpm0 <- as.data.frame(apply(counts,2,FPKM2TPM))
colSums(tpm0)
