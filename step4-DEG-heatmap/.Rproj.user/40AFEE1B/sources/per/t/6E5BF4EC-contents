## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2019-07-24 15:03:19
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2019-07-24  First version
###
### ---------------


rm(list=ls())
options(stringsAsFactors = F)
library(ROTS)
## 加载一系列自定义函数及必备R包，及配色
source('functions.R')
# 加载转录组上游数据分析得到的表达矩阵
load('../females_rpkm.Rdata')

tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]  

## 差异分析，一般需要counts矩阵, 这里演示  ROTS 包。
load('../female_count.Rdata')
# Extract embryonic stage data from cell names
head(colnames(female_count))
female_stages <- sapply(strsplit(colnames(female_count), "_"), `[`, 1)
names(female_stages) <- colnames(female_count)
table(female_stages)
# Extract capture date from cell names
female_captures <- sapply(strsplit(colnames(female_count), "_"), `[`, 3)
female_captures <- paste(female_stages, female_captures, sep="_")
names(female_captures) <- colnames(female_count)
table(female_captures)


###########################################
#                                         #
#         DEG by ROTS                     #
#                                         #
###########################################
## 分开做4次差异分析
## 运行速度非常慢，建议过夜等待，然后保存结果，方便重复使用。
 
library(ROTS)
library(plyr)
load('../females_rpkm.Rdata')
# 有趣的是，我看 ROTS 的一些例子，使用的是RPKM值。
RPKM.full=females
ROTS_input<-RPKM.full[rowMeans(RPKM.full)>=1,]
ROTS_input<-as.matrix(log2(ROTS_input+1))

if(!file.exists('ROTS_summary_pop.Rdata')){
  table(female_clustering)
  groups<-female_clustering
  groups<-as.numeric(substring(as.character(groups),2,2))
  table(groups)
  CAFgroups=groups
  groups[groups!=1]<-234
  
  Sys.time()
  results_pop1 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
  Sys.time()
  summary_pop1<-data.frame(summary(results_pop1, fdr=1))
  
  groups<-CAFgroups
  groups[groups!=2]<-134
  Sys.time()
  results_pop2 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
  Sys.time()
  summary_pop2<-data.frame(summary(results_pop2, fdr=1))
  
  
  groups<-CAFgroups
  groups[groups!=3]<-124
  Sys.time()
  results_pop3 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
  Sys.time()
  summary_pop3<-data.frame(summary(results_pop3, fdr=1))
  
  
  groups<-CAFgroups
  groups[groups!=4]<-123
  Sys.time()
  results_pop4 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
  Sys.time()
  summary_pop4<-data.frame(summary(results_pop4, fdr=1))
  
  save(summary_pop1,summary_pop2,summary_pop3,summary_pop4,
       file = 'ROTS_summary_pop.Rdata')
}

# 从 "2019-07-18 16:18:28 CST" 到  "2019-07-18 16:33:33 CST"
# 运行4次比较，耗时 45min，电脑电量耗尽！

## 以后直接载入 ROTS 的差异分析结果即可。
 
