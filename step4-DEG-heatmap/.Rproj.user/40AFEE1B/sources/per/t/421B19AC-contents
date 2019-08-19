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

## 加载一系列自定义函数及必备R包，及配色
source('functions.R')
# 加载转录组上游数据分析得到的表达矩阵
load('../females_rpkm.Rdata')

tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]  

## 差异分析，一般需要counts矩阵，使用 monocle 包
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


library(monocle)
# Prepare data to be loaded in Monocle
DE_female <- prepare_for_DE (
  female_count, 
  female_clustering, 
  female_stages
)

# Get genes diffet=rentially expressed between clusters (Monocle)
# 使用 monocle 做counts矩阵的差异分析通常耗时比较客观,至少一分钟。
# 1.12496 mins
if(F){
  start_time <- Sys.time()
  female_DE_genes <- findDEgenes(
    DE_female, 
    qvalue=0.05
  )
  end_time <- Sys.time()
  end_time - start_time
  save(female_DE_genes,file = 'female_DE_genes.Rdata')
}

load(file = 'female_DE_genes.Rdata')

# Get in which cluster the DE genes ar ethe most expressed (see analysis_functions.R for details)
load('../females_rpkm.Rdata')
# 后面可视化又开始使用 RPKM 值。
de_clusters <- get_up_reg_clusters(
  females, 
  female_clustering, 
  female_DE_genes
)

# Write the result in a file
write.csv(
  de_clusters, 
  quote = FALSE, 
  file= "female_DE_genes_per_clusters_4_groups.csv"
)
save(de_clusters,file = 'de_clusters.Rdata')


