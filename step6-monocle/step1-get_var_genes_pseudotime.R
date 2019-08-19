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
library(data.table)
library(densityClust)
library(destiny)
library(rgl)
library(monocle)
source('functions.R') 

# 载入rpkm表达量进行可视化，载入counts值找差异基因。
load('../females_rpkm.Rdata')
load('../female_count.Rdata')
# 载入上一步的细胞谱系信息，分2类：
load('../step5-Slingshot/female_pseudotime.Rdata')
females[1:4,1:4]
female_count[1:4,1:4]
head(female_pseudotime)
# 细胞谱系需要归一化成为百分比。

## 载入细胞的(PCA+tSNE+DBSCAN) 4类
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]
table(female_clustering)

######################################################
#										             #
#   Genes dynamically expressed through pseudotime   #
#										             #
######################################################

# 使用 Monocle 找到那些随着不同细胞谱系变化剧烈的基因
# 因为是 Monocle 所以耗时有点可观。
# 对 female_pseudotime 的两列分别进行差异分析

female_lineage1_sig_gene_pseudoT  <- get_var_genes_pseudotime(
  females, 
  female_count, 
  female_pseudotime, 
  lineageNb=1, 
  female_clustering
)



# 从12612个基因里面挑选出2858个差异显著的基因
female_lineage1_sig_gene_pseudoT <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
write.csv(female_lineage1_sig_gene_pseudoT, 
          file= "female_lineage1_pseudotime_DE_genes.csv")

#############################################################
female_lineage2_sig_gene_pseudoT <- get_var_genes_pseudotime(
  females, 
  female_count, 
  female_pseudotime, 
  lineageNb=2, 
  female_clustering
)
# 从11937个基因里面挑选出2182个差异显著的基因
female_lineage2_sig_gene_pseudoT <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]
write.csv(female_lineage2_sig_gene_pseudoT, 
          file="female_lineage2_pseudotime_DE_genes.csv")

save(female_lineage1_sig_gene_pseudoT,
     female_lineage2_sig_gene_pseudoT,
     file = 'lineage_sig_gene.Rdata')

