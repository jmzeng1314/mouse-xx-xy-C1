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
library(Seurat)
library(monocle)

# 载入表达矩阵及表型信息
load('../females_rpkm.Rdata') 
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]  
RPKM.full=females

## 载入上一步 ROTS 的差异分析结果
load(file = '../step3-DEG-annotation/ROTS_summary_pop.Rdata')
head(summary_pop1)
head(summary_pop2)
head(summary_pop3)
head(summary_pop4)



## 每个亚群挑选top18的差异基因绘制热图

population_subset<-c(rownames(summary_pop1[summary_pop1$ROTS.statistic<0,])[1:18],
                     rownames(summary_pop2[summary_pop2$ROTS.statistic<0,])[1:18],
                     rownames(summary_pop3[summary_pop3$ROTS.statistic<0,])[1:18],
                     rownames(summary_pop4[summary_pop4$ROTS.statistic<0,])[1:18])

RPKM_heatmap<-RPKM.full[population_subset,]

RPKM_heatmap<-RPKM_heatmap[,order(female_clustering)]
RPKM_heatmap<-log2(RPKM_heatmap+1)

popul.col<-sort(female_clustering);table(popul.col)
popul.col<-replace(popul.col, popul.col=='C1',"#1C86EE" )
popul.col<-replace(popul.col, popul.col=='C2',"#00EE00" )
popul.col<-replace(popul.col, popul.col=='C3',"#FF9912" )
popul.col<-replace(popul.col, popul.col=='C4',"#FF3E96" )
library(gplots)

pdf("heatmap_DEG_ROST.pdf")
heatmap.2(as.matrix(RPKM_heatmap),
          ColSideColors = as.character(popul.col), tracecol = NA, 
          dendrogram = "none",col=bluered, labCol = FALSE, 
          scale="none", key = TRUE, symkey = F, symm=F,  key.xlab = "", 
          key.ylab = "", density.info = "density", 
          key.title = "log2(RPKM+1)", keysize = 1.2, denscol="black", Colv=FALSE)
dev.off()







