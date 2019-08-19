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

# 加载上一步得到的4个细胞类群分布信息。
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]
table(female_clustering)
# 从文章里面拿到的基因列表。
markerGenes <- c(
  "Nr2f1",
  "Nr2f2",
  "Maf",
  "Foxl2",
  "Rspo1",
  "Lgr5",
  "Bmp2",
  "Runx1",
  "Amhr2",
  "Kitl",
  "Fst",
  "Esr2",
  "Amh",
  "Ptges"
)

# 首先把表达矩阵拆分，重新整理，方便绘图。
gene_subset <- as.matrix(log(females[rownames(females) %in% markerGenes,]+1))
gene_subset[1:4,1:4]
cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C1"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C2"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C3"])]
cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C4"])]

heatmap_gene_subset <- cbind(
  cl1_gene_subset, 
  cl2_gene_subset,
  cl3_gene_subset,
  cl4_gene_subset
)

heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]
heatmap_female_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)
table(heatmap_female_stages)
# 整理好的 heatmap_gene_subset 矩阵 后续绘制热图，共6个发育时期。
# 分成4组。

rowbreaks <- c(6, 15)

colbreaks <- c(
  ncol(cl1_gene_subset),
  ncol(cl1_gene_subset)+ncol(cl2_gene_subset), 
  ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset)
)

cluster_color <- c(
  C1="#560047",
  C2="#a53bad", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)

stage_color=c(
  E10.5="#2754b5", 
  E11.5="#8a00b0", 
  E12.5="#d20e0f", 
  E13.5="#f77f05", 
  E16.5="#f9db21",
  P6="#43f14b"
)


# 这里是自定义的热图，比较复杂，修改了pheatmap的非常多参数。 
library(pheatmap)
# tiff(file=paste0("../graph/", format(Sys.time(), "%Y%m%d"), "_female_de_genes_per_cell_clusters_small.tiff"), res = 300, height = 10, width = 16, units = 'cm')
png("female_marker.png")
plot_heatmap_2(
  heatmap_gene_subset, 
  female_clustering, 
  heatmap_female_stages, 
  rowbreaks, 
  colbreaks,
  cluster_color,
  stage_color
)
dev.off()
 


## 小提琴图显示表达量。
pdf("female_violin_markers.pdf", width=10, height=22)
require(gridExtra)

# Cluster color palette
female_clusterPalette <- c(
  "#560047", 
  "#a53bad", 
  "#eb6bac", 
  "#ffa8a0"
)
p <- list()
for (genes in markerGenes) {
  p[[genes]] <- violin_gene_exp(
    genes, 
    females, 
    female_clustering, 
    female_clusterPalette
  )
}

do.call(grid.arrange,c(p, ncol=3))

dev.off()

# t-SNE空间分布的表达量 等高线图。
 
pdf("female_tSNE_markers.pdf", width=16, height=28)
require(gridExtra)

load('../step1-tSNE-female-RPKM/step1-female_t_sne.Rdata')
p <- list()
for (genes in markerGenes) {
  p[[genes]] <- tsne_gene_exp(
    female_t_sne,
    genes, 
    females
  )
}

do.call(grid.arrange,c(p, ncol=4))

dev.off()
