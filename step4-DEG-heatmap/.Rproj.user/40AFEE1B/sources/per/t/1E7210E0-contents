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
library(pheatmap)

# 载入表达矩阵及表型信息
load('../females_rpkm.Rdata')
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]  
  
source('functions.R')

## 载入 monocle 的差异分析结果。
load(file = '../step3-DEG-annotation/de_clusters.Rdata')
# DE genes
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.0001)
gene_names <- gene_names$genes

gene_names <- get_top_up_reg_clusters(de_clusters, 20)
gene_subset <- as.matrix(log(females[rownames(females) %in% gene_names,]+1))

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

heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]
heatmap_female_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)
table(heatmap_female_stages)
female_stages=heatmap_female_stages
table(female_clustering)
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

tiff(file="female_de_genes_per_cell_clusters_small.tiff", 
res = 300, height = 10, width = 16, units = 'cm')
library(pheatmap)
plot_heatmap_2(
  heatmap_gene_subset, 
  female_clustering, 
  female_stages, 
  rowbreaks, 
  colbreaks,
  cluster_color,
  stage_color
)
dev.off()


