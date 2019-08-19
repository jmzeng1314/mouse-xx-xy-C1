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
library(pheatmap)
source('functions.R') 
# 载入rpkm表达量进行可视化，载入counts值找差异基因。
load('../females_rpkm.Rdata')
load('../female_count.Rdata')
# 载入上一步的细胞谱系分类 (2类) 
load('../step5-Slingshot/female_pseudotime.Rdata')
females[1:4,1:4]
female_count[1:4,1:4]
head(female_pseudotime)
## 载入细胞的(PCA+tSNE+DBSCAN) 4类
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]
table(female_clustering)

## 载入 细胞谱系分类，细胞的聚类，基因的分类：
load(file = 'lineage_sig_gene.Rdata') 
load(file = 'gene_clustering.Rdata')
load(file = 'for_heatmap.Rdata')
table(cellLin)
table(cellType)
table(clustering[,1])

# Palette for the heatmap (pheatmap does not take the color palette into account, I don't know why but I kept the code)
gene_cluster_palette <- c(
  '#a6cee3',
  '#1f78b4',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#ffff99',
  '#b15928', 
  '#49beaa', 
  '#611c35', 
  '#2708a0',
  '#fccde5',
  '#bc80bd'
)
gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
names(gene_cluster_colors) <- 1:max(clusters)


# Annotate the row of the hearmap with the gene clusters
annotation_row <- data.frame(clustering=clustering)

# Annotation of the colums of the heatmap with the cell lineage, the cell type and the embryonic stage of the cells
annotation_col <- data.frame(
  cellLineages=cellLin,
  cellType=cellType,
  Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)
rownames(annotation_col) <- colnames(data_heatmap)

# Colors for the cell lineages
cellLinCol <- c(
  "#3b3561", 
  "#c8c8c8", 
  "#ff6663"
)
names(cellLinCol) <- unique(cellLin)

# Colors for the cell clusters
cellTypeCol <- c(
  C2="#a53bad", 
  C1="#560047", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)
names(cellTypeCol) <- unique(cellType)

# Legend for all the aannotation to plot on pheatmap
annotation_colors <- list(
  cellType=cellTypeCol,
  cellLineages=cellLinCol,
  clustering=gene_cluster_colors,
  Stages=c(
    E10.5="#2754b5", 
    E11.5="#8a00b0", 
    E12.5="#d20e0f", 
    E13.5="#f77f05", 
    E16.5="#f9db21",
    P6="#43f14b"
  )
)

# Color palette for gene expression
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
breaksList = seq(-2.2, 2.5, by = 0.2)

# Plot the heatmap
library(pheatmap)
tiff(file="female_heatmap_DE_genes_granulosa_progenitors_k_17_pval_005.tiff", 
res = 300, height = 21, width = 18, units = 'cm')
gene_clustering <- pheatmap(
  data_heatmap, 
  scale="row",
  gaps_col=length(cellType_L2),
  show_colnames=FALSE, 
  show_rownames=FALSE, 
  cluster_cols=FALSE,
  clustering_method="ward.D",
  annotation_row=annotation_row,
  annotation_col=annotation_col,
  annotation_colors=annotation_colors,
  cutree_rows=17, 
  annotation_names_row=FALSE,
  color=mypalette
)
dev.off()



