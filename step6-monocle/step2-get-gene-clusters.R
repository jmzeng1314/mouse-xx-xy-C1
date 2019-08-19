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
# 主要是使用pheatmap对表达矩阵的基因进行聚类，17类
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
# 载入上一步的细胞谱系信息，分2类：
load('../step5-Slingshot/female_pseudotime.Rdata')
females[1:4,1:4]
female_count[1:4,1:4]
head(female_pseudotime)
## 载入细胞的(PCA+tSNE+DBSCAN) 4类
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]
table(female_clustering)

load(file = 'lineage_sig_gene.Rdata') 

head(female_lineage1_sig_gene_pseudoT)
head(female_lineage2_sig_gene_pseudoT)
#  为了保证后面的代码不改名，这里修一下
female_lineage1_clustering <- female_lineage1_sig_gene_pseudoT
female_lineage2_clustering <- female_lineage2_sig_gene_pseudoT


# Get gene names
gene_list <- unique(rownames(female_lineage1_clustering), 
                    rownames(female_lineage2_clustering))
length(gene_list)
dim(females)
# Get RPKM matrix of the dynamic genes
de_matrix <- log(females[rownames(females) %in% gene_list,]+1)
dim(de_matrix)

# Get cells attributed to the lineage 1
L1_lineage <- female_pseudotime[!is.na(female_pseudotime[,1]),1]
# Order the cells through the pseudotime
L1_ordered_lineage <- L1_lineage[order(L1_lineage, decreasing = FALSE)]
# Get expression values of the dynamic genes in the cells from lineage 1
L1_rpkm_exp_lineage <- de_matrix[,names(L1_ordered_lineage)]
# Order the cells (colums of the expression matrix) through the pseudotime
L1_cells <- L1_rpkm_exp_lineage[,order(match(names(L1_ordered_lineage), 
                                             L1_rpkm_exp_lineage))]

# Same as previously but for the linage 2 and but cells are ordered decreasingly (to be ploted on the left side)
L2_lineage <- female_pseudotime[!is.na(female_pseudotime[,2]),2]
L2_ordered_lineage <- L2_lineage[order(L2_lineage, decreasing = TRUE)]
L2_rpkm_exp_lineage <- de_matrix[,names(L2_ordered_lineage)]
L2_cells <- L2_rpkm_exp_lineage[,order(match(names(L2_ordered_lineage), 
                                             L2_rpkm_exp_lineage))]

# Extract cell names
L1_lineage_cells <- names(L1_ordered_lineage)
L2_lineage_cells <- names(L2_ordered_lineage)

# Compare nales in the two lisis
common_cells <- intersect(L1_lineage_cells, L2_lineage_cells)
L1_spe_cells <- setdiff(L1_lineage_cells,common_cells)
L2_spe_cells <- setdiff(L2_lineage_cells,common_cells)

# label the cells if they are part of the common branch or only in the lineage 1 branch
L1_cellLin <- c(rep("common cells",length(common_cells)),
                rep("L1 cells", length(L1_spe_cells))
)
 
names(L1_cellLin) <- c(common_cells, L1_spe_cells)
L1_cellLin <- L1_cellLin[order(match(names(L1_cellLin), colnames(L1_cells)))]

# label the cells if they are part of the common branch or only in the lineage 2 branch
L2_cellLin <- c(rep("common cells",length(common_cells)),
                rep("L2 cells", length(L2_spe_cells))
)

names(L2_cellLin) <- c(common_cells, L2_spe_cells)
L2_cellLin <- L2_cellLin[order(match(names(L2_cellLin), colnames(L2_cells)))]

# Get the cell cluster of the cells for each cell lineage
cellType_L1 <- female_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)
cellType_L2 <- female_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)

# Merge cell lineages annotation
cellLin <- c(
  L2_cellLin,
  L1_cellLin
)

# Merge the cell type annotation
cellType <- c(
  cellType_L2,
  cellType_L1
)

# Smooth the gene expression (Loess regression)
L2_cells_smooth <- smooth_gene_exp(
  L2_cells, 
  L2_ordered_lineage, 
  span=0.4
)
L1_cells_smooth <- smooth_gene_exp(
  L1_cells, 
  L1_ordered_lineage, 
  span=0.4
)

# Merge the smoothed expression matrices
data_heatmap <- data.frame(
  L2_cells_smooth,
  L1_cells_smooth
)

# Run pheatmap to cluster the genes per expression pattern
set.seed(123)
gene_clustering <- pheatmap::pheatmap(
  data_heatmap, 
  scale="row", 
  clustering_method="ward.D",
  silent=TRUE
)

# Extract the gene pattern clusters with k-means (17 expected clusters)
clusters <- cutree(gene_clustering$tree_row, k = 17)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
head(clustering)
table(clustering[,1])
# Save clustering results
write.csv(clustering, 
          file="female_mean_norm_female_lineage2_lineage1_DE_gene_pseudotime_qval_0.05_gene_clustering_kmeans_k17_scaled.csv")
save(clusters,clustering,file = 'gene_clustering.Rdata')
save(data_heatmap,cellLin,cellType,cellType_L2,
     file = 'for_heatmap.Rdata')
