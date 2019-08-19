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
library(monocle) 
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
# 提取embryonic stage信息
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)
table(female_stages)



## 首先构建 monocle 对象
library(monocle)
count_matrix=female_count
expr_matrix <- as.matrix(count_matrix)
sample_sheet <- data.frame(cells=names(count_matrix), 
                           stages=female_stages, 
                           cellType=female_clustering)
rownames(sample_sheet)<- names(count_matrix)
gene_annotation <- as.data.frame(rownames(count_matrix))
rownames(gene_annotation)<- rownames(count_matrix)
colnames(gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# Create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(
  as(expr_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit=0.5,
  expressionFamily=negbinomial.size()
)

HSMM <- detectGenes(HSMM, min_expr = 5)
# HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

cds=HSMM
# 单细胞转录组最重要的就是把细胞分群啦，这里可供选择的算法非常多，我们首先演示PCA结果。
# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
# plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 
plot_cell_clusters(cds, 1, 2, color = "cellType")
table(pData(cds)$Cluster,female_clustering)
plot_cell_clusters(cds, 1, 2 )


## 接着找差异基因

if(F){
  Sys.time()
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = "~cellType")
  Sys.time()
  # 可以看到运行耗时
  save(diff_test_res,file='diff_test_res.Rdata')
}
load(file = 'diff_test_res.Rdata')
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)
sig_genes$gene_short_name = rownames(sig_genes)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )

##  最后推断发育轨迹

## 首先挑选合适的基因
# 这里选取统计学显著的差异基因列表
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

# 然后降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# 降维是为了更好的展示数据。
# 降维有很多种方法, 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法

# 接着对细胞进行排序
cds <- orderCells(cds)

## 最后两个可视化函数 
plot_cell_trajectory(cds, color_by = "Cluster")  
plot_cell_trajectory(cds, color_by = "cellType")  
# 可以很明显看到细胞的发育轨迹

## 这里可以展现marker基因在发育轨迹推断的效果，本例子随便 选取了6个差异表达基因。
plot_genes_in_pseudotime(cds[head(sig_genes$gene_short_name),], 
                         color_by = "Cluster")






