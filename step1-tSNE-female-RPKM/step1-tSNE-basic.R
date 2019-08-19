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

# 提取embryonic stage信息
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
library(stringr)
female_stages=str_split(colnames(females), "_",simplify = T)[,1]
names(female_stages) <- colnames(females)
table(female_stages)
# 提取capture date 信息，可能没有用。
female_captures <- sapply(strsplit(colnames(females), "_"), `[`, 3)
female_captures <- paste(female_stages, female_captures, sep="_")
names(female_captures) <- colnames(females)
table(female_captures)


# 去除那些不在任何细胞表达的基因。
(dim(females))
females <- females[rowSums(females)>0,]
(dim(females))
# 去除了5千多个表达量恒定为0的基因
# 如果使用3大R包，这个步骤通常是被R包封装为参数或者函数。

## 首先探索基因的一些统计学指标，包括：mean,sd,mad,
exprSet=females
mean_per_gene <- apply(exprSet , 1, mean, na.rm = TRUE) #对表达矩阵每行求均值
sd_per_gene <- apply(exprSet, 1, sd, na.rm = TRUE) #对表达矩阵每行求标准差
mad_per_gene <-   apply(exprSet, 1, mad, na.rm = TRUE) #对表达矩阵每行求绝对中位差
cv = sd_per_gene/mean_per_gene

library(matrixStats)
var_per_gene <- rowVars(as.matrix(exprSet))
## 同样的apply函数，多次出现，请务必学透它！！！
cv2=var_per_gene/mean_per_gene^2
# 构造一个数据框来存放结果。
cv_per_gene <- data.frame(mean = mean_per_gene,
                          sd = sd_per_gene,
                          mad=mad_per_gene,
                          var=var_per_gene,
                          cv=cv,
                          cv2=cv2)
rownames(cv_per_gene) <- rownames(exprSet)
head(cv_per_gene)
# 进行表达量过滤。
cv_per_gene=cv_per_gene[cv_per_gene$mean>1,]

# pairs(cv_per_gene) 
with(cv_per_gene,plot(log10(mean),log10(cv2)))
# 可以继续美化

library(psych)
pairs.panels(cv_per_gene, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


# 这个函数摘抄自 Brennecke et al 2013 的文章：
# Extract genes with a squared coefficient of variation >2 times 
# the fit regression (Brennecke et al 2013 method)
# 下面的 getMostVarGenes 是自定义函数，在 functions.R 文件里面。
females_data <- getMostVarGenes(females, fitThr=2)
females_data <- log(females_data+1)
dim(females_data)
# 这里跟V3的seurat包得到的2000个var基因比较。
load(file = 'seurat_v3_2000_var_genes.Rdata')
length(intersect(rownames(females_data),seurat_v3_2000_var_genes))


# 挑选到的 822 个基因 后续做PCA分析，等等。
dim(females_data)
females_data[1:4,1:4]
save(females_data,file = 'females_high_sv_matrix.Rdata')
#pheatmap::pheatmap(females_data,show_rownames = F,show_colnames = F)

###########################################
#										  #
#			  RtSNE Analysis			  #
#										  #
###########################################

# Run a PCA on the highly variable genes
# 这里仅仅是 针对 822 个基因的表达矩阵做PCA分析。
female_sub_pca <- FactoMineR::PCA(
  t(females_data), 
  ncp = ncol(females_data), 
  graph=FALSE
)

# Estimate which PCs contain significant information with Jackstraw
# the result is 9 significant PCs
# 挑选显著的主成分数量，在seurat包有包装好的函数。
significant_pcs <- jackstraw::permutationPA(
  female_sub_pca$ind$coord, 
  B = 100, 
  threshold = 0.05, 
  verbose = TRUE, 
  seed = NULL
)$r
significant_pcs


# Compute and plot the t-SNE using the significant PCs
# 基于PCA的结果挑选显著的主成分进行tSNE
# 主要是因为 tSNE 耗费计算量。
female_t_sne <- run_plot_tSNE(
  pca=female_sub_pca,
  pc=significant_pcs,
  iter=5000,
  conditions=female_stages,
  colours=female_stagePalette
)
ggsave('tSNE_by_stage.pdf')
save(female_t_sne,file='step1-female_t_sne.Rdata')

# Recompute PCA with FactomineR package for the clustering
res.pca <- PCA(
  t(females_data), 
  ncp = significant_pcs, 
  graph=FALSE
)

# Clustering cells based on the PCA with 
# a minimum of 4 expected clusters
# Hierarchical Clustering On Principle Components (HCPC)
# 基于PCA的结果进行层次聚类，这里可以指定是4类。
res.hcpc <- HCPC(
  res.pca, 
  graph = FALSE,
  min=4
)

# Plot the hierarchical clustering result
plot(res.hcpc, choice ="tree", cex = 0.6)

# Extract clustering results
female_clustering <- res.hcpc$data.clust$clust
female_clustering <- paste("C", female_clustering, sep="")
names(female_clustering) <- rownames(res.hcpc$data.clust)
table(female_clustering)
# 作者这里想把调换两个元素的值，使用了非常丑陋的代码。
# Switch cluster names for convinience purpose (not very clean, quick and dirty method...)
female_clustering[female_clustering=="C1"] <- "C11"
female_clustering[female_clustering=="C2"] <- "C22"
female_clustering[female_clustering=="C22"] <- "C1"
female_clustering[female_clustering=="C11"] <- "C2"
table(female_clustering)
write.csv(female_clustering, file="female_clustering.csv")

# Cluster color palette
female_clusterPalette <- c(
  "#560047", 
  "#a53bad", 
  "#eb6bac", 
  "#ffa8a0"
)



# Plot the t-SNE colored by cell clusters
head(female_t_sne)
# 包装了一个ggplot绘图函数。
female_t_sne_new_clusters <- plot_tSNE(
  tsne=female_t_sne, 
  conditions=female_clustering, 
  colours= female_clusterPalette
)
ggsave('tSNE_cluster.pdf')


