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
names(female_stages) <- colnames(females)
table(female_stages)
# 提取capture date 信息，可能没有用。
female_captures <- sapply(strsplit(colnames(females), "_"), `[`, 3)
female_captures <- paste(female_stages, female_captures, sep="_")
names(female_captures) <- colnames(females)
table(female_captures)


# 下面的代码仅仅是对 Seurat(V3) 有效
# Create Seurat object
sce_xx <- CreateSeuratObject(counts = females, 
                             min.cells = 1, min.features = 0, 
                             project = 'c1_xy') # already normalized
sce_xx  

# Add meta.data (nUMI and timePoints)
sce_xx <- AddMetaData(object = sce_xx, 
                      metadata = apply(females, 2, sum), 
                      col.name = 'nUMI_raw')
sce_xx <- AddMetaData(object = sce_xx, 
                      metadata = female_stages, 
                      col.name = 'female_stages')

# 步骤 ScaleData 的耗时取决于电脑系统配置（保守估计大于一分钟）
sce_xx <- ScaleData(object = sce_xx, 
                    vars.to.regress = c('nUMI_raw'), 
                    model.use = 'linear', 
                    use.umi = FALSE)
sce_xx <- FindVariableFeatures(object = sce_xx, 
                               mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, 
                               x.high.cutoff = 4, 
                               y.cutoff = 0.5)
length(VariableFeatures(sce_xx))
seurat_v3_2000_var_genes=VariableFeatures(sce_xx)
save(seurat_v3_2000_var_genes,file = 'seurat_v3_2000_var_genes.Rdata')
sce_xx <- RunPCA(object = sce_xx, pc.genes = VariableFeatures(sce_xx))
# 下面只是展现不同降维算法而已，并不要求都使用一次。
sce_xx <- RunICA(sce_xx )
sce_xx <- RunTSNE(sce_xx )
#sce_xx <- RunUMAP(sce_xx,dims = 1:10)
#VizPCA( sce_xx, pcs.use = 1:2)
DimPlot(object = sce_xx, reduction = "pca") 
DimPlot(object = sce_xx, reduction = "ica")
DimPlot(object = sce_xx, reduction = "tsne")
#DimPlot(object = sce_xx, reduction = "umap")

# 针对PCA降维后的表达矩阵进行聚类 FindNeighbors+FindClusters 两个步骤。
sce_xx <- FindNeighbors(object = sce_xx, dims = 1:20, verbose = FALSE) 
sce_xx <- FindClusters(object = sce_xx, verbose = FALSE)
# 继续tSNE可视化
DimPlot(object = sce_xx, reduction = "tsne")
DimPlot(object = sce_xx, reduction = "tsne",
        group.by = 'female_stages')

clu1=as.data.frame(Idents(sce_xx))
clu2=read.csv('female_clustering.csv')
identical(rownames(clu1),clu2[,1])
table(clu1[,1],clu2[,2])
group4=as.character(clu2[,2])
sce_xx <- AddMetaData(object = sce_xx, 
                      metadata = group4, col.name = 'group4')
DimPlot(object = sce_xx, reduction = "tsne")
DimPlot(object = sce_xx, reduction = "tsne",
        group.by = 'group4')

# 聚类的时候可以控制 resolution 参数，调整细胞亚群的数量。
sce_xx <- FindClusters(object = sce_xx,resolution = 0.3)
DimPlot(object = sce_xx, reduction = "tsne")
clu1=as.data.frame(Idents(sce_xx))
clu2=read.csv('female_clustering.csv')
identical(rownames(clu1),clu2[,1])
table(clu1[,1],clu2[,2])

# 可以看到聚类并不是完全一致的，但是还算是比较稳定。



