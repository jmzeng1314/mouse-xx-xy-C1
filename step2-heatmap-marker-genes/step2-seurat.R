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

# Extract embryonic stage data from cell names
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)
table(female_stages)
# Extract capture date from cell names
female_captures <- sapply(strsplit(colnames(females), "_"), `[`, 3)
female_captures <- paste(female_stages, female_captures, sep="_")
names(female_captures) <- colnames(females)
table(female_captures)

library(Seurat)
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
                      metadata = female_stages, col.name = 'female_stages')

# Cluster sce_xx
sce_xx <- ScaleData(object = sce_xx, 
                    vars.to.regress = c('nUMI_raw'), 
                    model.use = 'linear', 
                    use.umi = FALSE)
sce_xx <- FindVariableFeatures(object = sce_xx, 
                               mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, 
                               x.high.cutoff = 3, 
                               y.cutoff = 0.5)
length(VariableFeatures(sce_xx))

sce_xx <- RunPCA(object = sce_xx, pc.genes = VariableFeatures(sce_xx))
sce_xx <- RunTSNE(sce_xx )
sce_xx <- FindNeighbors(object = sce_xx, dims = 1:20, verbose = FALSE) 
sce_xx <- FindClusters(object = sce_xx,resolution = 0.3)
DimPlot(object = sce_xx, reduction = "tsne")

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
# Visualize canonical marker genes as violin plots.
pdf('seurat_VlnPlot.pdf', width=10, height=15)
VlnPlot(object = sce_xx, features =  markerGenes , 
        pt.size = 0.2,ncol = 4)
dev.off()

# Visualize canonical marker genes on the sctransform embedding.
pdf('seurat_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = sce_xx, features = markerGenes ,
            pt.size = 0.2,ncol = 3)
dev.off()

## 如何自己绘制小提琴图，等等
g=Idents(sce_xx)
mat=GetAssayData(object = sce_xx, slot = "scale.data")
mat[1:4,1:4]
ge=as.numeric(females['Nr2f2',])
ge=log(ge+1)
boxplot(ge~g)
library(ggpubr)
df=data.frame(value=ge,
              group=g)
ggviolin(df, "group", "value", fill = "group",
         # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))
## 作业，加上显著性统计指标


