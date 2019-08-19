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
source('functions.R')
source('colorPalette.R')

# 载入RPKM形式的表达矩阵。
load('../females_rpkm.Rdata')
# Extract embryonic stage data from cell names
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)
table(female_stages)

### 载入第一步挑选到的800多个高变异基因
load(file = '../step1-tSNE-female-RPKM/females_high_sv_matrix.Rdata')
dim(females_data)

## 载入第一步的tSNE后细胞类群信息
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]
table(female_clustering)

# 下面运行 DiffusionMap 需要822个基因的表达矩阵
# 以及它们被tSNE+DBSCAN聚的4类

# 可以把 DiffusionMap 类比 PCA 分析

# Compute the Diffusion map
female_dm <- run_diffMap(
  females_data, 
  female_clustering,
  sigma=15
)
 
save(female_dm,female_clustering,female_stages,
     file = 'diffusionMap_output.Rdata')
# Plot the Eigen values per diffusion component (similar to screeplots for PCAs)
plot_eigenVal(
  dm=female_dm
)
# 主要看拐点。

# Plot the 3D diffusion map using DCs 1 to 3 colored by the cell clustering
plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_clustering, 
  colour=female_clusterPalette
)

# Plot the 3D diffusion map using DCs 1 to 3 colored by embryonic stages
plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_stages, 
  colour=female_stagePalette
)


