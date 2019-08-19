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

# 载入表达矩阵以及diffusion map结果
load('../females_rpkm.Rdata')
load( file = 'diffusionMap_output.Rdata')

# 使用 Slingshot 来进行 pseudotime 分析，基于 diffusion map的主要成分
# 前面的图看的 DC3到DC4之间有一个拐点，所以我们考虑前面的4个diffusion map的主要成分
# We also help the pseudotime calculation by giving the starting cliuster and the end clusters
female_lineage <- get_lineage(
  dm=female_dm, 
  dim=c(1:4), 
  condition=factor(female_clustering),
  start="C1",
  end=c("C2", "C4"),
  shrink.method="cosine"
)
# 主要是 把  slingshot 函数的用法包装起来了。

# Plot the 3D diffusion map using DCs 1 to 3 colored by embryonic stages
# The black line is the cell lineage reconstruction by Slinshot
plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_stages, 
  colour=female_stagePalette
)
plot3d(female_lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file (first orient the graph as you want, the svg save does a snapshot of the graph the way you orient it)
rgl.postscript("female_dm_stages_cols.svg", fmt = "svg", drawText = TRUE )


# Plot the 3D diffusion map using DCs 1 to 3 colored by cell clusters
plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_clustering, 
  colour=female_clusterPalette
)
plot3d(female_lineage, dim=c(1:3), add=TRUE, lwd=5)
# Save as a svg file
rgl.postscript( "female_dm_clusters_cols.svg", fmt = "svg", drawText = TRUE )

# Attribute cells to lineages with a distance of 0.9 from the smoothed line and compute the pseudotime
female_pseudotime <- get_pseudotime(female_lineage, wthres=0.9)
rownames(female_pseudotime) <- colnames(females)
save(female_pseudotime,file = 'female_pseudotime.Rdata')
write.csv(female_pseudotime, file="female_pseudotime.csv")


