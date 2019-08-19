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

# 安装及加载必备的R包，耗时至少 1 小时+

# 这里需要 Seurat 是 3.0 的，如果版本不够，需要升级
if(F){
  library(Seurat)
  remove.packages('Seurat')
  install.packages('Seurat')
  library(Seurat)
}

library(ROTS)
library(Seurat)
library("matrixStats")
library("ggplot2")
library("Rtsne")
library("fpc")
library("factoextra")
library("monocle")
library("viridis")
library("gplots")
library("RColorBrewer")
library("destiny")
library("slingshot")
library("rgl")
library("scatterplot3d")
library("made4")
library("pheatmap")
library("matrixStats")
library("statmod")

library("FactoMineR")
library("jackstraw")
 
library("org.Mm.eg.db")
library("clusterProfiler")
library("GOSemSim")
library("arulesViz")

library("ggpubr")




options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror
if(F){
  #BiocManager::install("Seurat",ask = F,update = F)
  BiocManager::install("monocle",ask = F,update = F) 
  BiocManager::install("destiny",ask = F,update = F)
  BiocManager::install("slingshot",ask = F,update = F)
  BiocManager::install("made4",ask = F,update = F)
  BiocManager::install("lfa",ask = F,update = F)
  BiocManager::install("ROTS",ask = F,update = F)
}


pkgs = c("taRifx"
         ,"matrixStats"
         ,"ggplot2"
         ,"Rtsne"
         ,"fpc"
         ,"factoextra"
         ,"monocle"
         ,"viridis"
         ,"gplots"
         ,"RColorBrewer"
         ,"destiny"
         ,"slingshot"
         ,"rgl"
         ,"scatterplot3d"
         ,"made4"
         ,"pheatmap"
         ,"matrixStats"
         ,"statmod"
         ,"FactoMineR"
         ,"jackstraw"
         ,"ReactomePA"
         ,"org.Mm.eg.db"
         ,"clusterProfiler"
         ,"GOSemSim"
         ,"arulesViz"
         ,"ggpubr")

# BiocManager::install(pkgs,ask = F,update = F)

###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

library("taRifx")
library("matrixStats")
library("ggplot2")
library("Rtsne")
library("fpc")
library("factoextra")
library("monocle")
library("viridis")
library("gplots")
library("RColorBrewer")
library("destiny")
library("slingshot")
library("rgl")
library("scatterplot3d")
library("made4")
library("pheatmap")
library("matrixStats")
library("statmod")

library("FactoMineR")
library("jackstraw")

library("ReactomePA")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GOSemSim")
library("arulesViz")

library("ggpubr")

