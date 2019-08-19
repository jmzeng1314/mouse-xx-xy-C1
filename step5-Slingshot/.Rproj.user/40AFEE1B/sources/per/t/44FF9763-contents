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

# 载入表达矩阵以及diffusion map结果，还有 Slingshot 结果
load('../females_rpkm.Rdata')
load(file = 'diffusionMap_output.Rdata')
load(file = 'female_pseudotime.Rdata')

# Make the pseudotime of each lineage between 0 and 100 (percentage of cell progression over the pseudotime)
# thos makes the two cell lineage progression comparable
pseudotime_lin <- female_pseudotime[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime

pseudotime_lin <- female_pseudotime[,"curve2"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin2_percent <- (pseudotime_lin*100)/max_pseudotime

female_pseudotime[,"curve1"] <- pseudotime_lin1_percent
female_pseudotime[,"curve2"] <- pseudotime_lin2_percent

# Lineage 相当于 细胞新的表型信息，跟之前的细胞类群，细胞取样的胚胎时期等同。

## Lineage result: 
# Lineage 1 -> progenitors to Granulisa
# Lineage 2 -> progenitors to progenitors

####################################################
#
# Plot gene expression along pseudotime in cell lineages 1 and 2
#
####################################################

# Twisted palette to plot gene expression by lineages along pseudotime
female_clusterPalette2 <- c(
  "#ff6663", 
  "#3b3561"
)

# Plot one gene expression over pseudotime for the lineage 1
plot_smoothed_gene_per_lineage(
  rpkm_matrix=females, 
  pseudotime=female_pseudotime, 
  lin=c(1),
  gene="Amhr2", 
  stages=female_stages, 
  clusters=female_clustering, 
  stage_colors=female_stagePalette,
  cluster_colors=female_clusterPalette,
  lineage_colors=female_clusterPalette2
)

# Plot one gene expression over pseudotime for the lineage 2
plot_smoothed_gene_per_lineage(
  rpkm_matrix=females, 
  pseudotime=female_pseudotime, 
  lin=c(2),
  gene="Amhr2", 
  stages=female_stages, 
  clusters=female_clustering, 
  stage_colors=female_stagePalette,
  cluster_colors=female_clusterPalette,
  lineage_colors=female_clusterPalette2
)

# Plot one gene expression over pseudotime for the lineage 1 and 2
plot_smoothed_gene_per_lineage(
  rpkm_matrix=females, 
  pseudotime=female_pseudotime, 
  lin=c(1,2),
  gene="Amhr2", 
  stages=female_stages, 
  clusters=female_clustering, 
  stage_colors=female_stagePalette,
  cluster_colors=female_clusterPalette,
  lineage_colors=female_clusterPalette2
)


# List of genes of interest to plot lineage by lineage through pseudotime
gene_list <- c(
  "Sall4",
  "Sox11",
  "Gata4",
  "Lgr5",
  "Runx1",
  "Foxl2",
  "Hey2",
  "Wnt5a",
  "Pdgfra",
  "Nr2f2",
  "Sfrp1",
  "Ifitm3",
  "Ptch1",
  "Wnt4",
  "Rspo1",
  "Cdkn1b",
  "Gli1",
  "Tcf21",
  "Nr0b1",
  "Nr0b2",
  "Nr5a1",
  "Nr6a1"
)

plot_smoothed_genes <- function(genes, lin){
  female_clusterPalette2 <- c("#ff6663", "#3b3561")
  for (gene in genes){
    plot_smoothed_gene_per_lineage(
      rpkm_matrix=females, 
      pseudotime=female_pseudotime, 
      lin=lin,
      gene=gene, 
      stages=female_stages, 
      clusters=female_clustering, 
      stage_colors=female_stagePalette,
      cluster_colors=female_clusterPalette,
      lineage_colors=female_clusterPalette2
    )
  }
}

pdf("female_interesting_profiles.pdf", width=4, height=4)
plot_smoothed_genes(gene_list, 1) # plot only lineage 1
plot_smoothed_genes(gene_list, 2) # plot only lineage 2
plot_smoothed_genes(gene_list, c(1,2)) # plot the two moleages in the same graph to see the divergence
dev.off()



