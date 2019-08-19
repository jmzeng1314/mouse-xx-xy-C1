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

load('../females_rpkm.Rdata')
load('../female_count.Rdata')
load('../step5-Slingshot/female_pseudotime.Rdata')
females[1:4,1:4]
female_count[1:4,1:4]
head(female_pseudotime)
tmp=read.csv('../step1-tSNE-female-RPKM/female_clustering.csv') 
female_clustering=tmp[,2];names(female_clustering)=tmp[,1]
table(female_clustering)

load(file = 'lineage_sig_gene.Rdata') 
load(file = 'gene_clustering.Rdata')
load(file = 'for_heatmap.Rdata')
table(cellLin)
table(cellType)
table(clustering[,1])
head(clustering)

# Load manually annotated gene clustering results (gene custers classified as a, b, c... according to their expression pattern redundandy, see paper fig. 2)
dyn_genes <- read.csv(file="../female_lineages_DE_gene_pseudotime_clustered_annotated.csv")
gene_names <- dyn_genes$Genes

library(clusterProfiler)
#convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# Subset genes with a correnc entrez ID
gene_clusters <- dyn_genes[dyn_genes$Genes %in% entrez_genes$SYMBOL,,drop=FALSE]

# Generate data frame for GoSemSim
de_gene_clusters <- data.frame(
  ENTREZID=entrez_genes[!duplicated(entrez_genes$SYMBOL),"ENTREZID"],
  Gene_Clusters=gene_clusters$Gene.categories
)

split(de_gene_clusters,de_gene_clusters$Gene_Clusters)

# Run enrichment analysis
formula_res <- compareCluster(
  ENTREZID~Gene_Clusters, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Run simplified GO enrochment analysis
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)

write.csv(formula_res@compareClusterResult, 
          file="female_compared_GO_term_DE_cluster.csv")
write.csv(lineage1_ego@compareClusterResult, 
          file="female_compared_symplified_GO_term_DE_cluster.csv")

pdf(file="female_GO_term_DE_genes_clusters.pdf", width=11, height=8)
dotplot(formula_res, showCategory=3)+ theme(aspect.ratio=0.8)
dotplot(lineage1_ego, showCategory=3)+ theme(aspect.ratio=2)
dev.off()


