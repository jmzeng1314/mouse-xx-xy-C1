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
library(Seurat)
library(monocle)
library(ROTS)
library(clusterProfiler)

load(file = 'de_clusters.Rdata')
load(file = 'ROTS_summary_pop.Rdata')
head(summary_pop1)
head(summary_pop2)
head(summary_pop3)
head(summary_pop4)



# Extract DE gene names
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.05)
gene_names <- gene_names$genes

# Convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]

de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,
                             c("genes", "cluster")]
de_gene_clusters <- data.frame(
  ENTREZID=entrez_genes$ENTREZID[entrez_genes$SYMBOL %in% de_gene_clusters$genes],
  cluster=de_gene_clusters$cluster
)
table(de_gene_clusters$cluster)
list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, 
                               de_gene_clusters$cluster)


# Run full GO enrichment test
formula_res <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)

# Plot both analysis results
pdf('female_compared_GO_term_DE_cluster.pdf',width = 11,height = 6)
dotplot(formula_res, showCategory=5)
dev.off()
pdf('female_compared_GO_term_DE_cluster_simplified.pdf',width = 11,height = 6)
dotplot(lineage1_ego, showCategory=5)
dev.off()


# Save results
write.csv(formula_res@compareClusterResult, 
          file="female_compared_GO_term_DE_cluster.csv")
write.csv(lineage1_ego@compareClusterResult, 
          file="female_compared_GO_term_DE_cluster_simplified.csv")

## 扩展：enric其它数据库，KEGG, RECTOME, GO(CC,MF) 

library(GO.db)
ls("package:GO.db")
library(org.Mm.eg.db)
go2gene=toTable(org.Mm.egGO2ALLEGS)
this_go_this_gene=go2gene[go2gene$go_id=='GO:0140014',]
table(this_go_this_gene$Evidence)
# this_go_this_gene=this_go_this_gene[this_go_this_gene$Evidence %in% c(),]
this_go_this_gene=this_go_this_gene[,1]
length(this_go_this_gene)
c1_diff_gene=list_de_gene_clusters[['C1']]
length(c1_diff_gene)
length(intersect(this_go_this_gene,c1_diff_gene))




