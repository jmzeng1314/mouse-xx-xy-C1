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


load(file = 'de_clusters.Rdata')
load(file = 'ROTS_summary_pop.Rdata')
head(summary_pop1)
head(summary_pop2)
head(summary_pop3)
head(summary_pop4)



# Extract DE gene names
de_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.05)
# sgene_names <- gene_names$genes


g1=rownames(summary_pop1[summary_pop1$FDR>0.05,])
g2=rownames(gene_names[gene_names$cluster=='C1',])
length(intersect(g1,g2));length(g1);length(g2)

g1=rownames(summary_pop2[summary_pop2$FDR>0.05,])
g2=rownames(gene_names[gene_names$cluster=='C2',])
length(intersect(g1,g2));length(g1);length(g2)

g1=rownames(summary_pop3[summary_pop3$FDR>0.05,])
g2=rownames(gene_names[gene_names$cluster=='C3',])
length(intersect(g1,g2));length(g1);length(g2)

g1=rownames(summary_pop4[summary_pop4$FDR>0.05,])
g2=rownames(gene_names[gene_names$cluster=='C4',])
length(intersect(g1,g2));length(g1);length(g2)


g1=head(rownames(summary_pop4[summary_pop4$FDR>0.05,]),20)
g2=head(rownames(gene_names[gene_names$cluster=='C4',]),20)
length(intersect(g1,g2));length(g1);length(g2)

## limma+voom, deseq2, edgeR



