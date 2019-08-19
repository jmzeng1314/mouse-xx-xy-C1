rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(Seurat))
library(ggfortify)

load('../females_rpkm.Rdata')
count_matrix=females
count_matrix[1:4,1:4]
dim(count_matrix)
clu2=read.csv('female_clustering.csv')
cluster=clu2[,2]
table(cluster) 
table(cluster)
col=rainbow(4)[as.factor(cluster)]
table(col)
## only tSNE
if(T){
  choose_T_counts=count_matrix
  choose_T_counts=choose_T_counts[apply(choose_T_counts,1, sd)>0,]
  choose_T_counts=choose_T_counts[names(tail(sort(apply(choose_T_counts,1, sd)),1000)),]
  pca_dat <- prcomp(t(choose_T_counts), scale. = TRUE)
  p=autoplot(pca_dat,col=col) + theme_classic() + ggtitle('PCA plot')
  print(p)
  
  str(pca_dat)
  pca_dat$x[1:4,1:4]
  new_pca <- FactoMineR::PCA(
    t(choose_T_counts), 
    ncp = ncol(choose_T_counts), 
    graph=FALSE
  )
  set.seed(42)
  # TSNE即t-distributed Stochastic Neighbor Embedding.
  library(Rtsne)
  tsne_out <- Rtsne(pca_dat$x[,1:5], perplexity = 10,
                    pca=FALSE, 
                    max_iter=2000, 
                    verbose=TRUE) # Run TSNE
  tsnes=tsne_out$Y
  colnames(tsnes) <- c("tSNE1", "tSNE2")
  ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(col=col)
  # 拓展，使用 kmeans或者DBSCAN 进行分群
  # 这里仍然使用 Hierarchical Clustering On Principle Components (HCPC)
  res.hcpc <- HCPC(
    new_pca, 
    graph = FALSE,
    min=4
  )
  new_c=res.hcpc$data.clust$clust
  table(new_c,cluster)
  
}