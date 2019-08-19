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


# 使用monocle包的示例标准代码：
prepare_for_DE <- function(count_matrix=count_matrix, clustering=clustering, stages=stages){
  # Prepare tables for monocle object
  expr_matrix <- as.matrix(count_matrix)
  sample_sheet <- data.frame(cells=names(count_matrix), 
                             stages=stages, 
                             cellType=clustering)
  rownames(sample_sheet)<- names(count_matrix)
  gene_annotation <- as.data.frame(rownames(count_matrix))
  rownames(gene_annotation)<- rownames(count_matrix)
  colnames(gene_annotation)<- "genes"
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  
  # Create a CellDataSet from the relative expression levels
  HSMM <- newCellDataSet(
    as(expr_matrix, "sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit=0.5,
    expressionFamily=negbinomial.size()
  )
  
  HSMM <- detectGenes(HSMM, min_expr = 5)
  # HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
  HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  return(HSMM)
}

# Compute DE genes
findDEgenes <- function(HSMM=HSMM, qvalue=qvalue){
  diff_test_res <- differentialGeneTest(
    HSMM,
    fullModelFormulaStr="~cellType",
    cores = 3
  )
  
  sig_genes_0.05 <- subset(diff_test_res, qval < 0.05)
  sig_genes_0.01 <- subset(diff_test_res, qval < 0.01)
  
  print(paste(nrow(sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
  print(paste(nrow(sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))
  
  diff_test_res <- subset(diff_test_res, qval< qvalue)
  
  return(diff_test_res)
}

plot_DE_heatmap <- function(rpkm_matrix=rpkm_matrix, diff_test_res=diff_test_res, qvalue=qvalue, condition=condition) {
  sig_genes <- subset(diff_test_res, qval < qvalue)
  dim(sig_genes)
  de_genes <- as.matrix(log(rpkm_matrix[rownames(rpkm_matrix) %in% rownames(sig_genes),]+1))
  plot_heatmap(de_genes, condition)
}

get_up_reg_clusters <- function(count, clustering, DE_genes){
  cluster_nb <- unique(clustering)
  mean_per_cluster <- vector()
  DE_genes <- DE_genes[order(rownames(DE_genes)),]
  count <- count[order(rownames(count)),]
  count_de_genes <- count[rownames(count) %in% DE_genes$genes,]
  print(dim(count_de_genes))
  for (clusters in cluster_nb) {
    # print(head(count_de_genes[,
    # 		colnames(count_de_genes) %in% names(clustering[clustering==clusters])
    # 	]))
    mean <- rowMeans(
      as.matrix(count_de_genes[,
                               colnames(count_de_genes) %in% names(clustering[clustering==clusters])
                               ])
    )
    names(mean) <- clusters
    mean_per_cluster <- cbind(
      mean_per_cluster,
      mean
    )
  }
  colnames(mean_per_cluster) <- cluster_nb
  up_reg_cluster <- colnames(mean_per_cluster)[apply(mean_per_cluster,1,which.max)]
  de_genes_table <- data.frame(
    DE_genes,
    mean_per_cluster,
    cluster=up_reg_cluster
  )
  
  return(de_genes_table)
}

get_top_up_reg_clusters <- function(DE_genes, gene_nb){
  cluster_nb <- unique(DE_genes$cluster)
  cluster_nb <- cluster_nb[order(cluster_nb)]
  DE_genes <- DE_genes[order(DE_genes$qval),]
  DE_genes <- DE_genes[DE_genes$qval<1.10e-2,]
  
  top_DE_genes <- vector()
  for (clusters in cluster_nb) {
    genes <- rownames(DE_genes[DE_genes$cluster==clusters,])[1:gene_nb]
    top_DE_genes <- c(top_DE_genes, genes)
  }
  
  return(top_DE_genes)
}
