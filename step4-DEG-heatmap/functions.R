# Code taken and adapted from doi:10.1038.nmeth.2645
# Select the highly variable genes based on the squared coefficient of variation and the mean gene expression and return the RPKM matrix the the HVG
library(matrixStats)
library(statmod)
library(ggplot2)
library("Rtsne")
library("factoextra")
library(FactoMineR)
library(viridis)

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


plot_heatmap <- function(matrix=matrix, clusters=clusters, stages=stages){
  annotation_col <- data.frame(
    Cell_Clusters=clusters,
    Stages=stages
  )
  
  cell_cluster_colors <- c("#bd66d4", "#77b0f3", "#8dcf38", "#ffd65c", "#fb7072")
  names(cell_cluster_colors) <- unique(clusters)
  
  annotation_colors <- list(
    Stages=c(P6="#b0f137", E16.5="#fc8d59", E13.5="#fee090", E12.5="#e0f3f8", E11.5="#91bfdb", E10.5="#4575b4"),
    Cell_Clusters=cell_cluster_colors
  )
  
  pheatmap(
    matrix,  
    show_colnames=FALSE, 
    show_rownames=FALSE, 
    clustering_method="ward.D2",
    annotation_col=annotation_col,
    annotation_colors=annotation_colors,
    color=viridis(100)
  )
  
}


plot_heatmap_2 <- function(matrix=matrix, clusters=clusters, stages=stages, rowbreaks, colbreaks, cluster_color, stage_color){
  annotation_col <- data.frame(
    Cell_Clusters=clusters,
    Stages=stages
  )
  
  annotation_colors <- list(
    Stages=stage_color,
    Cell_Clusters=cluster_color
  )
  
  # Color palette for the heatmap
  cold <- colorRampPalette(c('#41b6c4','#253494','#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58', '#081d58'))
  warm <- colorRampPalette(c('#fecc5c','#e31a1c','#800026','#800026','#800026','#800026','#800026','#800026'))
  mypalette <- c(rev(cold(20)), warm(20))
  # mypalette <- c(rev(cold(15)), warm(16))
  breaksList = seq(0, 5, by = 0.5)
  
  
  pheatmap(
    matrix, 
    # scale="row",
    show_colnames=FALSE, 
    # show_rownames=FALSE, 
    cluster_cols=FALSE,
    # cluster_rows=FALSE,
    clustering_method="ward.D",
    annotation_col=annotation_col,
    annotation_colors=annotation_colors,
    color=viridis(10),
    # color=mypalette,
    # gaps_row=rowbreaks,
    gaps_col=colbreaks,
    border_color = FALSE
    # breaks=breaksList
  )
}



