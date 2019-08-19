# Code taken and adapted from doi:10.1038.nmeth.2645
# Select the highly variable genes based on the squared coefficient of variation and the mean gene expression and return the RPKM matrix the the HVG
library(matrixStats)
library(statmod)
library(ggplot2)
library("Rtsne")
library("factoextra")
library(FactoMineR)
library(viridis)


###########################################
#                                         #
#      一些 monocle 相关的函数            #
#                                         #
###########################################
library(monocle)
# 使用 monocle 根据 read counts (size-factor) 来构建对象
prepare_for_monocle <- function(count_matrix=count_matrix, stages=stages){
  # Prepare tables for monocle object
  expr_matrix <- as.matrix(count_matrix)
  sample_sheet <- data.frame(cells=names(count_matrix), stages=stages)
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

# 同样是构建 monocle 对象，但是增加了两个属性，可供差异分析
prepare_for_DE <- function(count_matrix=count_matrix, clustering=clustering, stages=stages){
  # Prepare tables for monocle object
  expr_matrix <- as.matrix(count_matrix)
  sample_sheet <- data.frame(cells=names(count_matrix), stages=stages, cellType=clustering)
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

# 使用 monocle 根据 cellType 来找差异分析
# 所以之前构建的monocle对象里面需要有cellType参数
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



get_var_genes_pseudotime <- function(rpkm, count, pseudotime, 
                                     lineageNb=0, clusters){
  if(F){
    rpkm=females
    count=female_count
    pseudotime=female_pseudotime
    lineageNb=1
    clusters=female_clustering
  }
  
  lineage_pseudotime_all <- pseudotime[,lineageNb]
  lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]
  
  rpkm_exp_lineage <- rpkm[,rownames(lineage)]
  count_exp_lineage <- count[,rownames(lineage)]
  clusters_in_lineage <- clusters[names(clusters) %in% rownames(lineage)]
  table(clusters_in_lineage)
  
  print("Prepare for DE...")
  genes_pseudoT <- prepare_for_DE (
    count, 
    lineage_pseudotime_all, 
    lineage_pseudotime_all
  )
  print("Done.")
  
  genes_pseudoT$Pseudotime <- lineage_pseudotime_all
  
  print("Compute DE genes...")
  
  #Select cells from the lineage of interest
  genes_pseudoT <- genes_pseudoT[,rownames(lineage)]
  # Re-detect how many cells express each genes in the subset of cells
  genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)
  # Remove genes expressed in less than 10 cells
  genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 10, ]
  #  monocle这个包提供 differentialGeneTest 函数来做差异分析
  # 作用就是挑选那些在某些类别细胞里面高表达的基因，假设其为那一组细胞的marker基因。
  DE_genes_pseudoT <- differentialGeneTest(
    genes_pseudoT, 
    fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
    cores = 3
  )
  print("Done.")
  
  return(DE_genes_pseudoT)
}

###########################################
#                                         #
#            Loess smoothing              #
#                                         #
###########################################

smooth_gene_exp <- function(data=data, pseudotime=pseudotime, span=0.75){
  smooth_data <- data
  for (gene in 1:nrow(data)){
    gene_exp <- t(data[gene,])
    smooth <- loess(formula=gene_exp~pseudotime, span=span)
    smooth_data[gene,] <- predict(smooth, newdata=pseudotime)
  }
  return(smooth_data)
}


get_var_genes_pseudotime <- function(rpkm, count, pseudotime, lineageNb=0, clusters){
  
  
  lineage_pseudotime_all <- pseudotime[,lineageNb]
  lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]
  
  rpkm_exp_lineage <- rpkm[,rownames(lineage)]
  count_exp_lineage <- count[,rownames(lineage)]
  clusters_in_lineage <- clusters[names(clusters) %in% rownames(lineage)]
  
  print("Prepare for DE...")
  genes_pseudoT <- prepare_for_DE (
    count, 
    lineage_pseudotime_all, 
    lineage_pseudotime_all
  )
  print("Done.")
  
  genes_pseudoT$Pseudotime <- lineage_pseudotime_all
  
  print("Compute DE genes...")
  
  #Select cells from the lineage of interest
  genes_pseudoT <- genes_pseudoT[,rownames(lineage)]
  # Re-detect how many cells express each genes in the subset of cells
  genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)
  # Remove genes expressed in less than 10 cells
  genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 10, ]
  
  DE_genes_pseudoT <- differentialGeneTest(
    genes_pseudoT, 
    fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
    cores = 3
  )
  print("Done.")
  
  return(DE_genes_pseudoT)
}


get_var_genes_pseudotime_special <- function(rpkm, count, pseudotime, lineageNb=0, clusters){
  
  if (lineageNb > 0) {
    lineage_pseudotime_all <- pseudotime[,lineageNb]
    lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]
  } else {
    lineage_pseudotime_all <- pseudotime
    lineage <- pseudotime[!is.na(pseudotime),drop=FALSE]
  }
  
  
  rpkm_exp_lineage <- rpkm[,names(lineage)]
  count_exp_lineage <- count[,names(lineage)]
  clusters_in_lineage <- clusters[names(clusters) %in% names(lineage)]
  
  print("Prepare for DE...")
  genes_pseudoT <- prepare_for_DE (
    count, 
    lineage_pseudotime_all, 
    lineage_pseudotime_all
  )
  print("Done.")
  
  genes_pseudoT$Pseudotime <- lineage_pseudotime_all
  
  print("Compute DE genes...")
  
  #Select cells from the lineage of interest
  genes_pseudoT <- genes_pseudoT[,names(lineage)]
  # Re-detect how many cells express each genes in the subset of cells
  genes_pseudoT <- detectGenes(genes_pseudoT, min_expr = 5)
  # Remove genes expressed in less than 10 cells
  genes_pseudoT <- genes_pseudoT[fData(genes_pseudoT)$num_cells_expressed >= 10, ]
  
  DE_genes_pseudoT <- differentialGeneTest(
    genes_pseudoT, 
    fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
    cores = 3
  )
  print("Done.")
  
  return(DE_genes_pseudoT)
}


get_gene_clustering <- function(DE_genes_pseudoT, rpkm, pseudotime, lineageNb, cell_clusters, clusterNb=10, qvalue=0.1){
  lineage <- pseudotime[!is.na(pseudotime[,lineageNb]),lineageNb,drop=FALSE]
  clusters_in_lineage <- cell_clusters[names(cell_clusters) %in% rownames(lineage)]
  
  print(paste("Select genes with q-value<", qvalue, "...", sep=""))
  sig_gene_pseudoT <- row.names(subset(DE_genes_pseudoT, qval < qvalue))
  sig_var_gene_pseudoT <- rpkm[sig_gene_pseudoT,rownames(lineage)]
  
  print(paste(length(sig_gene_pseudoT), " selected genes.", sep=""))
  
  print("Smooth gene expression...")
  smooth <- smooth_gene_exp(
    log(sig_var_gene_pseudoT+1), 
    lineage[,1], 
    span=0.5
  )
  print("Done.")
  
  smooth <- smooth[ , order(lineage[,1])]
  
  print("Cluster the genes...")
  
  
  # cell_cluster_palette <- c(
  # 	C1="#560047", 
  # 	C2="#a53bad", 
  # 	C3="#eb6bac", 
  # 	C4="#ffa8a0"
  # )
  
  cell_cluster_palette <- c(
    C2="#aeff01", 
    C4="#009900", 
    C1="#457cff", 
    C3="#00c5ec",
    C6="#8dfaff",
    C5="#a38cff"
  )
  
  
  gene_clustering <- pheatmap(
    smooth, 
    scale="row", 
    cutree_rows=clusterNb,
    clustering_method="ward.D",
    silent=TRUE
  )
  
  print("Plot heatmap...")
  
  clusters <- cutree(gene_clustering$tree_row, k = clusterNb)
  clustering <- data.frame(clusters)
  clustering[,1] <- as.character(clustering[,1])
  colnames(clustering) <- "Gene_Clusters"
  
  gene_cluster_palette <- c(
    '#a6cee3',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
    '#fb9a99',
    '#e31a1c',
    '#fdbf6f',
    '#ff7f00',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928',
    '#49beaa',
    '#611c35',
    '#2708a0'
  )
  gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
  names(gene_cluster_colors) <- 1:max(clusters)
  
  
  annotation_row <- clustering
  
  annotation_col <- data.frame(
    Cell_Clusters=clusters_in_lineage,
    Stages=sapply(strsplit(names(clusters_in_lineage), "_"), `[`, 1)
  )
  
  annotation_colors <- list(
    Stages=c(
      E10.5="#2754b5", 
      E11.5="#8a00b0", 
      E12.5="#d20e0f", 
      E13.5="#f77f05", 
      E16.5="#f9db21",
      P6="#b0f137"
    ),
    Cell_Clusters=cell_cluster_palette,
    Gene_Clusters=gene_cluster_colors
  )
  
  # Color palette for the heatmap
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58', '#081d58'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
  mypalette <- c(rev(cold(20)), warm(21))
  # mypalette <- c(rev(cold(15)), warm(16))
  breaksList = seq(-10, 10, by = 0.5)
  
  pheatmap(
    smooth, 
    scale="row", 
    cluster_cols=FALSE, 
    show_colnames=FALSE, 
    show_rownames=FALSE, 
    cutree_rows=clusterNb, 
    clustering_method="ward.D",
    annotation_col=annotation_col,
    annotation_row=annotation_row,
    annotation_colors=annotation_colors,
    annotation_names_row=FALSE,
    color=mypalette
  )
  
  return(clustering)
}


