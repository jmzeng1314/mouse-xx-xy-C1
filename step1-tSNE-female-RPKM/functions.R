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

# 加载必备的R包
library('data.table')
library('Seurat')
library('matrixStats')
library('statmod')
library('ggplot2')
library("Rtsne")
library("factoextra")
library('FactoMineR')

# 调色
female_stagePalette <- 	c(
  "#2754b5", 
  "#8a00b0", 
  "#d20e0f", 
  "#f77f05", 
  "#f9db21",
  "#43f14b"
)


###########################################
#										  #
#			Var Gene Selection		      #
#										  #
###########################################
# 这个函数摘抄自 Brennecke et al 2013 的文章：
# Extract genes with a squared coefficient of variation >2 times 
# the fit regression (Brennecke et al 2013 method)

# Code taken and adapted from doi:10.1038.nmeth.2645
# Select the highly variable genes based on the 
# squared coefficient of variation and the mean gene expression 
# and return the RPKM matrix the the HVG

getMostVarGenes <- function(
  data=data,				# RPKM matrix
  fitThr=1.5, 			# Threshold above the fit to select the HGV
  minMeanForFit=1			# Minimum mean gene expression level
){
 # data=females;fitThr=2;minMeanForFit=1	
  # Remove genes expressed in no cells
  data_no0 <- as.matrix(
    data[rowSums(data)>0,]
  )
  # Compute the mean expression of each genes
  meanGeneExp <- rowMeans(data_no0)
  names(meanGeneExp)<- rownames(data_no0)
  
  # Compute the squared coefficient of variation
  varGenes <- rowVars(data_no0)
  cv2 <- varGenes / meanGeneExp^2
  
  # Select the genes which the mean expression is above the expression threshold minMeanForFit
  useForFit <- meanGeneExp >= minMeanForFit
  
  # Compute the model of the CV2 as a function of the mean expression using GLMGAM
  fit <- glmgam.fit( cbind( a0 = 1, 
                            a1tilde = 1/meanGeneExp[useForFit] ), 
                     cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  
  # Get the highly variable gene counts and names
  fit_genes <- names(meanGeneExp[useForFit])
  cv2_fit_genes <- cv2[useForFit]
  fitModel <- fit$fitted.values
  names(fitModel) <- fit_genes
  HVGenes <- fitModel[cv2_fit_genes>fitModel*fitThr]
  print(length(HVGenes))
  
  # Plot the result
  plot_meanGeneExp <- log10(meanGeneExp)
  plot_cv2 <- log10(cv2)
  plotData <-  data.frame(
    x=plot_meanGeneExp[useForFit],
    y=plot_cv2[useForFit],
    fit=log10(fit$fitted.values),
    HVGenes=log10((fit$fitted.values*fitThr))
  )
  p <- ggplot(plotData, aes(x,y)) +
    geom_point(size=0.1) +
    geom_line(aes(y=fit), color="red") +
    geom_line(aes(y=HVGenes), color="blue") +
    theme_bw() +
    labs(x = "Mean expression (log10)", y="CV2 (log10)")+
    ggtitle(paste(length(HVGenes), " selected genes", sep="")) +
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=16),
      legend.text = element_text(size =16),
      legend.title = element_text(size =16 ,face="bold"),
      legend.position= "none",
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      aspect.ratio=1
    )+
    scale_color_manual(
      values=c("#595959","#5a9ca9")
    )
  print(p)
  
  # Return the RPKM matrix containing only the HVG
  HVG <- data_no0[rownames(data_no0) %in% names(HVGenes),]
  return(HVG)
}


###########################################
#                                         #
#              RtSNE Analysis             #
#                                         #
###########################################

# Compute the t-SNE
run_tSNE <- function(
  pca=pca, 					# PCA computed with FactoMineR package
  pc=20, 						# Number of PCs to use in the t-SNE, default is 20
  iter=2000, 					# Number of iteration for th rt-SNE, default is 2000
  perplexity=30				# Set the t-SNE perplexity parameter if needed, default is 30
){
  set.seed(1)
  # Run t-SNE
  rtsne_out <- Rtsne(
    pca$ind$coord[,pc], 
    pca=FALSE, 
    max_iter=iter, 
    verbose=TRUE
  )
  tSNE <- data.frame(
    rtsne_out$Y
  )
  return(tSNE)
}

# Compute and plot the t-SNE
run_plot_tSNE <- function(
  pca=pca, 
  pc=pc, 
  iter=2000,
  perplexity=30, 
  conditions=conditions, 
  colours=colours, 
  title=FALSE
){
  tsne <- run_tSNE(
    pca, 
    1:pc, 
    iter,  
    perplexity
  )
  tSNE <- data.frame(
    tsne,
    conditions
  )
  colnames(tSNE)<- c("tSNE_1", "tSNE_2", "cond")
  if (title==TRUE){
    g<-ggplot(tSNE, aes(tSNE_1,tSNE_2)) +
      geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
      theme_bw() +
      scale_fill_manual(
        values=colours,
        name=""
      ) +
      scale_color_manual(
        values=colours,
        name=""
      ) +
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      ggtitle(paste("t-SNE plot (", pc," PCs)", sep="")) +
      theme(
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_text(size =16 ,face="bold"),
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        aspect.ratio=1#,
        #legend.position="none"
      )
  } else {
    g<-ggplot(tSNE, aes(tSNE_1,tSNE_2)) +
      geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
      theme_bw() +
      scale_fill_manual(
        values=colours,
        name=""
      ) +
      scale_color_manual(
        values=colours,
        name=""
      ) +
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      theme(
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_text(size =16 ,face="bold"),
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        aspect.ratio=1#,
        #legend.position="none"
      )
  }
  print(g)
  return(tSNE)
}


# Plot the pre-computed t-SNE
plot_tSNE <- function(tsne=tsne, conditions=conditions, colours=colours){
  
  tSNE <- data.frame(
    tSNE_1=tsne[,1],
    tSNE_2=tsne[,2],
    c(conditions)
  )
  
  colnames(tSNE)<- c("tSNE_1", "tSNE_2", "cond")
  
  print(head(tSNE))
  
  
  g<-ggplot(tSNE, aes(tSNE_1, tSNE_2)) +
    geom_point(shape = 21, size = 2, stroke=0.5, aes(fill=cond, color=cond),alpha=6/10) +
    theme_bw() +
    scale_fill_manual(
      values=colours,
      name=""
    ) +
    scale_color_manual(
      values=colours,
      name=""
    ) +
    # ggtitle(paste("t-SNE plot (", pc," PCs)", sep="")) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=16),
      legend.text = element_text(size =16),
      legend.title = element_text(size =16 ,face="bold"),
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      aspect.ratio=1#,
      #legend.position="none"
    )
  print(g)
  return(tSNE)
}


# t-SNE colored by gene expression level
tsne_gene_per_cell <- function(tsne_result, rpkm){
  gene_per_cell <- rpkm
  gene_per_cell[gene_per_cell>0] <- 1
  gene_per_cell <- colSums(gene_per_cell)
  
  colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
  mypalette <- c(warm(20))
  
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  
  p <- ggplot(tsne_result, aes(tSNE_1, tSNE_2)) +
    geom_point(shape = 21, stroke=0.25, aes(colour=as.numeric(gene_per_cell), fill=as.numeric(gene_per_cell)), size = 2) +
    # scale_fill_gradient2(high="darkred", low="yellow")+
    # scale_fill_gradientn(colours = mypalette)+
    # scale_colour_gradientn(colours = mypalette)+
    scale_fill_viridis_c()+
    scale_colour_viridis_c()+
    theme_bw() +
    # ggtitle(gene) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    theme(
      plot.title = element_text(size=28, face="bold.italic", hjust = 0),
      # axis.text=element_text(size=12),
      # axis.title=element_text(size=16),
      axis.title=element_blank(),
      axis.text=element_blank(),
      legend.text = element_text(size =12),
      legend.title=element_blank(),
      aspect.ratio=1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  # print(p)
  
}



# t-SNE colored by gene expression level
tsne_gene_exp <- function(tsne_result, gene, rpkm){
  colnames(tsne_result)<- c("tSNE_1", "tSNE_2")
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
  mypalette <- c(warm(20))
  
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  gene_exp <- rpkm[gene,,drop=FALSE]
  gene_exp <- log(gene_exp+1)
  density_exp <- tsne_result[rownames(tsne_result) %in% names(gene_exp[,gene_exp>rowMeans(gene_exp),drop=FALSE]),]
  # print(density_exp)
  
  title <- grobTree(textGrob(gene, x=0.05, y=0.93, , hjust = 0, gp=gpar(fontsize=32, fontface="bold.italic")))
  
  p <- ggplot(tsne_result, aes(tSNE_1, tSNE_2)) +
    geom_point(shape = 21, stroke=0.25, aes(colour=as.numeric(log(rpkm[gene,]+1)), fill=as.numeric(log(rpkm[gene,]+1))), size = 2) +
    geom_density_2d(data=density_exp, aes(x=tSNE_1, y=tSNE_2), bins = 5, colour="black") +
    # scale_fill_gradient2(high="darkred", low="yellow")+
    scale_fill_gradientn(colours = mypalette)+
    scale_colour_gradientn(colours = mypalette)+
    # scale_fill_viridis_c()+
    theme_bw() +
    # ggtitle(gene) +
    annotation_custom(title) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    theme(
      plot.title = element_text(size=28, face="bold.italic", hjust = 0),
      # axis.text=element_text(size=12),
      # axis.title=element_text(size=16),
      axis.title=element_blank(),
      axis.text=element_blank(),
      legend.text = element_text(size =12),
      legend.title=element_blank(),
      aspect.ratio=1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none",
      axis.ticks = element_blank()
    )
  # print(p)
  
}

violin_gene_exp <- function(gene, rpkm, conditions=conditions, colours=colours, test=TRUE){
  exp <- as.numeric(log(rpkm[gene,]+1))
  
  gene_exp <- data.frame(
    cell=colnames(rpkm),
    clusters=conditions,
    gene=exp
  )
  
  if (test==TRUE){
    # Perform pairwise comparisons
    test <- compare_means(gene ~ clusters, data = gene_exp, method = "kruskal.test")
    # print(test)
    # print(compare_means(gene ~ clusters, data = gene_exp, , method = "wilcox.test", ref.group = ".all."))
    
    p <- ggboxplot(gene_exp, x = "clusters", y = "gene", color = "white")+
      geom_violin(scale = "width", width=0.7, adjust = .5,aes(fill=clusters)) +
      stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke = 1, fill="white")+
      
      # geom_jitter(size=0.3)+  	  
      geom_hline(yintercept = mean(gene_exp$gene), linetype = 2)+
      
      # geom_boxplot(fill="white", outlier.shape = NA, width = 0.2)+
      scale_fill_manual(
        values=colours
      ) +	  
      theme_bw()+
      ggtitle(gene)+
      expand_limits(y = c(0, max(gene_exp$gene)+1.5)) +
      # stat_compare_means(method = "kruskal.test", label.y = max(gene_exp$gene)+1)+      # Add global p-value
      stat_compare_means(
        label = "p.signif", 
        method = "wilcox.test", 
        ref.group = ".all.", 
        label.y = max(gene_exp$gene)+0.75, 
        size=6
      )+ # Pairwise comparison against all
      theme(
        plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
        axis.text=element_text(size=16),
        # axis.title=element_text(size=16),
        axis.title=element_blank(),
        legend.text = element_text(size =16),
        legend.title=element_blank(),
        aspect.ratio=0.5,
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
  }else{
    p <- ggplot(gene_exp, aes(clusters, gene))+
      geom_violin(scale = "width", width=0.7, adjust = .5,aes(fill=clusters)) +
      stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke = 1, fill="white")+
      # geom_jitter(size=0.3)+ 
      scale_fill_manual(
        values=colours
      ) +	  
      theme_bw()+
      ggtitle(gene)+
      theme(
        plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
        axis.text=element_text(size=16),
        # axis.title=element_text(size=16),
        axis.title=element_blank(),
        legend.text = element_text(size =16),
        legend.title=element_blank(),
        aspect.ratio=0.5,
        legend.position="none"
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
      )
  }
  
  
  print(p)
}


