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



library(pheatmap)
library(viridisLite)
library(ggpubr)
library(grid)
library(Seurat)


# 这里是自定义的热图，比较复杂，修改了pheatmap的非常多参数。 
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


# 这里是自定义的tSNE等高线图，同样是看表达量。 
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

# 这里是自定义的小提琴图，看表达量 
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
 


