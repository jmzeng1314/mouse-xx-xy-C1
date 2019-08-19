# Code taken and adapted from doi:10.1038.nmeth.2645
# Select the highly variable genes based on the squared coefficient of variation and the mean gene expression and return the RPKM matrix the the HVG
library(matrixStats)
library(statmod)
library(ggplot2)
library("Rtsne")
library("factoextra")
library(FactoMineR)
library(viridis)
library("scatterplot3d")
library(slingshot)

###########################################
#                                         #
#             Diffusion map               #
#                                         #
###########################################

#  DiffusionMap class {destiny}
run_diffMap <- function(data=data, condition=condition, sigma="local"){
  # data=females_data
  # condition=female_clustering
  # sigma=15
  destinyObj <- as.ExpressionSet(as.data.frame(t(data)))
  destinyObj$condition <- factor(condition)
  dm <- DiffusionMap(destinyObj, sigma, rotate = TRUE)
  return(dm)
}

plot_eigenVal <- function(dm=dm){
  linepad <- .5
  plot(
    eigenvalues(dm), 
    ylim = 0:1, 
    pch = 20, 
    xlab ='Diffusion component (DC)', 
    ylab ='Eigenvalue'
  )
}

plot_dm_2D <- function(dm=dm, dc=2, condition=condition, colours=colours){
  DCs <- paste("DC",1:dc, sep="")
  
  dm_eigen <- data.frame(
    dm@eigenvectors
  )
  
  DCs.combinations <- combn(DCs,2)
  g <- apply(
    DCs.combinations,
    2,
    function(combination)
    {
      p1 <- ggplot(dm_eigen, aes_string(x=combination[1], y=combination[2])) +
        geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=condition)) +
        theme_bw() +
        scale_fill_manual(
          values=colours,
          name=""
        ) +
        xlab(combination[1])+
        ylab(combination[2])+
        ggtitle("Diffusion Map")+
        theme(
          axis.text=element_text(size=16),
          axis.title=element_text(size=16),
          legend.text = element_text(size =16),
          legend.title = element_text(size =16 ,face="bold"),
          plot.title = element_text(size=18, face="bold", hjust = 0.5),
          aspect.ratio=1
        )
      print(p1)
    }
  )
  
}


plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
  cond <- factor(condition)
  col <- factor(condition)
  levels(col) <- colours
  col <- as.vector(col)
  DCs <- paste("DC",dc, sep="")
  
  data <- data.frame(
    dm@eigenvectors[,DCs[1]], 
    dm@eigenvectors[,DCs[2]], 
    dm@eigenvectors[,DCs[3]]
  )
  colnames(data) <- DCs
  
  plot3d(
    data,
    col=col,
    size=6.5,
    box = FALSE
  )
 
}




plot_dm_3D_transparent <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
  cond <- factor(condition)
  col <- factor(condition)
  levels(col) <- colours
  col <- as.vector(col)
  DCs <- paste("DC",dc, sep="")
  
  data <- data.frame(
    dm@eigenvectors[,DCs[1]], 
    dm@eigenvectors[,DCs[2]], 
    dm@eigenvectors[,DCs[3]]
  )
  colnames(data) <- DCs
  
  plot3d(
    data,
    col=col,
    size=6.5,
    alpha=0,
    box = FALSE
  )
  
  legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
  
}




###########################################
#                                         #
#          Lineage prediction             #
#                                         #
###########################################

get_lineage <- function(dm=dm, 
                        dim=c(0:20), 
                        condition=condition, 
                        start=start, end=end, 
                        shrink.method="cosine"){
  
  data <- data.frame(
    dm@eigenvectors[,dim]
  )
  crv <- slingshot(
    data, 
    condition, 
    start.clus = start, 
    end.clus=end,
    maxit=100000,
    shrink.method=shrink.method
    # shrink.method="tricube"
  )
  
  return(crv)
}


# Thanks to Guillaume Devailly for this function
rankKeepNA <- function(x) {
  return(
    ifelse(
      is.na(x),
      NA,
      rank(
        x, 
        na.last = TRUE, 
        ties.method="random"
      )
    )
  )
}


get_pseudotime <- function(pseudotime, wthres=wthres, ranked=TRUE){
  pseudoT <- list()
  for(lineage in 1:length(pseudotime@curves))local({
    curve <- pseudotime@curves[[lineage]]
    lambda <- curve$lambda
    weight <- curve$w
    ps <- curve$lambda
    ps[weight < wthres] <- NA
    if (ranked==TRUE){
      ps <- rankKeepNA(ps)
    }
    pseudoT[[lineage]] <<- ps
  })
  df <- t(do.call("rbind",pseudoT))
  colnames(df) <- names(pseudotime@curves)
  return(df)
}



plot_gene_per_lineage <- function(
  rpkm_matrix=rpkm_matrix, 
  pseudotime=pseudotime, 
  geneName=geneName, 
  stages=stages, 
  clusters=clusters, 
  stage_colors=stage_colors,
  cluster_colors=cluster_colors
){
  
  myplots <- list()
  total_pseudotime <- vector()
  for (lineages in 1:ncol(pseudotime)){
    lineage <- as.vector(pseudotime[,lineages])
    total_pseudotime <- c(total_pseudotime, lineage)
    total_pseudotime <- na.omit(total_pseudotime)
  }
  max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
  max_pseudotime <- max(total_pseudotime)
  
  pseudotime_data <- data.frame(
    pseudotime=numeric(),
    lineage=numeric(),
    stages=character(),
    clusters=character(),
    gene=numeric()
  )
  
  for (lineages in 1:ncol(pseudotime)){
    
    i <- lineages
    
    lineage <- pseudotime[,lineages]
    sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)
    
    data <- data.frame(
      pseudotime=lineage,
      lineage=paste("Lineage ",lineages, sep=""),
      stages=stages,
      clusters=clusters,
      gene=t(sub_data)
    )
    colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")
    
    pseudotime_data <- rbind(pseudotime_data, data)
    
  }
  
  
  #
  # Old plot I keep for record in case I need it
  #
  
  # p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
  # scale_shape_manual(values = 21:25) +
  # geom_point(size = 3, stroke=0.7, aes(shape=condition1, fill=condition2), color="white", na.rm = TRUE)+
  # geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.75)+
  # ylab("log(RPKM+1)") +
  # theme_bw() +
  # ggtitle(geneName)+
  # # scale_color_manual(
  # # 	values=colours
  # # ) +
  # scale_fill_manual(
  # 	values=colours
  # ) +
  # guides(
  # 	fill = guide_legend(override.aes = list(color=colours, fill=colours))
  #    ) +
  # theme(
  # 	axis.text=element_text(size=16),
  # 	axis.title=element_text(size=16),
  # 	legend.text = element_text(size =16),
  # 	legend.title=element_blank(),
  # 	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
  # 	aspect.ratio=0.5#,
  # 	# legend.position="bottom"
  # ) +
  # facet_grid(.~lineage) +
  # theme(
  # 	strip.text.x = element_text(size = 16)
  # )
  
  # pseudotime_data[pseudotime_data=="Lineage 1"] <- "Sertoli lineage"
  # pseudotime_data[pseudotime_data=="Lineage 2"] <- "Leydig lineage"
  # pseudotime_data[pseudotime_data=="Lineage 3"] <- "Progenitor lineage"
  
  # p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
  # # geom_point(shape=21, size = 3, stroke=0.7, aes(fill=clusters), color="white", na.rm = TRUE)+
  # # geom_point(shape=108, size = 8,  aes(y=-1, color=stages), na.rm = TRUE)+
  # # geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.75)+
  # geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.5)+
  # ylab("log(RPKM+1)") +
  # theme_bw() +
  # ggtitle(geneName)+
  # scale_color_manual(
  # 	values=cluster_colors
  # ) +
  # scale_fill_manual(
  # 	values=cluster_colors
  # ) +
  # # ylim(0, max(pseudotime_data$gene))+
  # # guides(
  # # 	fill = guide_legend(override.aes = list(color=colours, fill=colours))
  # #    ) +
  # theme(
  # 	axis.text=element_text(size=16),
  # 	axis.title=element_text(size=16),
  # 	legend.text = element_text(size =16),
  # 	legend.title=element_blank(),
  # 	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
  # 	aspect.ratio=0.5,
  # 	legend.position="bottom"
  # ) +
  # # facet_grid(.~lineage) +
  # theme(
  # 	strip.text.x = element_text(size = 16)
  # )
  
  p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
    geom_point(shape=21, size = 2.5, stroke=0.5, aes(fill=clusters), color="white", na.rm = TRUE)+
    geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE)+
    geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.5)+
    ylab("log(RPKM+1)") +
    theme_bw() +
    ggtitle(geneName)+
    scale_color_manual(
      values=stage_colors
    ) +
    scale_fill_manual(
      values=cluster_colors
    ) +
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=16),
      legend.text = element_text(size =16),
      legend.title=element_blank(),
      plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
      aspect.ratio=0.8,
      legend.position="bottom"
    ) +
    facet_grid(.~lineage) +
    theme(
      strip.text.x = element_text(size = 16)
    )
  print(p)
}


plot_smoothed_gene_per_lineage <- function(
  rpkm_matrix=rpkm_matrix, 
  pseudotime=pseudotime, 
  lin=lin,
  geneName=geneName, 
  stages=stages, 
  clusters=clusters, 
  stage_colors=stage_colors,
  cluster_colors=cluster_colors,
  lineage_colors=lineage_colors
){
  max_time <- max(pseudotime[,c(1,2)], na.rm = TRUE)
  pseudotime <- pseudotime[,lin]
  
  lineage_colors <- lineage_colors[lin]
  
  myplots <- list()
  total_pseudotime <- vector()
  
  if (length(lin)==1){
    
    lineage <- pseudotime
    total_pseudotime <- c(total_pseudotime, lineage)
    total_pseudotime <- na.omit(total_pseudotime)
    max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
    max_pseudotime <- max(total_pseudotime)
    
    pseudotime_data <- data.frame(
      pseudotime=numeric(),
      lineage=numeric(),
      stages=character(),
      clusters=character(),
      gene=numeric()
    )
    
    
    lineage <- pseudotime
    sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)
    
    data <- data.frame(
      pseudotime=lineage,
      lineage=paste("Lineage ",lin, sep=""),
      stages=stages,
      clusters=clusters,
      gene=t(sub_data)
    )
    
    colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")
    
    # print(lin)
    
    pseudotime_data <- rbind(pseudotime_data, data)
    
    p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
      geom_point(shape=21, size = 3,  aes(fill=clusters), color="white", na.rm = TRUE)+
      geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE)+
      geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.5)+
      ylab("log(RPKM+1)") +
      theme_bw() +
      ggtitle(geneName)+
      # scale_color_manual(
      # 	values=cluster_colors
      # ) +
      scale_fill_manual(
        values=cluster_colors
      ) +
      scale_color_manual(
        values=stage_colors
      ) +
      expand_limits(x = c(0,max_time))+
      expand_limits(y = c(0,2))+
      theme(
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title=element_blank(),
        plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
        aspect.ratio=0.5,
        legend.position="bottom",
        strip.text.x = element_text(size = 16)
      )
    
    
  } else {
    
    for (lineages in 1:ncol(pseudotime)){
      lineage <- as.vector(pseudotime[,lineages])
      total_pseudotime <- c(total_pseudotime, lineage)
      total_pseudotime <- na.omit(total_pseudotime)
    }
    
    max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
    max_pseudotime <- max(total_pseudotime)
    
    pseudotime_data <- data.frame(
      pseudotime=numeric(),
      lineage=numeric(),
      stages=character(),
      clusters=character(),
      gene=numeric()
    )
    
    
    for (lineages in 1:ncol(as.data.frame(pseudotime))){
      lineage <- pseudotime[,lineages]
      sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)
      
      data <- data.frame(
        pseudotime=lineage,
        lineage=paste("Lineage ",lineages, sep=""),
        stages=stages,
        clusters=clusters,
        gene=t(sub_data)
      )
      
      colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")
      
      pseudotime_data <- rbind(pseudotime_data, data)
    }
    
    p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
      geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.5)+
      ylab("log(RPKM+1)") +
      theme_bw() +
      ggtitle(geneName)+
      scale_color_manual(
        values=lineage_colors
      ) +
      scale_fill_manual(
        values=lineage_colors
      ) +
      expand_limits(y = c(0,2))+
      theme(
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title=element_blank(),
        plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
        aspect.ratio=0.5,
        legend.position="bottom",
        strip.text.x = element_text(size = 16)
      )
    
  }
  
  print(p)
}


plot_smoothed_gene_per_lineage_sex <- function(
  rpkm_matrix=rpkm_matrix, 
  pseudotime=pseudotime, 
  lin=lin,
  geneName=geneName, 
  stages=stages, 
  clusters=clusters, 
  stage_colors=stage_colors,
  cluster_colors=cluster_colors,
  lineage_colors=lineage_colors
){
  max_time <- max(pseudotime[,c(1,2)], na.rm = TRUE)
  pseudotime <- pseudotime[,lin]
  
  lineage_colors <- lineage_colors[lin]
  
  myplots <- list()
  total_pseudotime <- vector()
  
  if (length(lin)==1){
    
    # lineage <- pseudotime
    # total_pseudotime <- c(total_pseudotime, lineage)
    # total_pseudotime <- na.omit(total_pseudotime)
    # max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
    # max_pseudotime <- max(total_pseudotime)
    
    lineage <- as.vector(pseudotime)
    total_pseudotime <- c(total_pseudotime, lineage)
    total_pseudotime <- na.omit(total_pseudotime)
    
    max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
    max_pseudotime <- max(total_pseudotime)
    
    pseudotime_data <- data.frame(
      pseudotime=numeric(),
      lineage=numeric(),
      stages=character(),
      clusters=character(),
      gene=numeric(),
      sex= character()
    )
    
    
    lineage <- pseudotime
    sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)
    data <- data.frame(
      pseudotime=lineage,
      lineage=paste("Lineage ",lin, sep=""),
      stages=stages,
      clusters=clusters,
      gene=t(sub_data),
      sex=sapply(strsplit(clusters, "_"), `[`, 1)
    )
    
    colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene", "sex")
    pseudotime_data <- rbind(pseudotime_data, data)
    
    print(head(pseudotime_data))
    
    pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 1" & pseudotime_data[,"sex"]=="XX"),]
    pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 2" & pseudotime_data[,"sex"]=="XY"),]
    
    rect1 <- data.frame(xmin=240, xmax=290, ymin=-Inf, ymax=Inf)
    rect2 <- data.frame(xmin=380, xmax=410, ymin=-Inf, ymax=Inf)
    
    p3 <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
      # geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
      geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
      geom_point(shape=21, size = 2,  aes(fill=clusters), color="white", stroke = 0.3, na.rm = TRUE, position="dodge", alpha=0.8)+
      geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE, position="dodge")+
      geom_smooth(aes(group=sex, color=sex), size=2, se=FALSE, na.rm = TRUE, method="loess", span=0.5)+
      ylab("log(RPKM+1)") +
      theme_bw()+
      ggtitle(geneName)+
      # scale_color_manual(
      # 	values=cluster_colors
      # ) +
      scale_fill_manual(
        values=cluster_colors
      ) +
      scale_color_manual(
        values=stage_colors
      ) +
      expand_limits(x = c(0,max_time))+
      expand_limits(y = c(0,2))+
      theme(
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size =12),
        legend.title=element_blank(),
        plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
        aspect.ratio=0.5,
        legend.position="none",
        strip.text.x = element_text(size = 12)
      )
    
    
  } else {
    
    for (lineages in 1:ncol(pseudotime)){
      lineage <- as.vector(pseudotime[,lineages])
      total_pseudotime <- c(total_pseudotime, lineage)
      total_pseudotime <- na.omit(total_pseudotime)
    }
    
    max_exp <- max(log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1))
    max_pseudotime <- max(total_pseudotime)
    
    pseudotime_data <- data.frame(
      pseudotime=numeric(),
      lineage=numeric(),
      stages=character(),
      clusters=character(),
      gene=numeric(),
      sex= character()
    )
    
    
    for (lineages in 1:ncol(as.data.frame(pseudotime))){
      lineage <- pseudotime[,lineages]
      sub_data <- log(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]+1)
      
      data <- data.frame(
        pseudotime=lineage,
        lineage=paste("Lineage ",lineages, sep=""),
        stages=stages,
        clusters=clusters,
        gene=t(sub_data),
        sex=sapply(strsplit(clusters, "_"), `[`, 1)
      )
      
      colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene", "sex")
      
      pseudotime_data <- rbind(pseudotime_data, data)
    }
    
    
    pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 1" & pseudotime_data[,"sex"]=="XX"),]
    pseudotime_data<-pseudotime_data[!(pseudotime_data[,"lineage"]=="Lineage 2" & pseudotime_data[,"sex"]=="XY"),]
    
    rect1 <- data.frame(xmin=240, xmax=290, ymin=-Inf, ymax=Inf)
    rect2 <- data.frame(xmin=380, xmax=410, ymin=-Inf, ymax=Inf)
    
    p3 <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
      geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
      geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
      geom_point(shape=21, size = 2,  aes(fill=clusters), color="white", stroke = 0.3, na.rm = TRUE, position="dodge", alpha=0.8)+
      geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE, position="dodge")+
      geom_smooth(aes(group=lineage, color=lineage), size=2, se=FALSE, na.rm = TRUE, method="loess", span=0.5)+
      ylab("log(RPKM+1)") +
      theme_bw()+
      ggtitle(geneName)+
      # scale_color_manual(
      # 	values=cluster_colors
      # ) +
      scale_fill_manual(
        values=cluster_colors
      ) +
      scale_color_manual(
        values=stage_colors
      ) +
      expand_limits(x = c(0,max_time))+
      expand_limits(y = c(0,2))+
      theme(
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size =12),
        legend.title=element_blank(),
        plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
        aspect.ratio=0.5,
        legend.position="none",
        strip.text.x = element_text(size = 12)
      )
    
  }
  # print(p1)
  # print(p3)
  return(p3)
}




plot_aligned_lineage_sex <- function(
  rpkm_matrix=rpkm_matrix, 
  pseudotime=pseudotime, 
  lin=lin,
  stage_colors=stage_colors
){
  
  sex <- sapply(strsplit(rownames(pseudotime), "_"), `[`, 2)
  
  lin_female <- pseudotime[grep("XX",rownames(pseudotime)), lin]
  lin_male <- pseudotime[grep("XY",rownames(pseudotime)), lin]
  lin_pseudotime <- c(lin_female, lin_male)
  
  cells <- c(names(lin_female), names(lin_male))
  sex <- sapply(strsplit(cells, "_"), `[`, 2)
  stages <- sapply(strsplit(cells, "_"), `[`, 1)
  
  
  pseudotime_data <- data.frame(
    pseudotime=lin_pseudotime,
    stages=stages,
    sex= sex
  )
  
  rownames(pseudotime_data) <- cells
  
  print(head(pseudotime_data))
  
  # rect1 <- data.frame(xmin=240, xmax=290, ymin=-Inf, ymax=Inf)
  # rect2 <- data.frame(xmin=380, xmax=410, ymin=-Inf, ymax=Inf)
  
  g <- ggplot(pseudotime_data, aes(x=pseudotime, y=sex))+
    # geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
    geom_point(size = 4,  aes(color=stages), stroke = 0.3, na.rm = TRUE, position="dodge", alpha=0.8)+
    # geom_point(shape=108, size = 4,  aes(y=-1, color=stages), na.rm = TRUE, position="dodge")+
    # geom_smooth(aes(group=sex, color=sex), size=2, se=FALSE, na.rm = TRUE, method="loess", span=0.5)+
    ylab("pseudotime") +
    theme_bw()+
    
    scale_color_manual(
      values=stage_colors
    ) +
    theme(
      axis.text=element_text(size=12),
      axis.title=element_text(size=12),
      legend.text = element_text(size =12),
      legend.title=element_blank(),
      plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
      aspect.ratio=0.05,
      legend.position="none",
      strip.text.x = element_text(size = 12)
    )
  
  print(g)
}




compare_lineage <- function(pseudotime=pseudotime, condition=condition){
  myLineages <- list()
  for (lineages in 1:length(pseudotime))local({
    i <- lineages
    lineage <- pseudotime[[lineages]]$pseudotime
    myLineages[[i]] <<- lineage
  })
  df <- t(do.call("rbind",myLineages))
  df[df>=0] <- 1
  df[is.na(df)] <- 0
  sdf <- as.character(rowSums(df))
  return(sdf)
  
}
