######################################################################################################################
#' MLG clustering 
#' This function takes a list of low dimensional embedding data as input, and performs MLG clustering.
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph
#' @param prune.param The prune parameter for SNN graph. There is an edge between cell i and j in the SNN graph, if the number of common neighbors between i and j exceeds the product of knn.param and knn.param.
#' @param cluster.resolution The resolution number of modularity maximization.
#' @param cluster.algorithm The clustering algorithm. 1--Louvain algorithm; 2--Louvain algorithm with multilevel refinement; 3--SLM algorithm.
#' @return A vector of cluster labels.
#' @export 
mlg_cluster <- function(factor.list, knn.param = 20, prune.param = 1/5, cluster.resolution, cluster.algorithm = 1){
  snn.list = .generate_snn_graph(factor.list, knn.param, prune.param)
  mlg = .generate_mlg(snn.list)
  cluster = Seurat::FindClusters(mlg, resolution= cluster.resolution, algorithm =cluster.algorithm)[[1]]
  return(cluster)
}

######################################################################################################################
######################################################################################################################
#' This function takes a list of low dimensional embedding data as input, construct a SNN graph for each low dimensional embeddings, and then plot a heatmap of the proportion of overlapping edges between each pair of SNN graphs.
#'
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph.
#' @param prune.param The prune parameter for SNN graph. There is an edge between cell i and j in the SNN graph, if the number of common neighbors between i and j exceeds the product of knn.param and knn.param.
#' @return A heatmap of the proportion of overlapping edges between each pair of SNN graphs constructed from low dimension embeddings.
#' @export 
prop.overlap.edges <- function(factor.list, knn.param=20, prune.param=1/5){
  snn.list = .generate_snn_graph(factor.list, knn.param, prune.param)
  n = length(snn.list)
  mat = matrix(NA, n, n)
  for (i in 1:n){
    for (j in i:n){
      mat[i,j]=sum(snn.list[[i]]*snn.list[[j]])*2/(sum(snn.list[[i]])+sum(snn.list[[j]]))
    }
  }
  rownames(mat)= names(snn.list)
  colnames(mat)= names(snn.list)
  melt_matrix =  reshape2::melt(mat)
  melt_matrix=melt_matrix[!is.na(melt_matrix$value),]
  
  p=ggplot2::ggplot(data = melt_matrix, ggplot2::aes(x=Var1, y=Var2, fill=value)) +
    ggplot2::geom_tile()+ggplot2::scale_fill_continuous(low="thistle2", high="darkred", 
                                                                                guide="colorbar",na.value="white",limits=c(0,1))+
    ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::labs(fill="Proportion of \noverlapping edges    ")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold",
                                     size=15),
          axis.text.y = ggplot2::element_text(face="bold",
                                     size=15),
          legend.text=ggplot2::element_text(size=10),
          legend.title=ggplot2::element_text(size=10, face="bold"),
          legend.position = "bottom",
          panel.background = ggplot2::element_blank())+ 
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size=5) +
    guides(fill = guide_colourbar(barwidth = 7, barheight = .5))
  return(p)
  
}

######################################################################################################################
######################################################################################################################
#' This function takes a list of low dimensional embedding data, and true cell type labels as input, and 
#' compute the graph signal to noise ratio of the SNN graphs constructed with each low dimension embeddding
#' and the MLG graph.
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph.
#' @param prune.param The prune parameter for SNN graph. There is an edge between cell i and j in the SNN graph, if the number of common neighbors between i and j exceeds the product of knn.param and knn.param.
#' @param cell_label True cell type label.
#' @return A vector of signal to noise ratio, labeled by graph name.
#' @export 
graph_signal_noise_ratio<- function(factor.list, knn.param=20, prune.param=1/5, cell_label){
  cell_label = as.integer(as.factor(cell_label))
  n_cells = length(cell_label)
  n_levels=length(unique(cell_label))
  snn.list = .generate_snn_graph(factor.list, knn.param, prune.param)
  M = matrix(0, nrow=n_levels, ncol = n_levels)
  n_layers = length(snn.list)
  SNR = rep(0, n_layers)
  for (l in 1:n_layers){
    for (i in 1:n_levels){
      for (j in i:n_levels){
        M[i,j] = Matrix::mean(snn.list[[l]][seq(n_cells)[cell_label == i], seq(n_cells)[cell_label == j]])
        M[j,i] = M[i,j]
      }
    }
    a=min(diag(M))
    b=max(M-diag(diag(M)))
    SNR[l]=(a-b)^2/a
  }
  if (n_layers==1){
    names(SNR)=names(snn.list)[1]
  }else{
    mlg= .generate_mlg(snn.list)
    for (i in 1:n_levels){
      for (j in i:n_levels){
        M[i,j] = Matrix::mean(mlg[seq(n_cells)[cell_label == i], seq(n_cells)[cell_label == j]])
        M[j,i] = M[i,j]
      }
    }
    a=min(diag(M))
    b=max(M-diag(diag(M)))
    SNR=c(SNR,(a-b)^2/a)
    names(SNR)=c(names(snn.list),"MLG")
  }
  return(SNR)
}

######################################################################################################################
######################################################################################################################
#' MLG  visualization
#' This function takes a list of low dimensional embedding data as input, and computes coordinates of the force-directed layout for MLG.
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph
#' @param prune.param The prune parameter for SNN graph. There is an edge between cell i and j in the SNN graph, if the number of common neighbors between i and j exceeds the product of knn.param and knn.param.
#' @param label Labels to be superimposed on the force directed layout of MLG.
#' @return A ggplot object.
#' @export 
mlg_visualization <- function(factor.list, knn.param = 20, prune.param = 1/5, label, label_title = "label"){
  snn.list = .generate_snn_graph(factor.list, knn.param, prune.param)
  mlg = .generate_mlg(snn.list)
  set.seed(1)
  net=igraph::graph_from_adjacency_matrix(mlg, mode =  "undirected",
                                  diag = F, add.colnames = NULL, add.rownames = NA)
  layer_layout=igraph::layout.fruchterman.reingold(net)
  layer_layout= apply(layer_layout,2, function(i){(i-mean(i))/sd(layer_layout[,2])})
  coord_mlg= data.frame(Coordinate_1 = layer_layout[,1], 
                        Coordinate_2 = layer_layout[,2])
  p<-ggplot(data=data.frame(coord_mlg), 
            ggplot2::aes(x=Coordinate_1, y=Coordinate_2, color=label)) +
            ggplot2::geom_point()+
            ggplot2::labs(color=label_title)+
            ggplot2::theme_bw()+
            ggplot2::theme(
          strip.text.x = ggplot2::element_text(size = 15, face="bold"),
          axis.title = ggplot2::element_text(size=15, face="bold"),
          axis.text = ggplot2::element_text(size=15, face="bold"),
          legend.text=ggplot2::element_text(size=15),
          legend.title=ggplot2::element_text(size=15, face="bold"),
          legend.position = "bottom",
          panel.background = ggplot2::element_blank())+
            ggplot2::xlab("Force-directed layout_1")+ggplot2::ylab("Force-directed layout_2")
  return(p)
}
######################################################################################################################
################################################ internal functions ##################################################
######################################################################################################################
.generate_snn_graph <- function(factor.list, knn.param=20, prune.param=1/5){
  n_layer = length(factor.list)
  # check the dimension of each of the low-dimensional embedding
  ncells = unlist(lapply(factor.list, nrow))
  ncell=ncells[1]
  if (sum(abs(ncells-ncell))>0){
    stop("Number of cells in each layer are not the same.")
  }
  Jaccard_Index =  prune.param/(2 - prune.param)
  
  if(is.null(names(factor.list))){
    names(factor.list)=sprintf("layer%s", 1:n_layer)
  }
  snn.list = list()
  for (i in 1:n_layer){
    if (is.null(rownames(factor.list[[i]]))){
      rownames(factor.list[[i]])=sprintf("Cell%s", 1:ncells[1])
    }
    nn_graph <- Seurat::FindNeighbors(factor.list[[i]],  k.param=knn.param, prune.SNN = Jaccard_Index, verbose = F)
    snn.list[[i]] = nn_graph$snn
    snn.list[[i]][snn.list[[i]]>0]=1
    diag(snn.list[[i]])=0
  }
  names(snn.list) = names(factor.list)
  return(snn.list)
}

.generate_mlg<-function(snn.list){
  n_layer = length(snn.list)
  ncell=nrow(snn.list[[1]])
  mlg = Matrix::sparseMatrix(dims = c(ncell, ncell), i={}, j={})
  for (i in 1:n_layer){
    mlg = mlg + snn.list[[i]]
  }
  mlg[mlg>0]=1
  return(mlg)
}
