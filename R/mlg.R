######################################################################################################################
######################################################################################################################
#' MLG clustering 
#' This function takes a list of low dimensional embedding data as input, and performs MLG clustering.
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph
#' @param prune.param The prune parameter forSNN graph ??
#' @param cluster.resolution The resolution number of modularity maximization
#' @param cluster.algorithm The clustering algorithm. 1--Louvain algorithm; 2--Louvain algorithm with multilevel refinement; 3--SLM algorithm.
#' @return a vector of cluster labels
#' @export 
mlg_cluster <- function(factor.list, knn.param = 20, prune.param = 1/5, cluster.resolution, cluster.algorithm = 1){
  snn.list = .generate_snn_graph(factor.list, knn.param, prune.param)
  mlg = .generate_mlg(snn.list)
  cluster = Seurat::FindClusters(mlg, resolution= cluster.resolution, algorithm =cluster.algorithm)[[1]]
  return(cluster)
}

#cc = mlg_cluster(factor.list = factor.list, knn.param = 20, prune.param = 1/5, 
#                       cluster.algorithm = 1, cluster.resolution = .3)

######################################################################################################################
######################################################################################################################
#' This function takes a list of low dimensional embedding data as input, construct a SNN graph for each low dimensional embeddings, and then plot a heatmap of the proportion of overlapping edges between each pair of SNN graphs.
#'
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph
#' @param prune.param The prune parameter forSNN graph ??
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
    ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::labs(fill="Proportion of \noverlapped edges    ")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold",
                                     size=10),
          axis.text.y = ggplot2::element_text(face="bold",
                                     size=10),
          legend.text=ggplot2::element_text(size=10),
          legend.title=ggplot2::element_text(size=12, face="bold"),
          legend.position = "bottom",
          panel.background = ggplot2::element_blank())+ 
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)))
  return(p)
  
}
#prop.overlap.edges(factor.list, knn.param, prune.param)
######################################################################################################################
######################################################################################################################
#' This function takes a list of low dimensional embedding data, and true cell type labels as input, and 
#' compute the graph signal to noise ratio of the SNN graphs constructed with each low dimension embeddding
#' and the MLG graph.
#'
#' @param factor.list A list that contains several low dimensional embedding data.
#' @param knn.param The number of neighbors used to construct k-nearest-neighbot graph
#' @param prune.param The prune parameter forSNN graph ??
#' @param cell_label True cell type label
#' @return a vector of signal to noise ratio, labeled by graph name
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
#graph_signal_noise_ratio(factor.list , cell_label =cell_label )
######################################################################################################################
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

#############################################################
## compute optimal permutations that most cluster labels match
# .recode<-function(vec1, level1, level2){
#   vec2=rep(0,length(vec1))
#   for(i in 1:length(level1)){
#     vec2[vec1==level1[i]]=level2[i]
#   }
#   return(vec2)
# }
# 
# opt.label.permutation <- function(true.label, predict.label){
#   true.label = as.character(true.label)
#   predict.label = as.character(predict.label)
#   num_levels_true = length(unique(true.label))
#   num_levels_pred = length(unique(predict.label))
#   
#   
#   nperm = min(500, factorial(num_levels_pred))
#   permute_label = matrix("NA", nrow=length(predict.label), nperm)
#   if (nperm<500){
#     perm_matrix = gtools::permutations(num_levels_pred, num_levels_pred)
#   }
#   
#   for ( j in 1:nperm){
#     if (nperm<500){
#       recode_pred = .recode(predict.label, unique(predict.label), perm_matrix[j,])
#     }else{
#       recode_pred = .recode(predict.label, unique(predict.label), sample(num_levels_pred))
#     }
#   
#     M=matrix(0, nrow = num_levels_true, ncol = num_levels_pred)
#     
#     for (i in 1:num_levels_true){
#       for (k in 1:num_levels_pred){
#         M[i,k] = mean(recode_pred[true.label==unique(true.label)[i]]==k)
#       }
#     }
# 
#     idx = 1:num_levels_pred
#     for (i in 1:num_levels_pred){
#       idx[i] = which.max(M[,i])
#       M[idx[i],]=-1
#     }
#     permute_label[,j] = rep("NA",length(predict.label))
#     for (i in 1:num_levels_pred){
#       permute_label[recode_pred==i,j]=unique(true.label)[idx[i]]
#     }
#   }
#   
#   match.ratio.vector<-apply(permute_label, 2, function(i){mean(i==true.label)})
#     
#   return(list(match.ratio=max(match.ratio.vector), label=permute_label[, which.max(match.ratio.vector)]))
# }
# 
# #opt.label.permutation(cluster_result[[2]], cluster_result[[1]])[[1]]
