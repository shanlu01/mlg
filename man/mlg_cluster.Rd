% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlg.R
\name{mlg_cluster}
\alias{mlg_cluster}
\title{MLG clustering
This function takes a list of low dimensional embedding data as input, and performs MLG clustering.}
\usage{
mlg_cluster(
  factor.list,
  knn.param = 20,
  prune.param = 1/5,
  cluster.resolution,
  cluster.algorithm = 1
)
}
\arguments{
\item{factor.list}{A list that contains several low dimensional embedding data.}

\item{knn.param}{The number of neighbors used to construct k-nearest-neighbot graph}

\item{prune.param}{The prune parameter for SNN graph. There is an edge between cell i and j in the SNN graph, if the number of common neighbors between i and j exceeds the product of knn.param and knn.param.}

\item{cluster.resolution}{The resolution number of modularity maximization.}

\item{cluster.algorithm}{The clustering algorithm. 1--Louvain algorithm; 2--Louvain algorithm with multilevel refinement; 3--SLM algorithm.}
}
\value{
A vector of cluster labels.
}
\description{
MLG clustering
This function takes a list of low dimensional embedding data as input, and performs MLG clustering.
}
