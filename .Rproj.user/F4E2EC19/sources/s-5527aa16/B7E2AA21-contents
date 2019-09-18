#' Create similarity matrix based on each sbm result
#'
#' \code{create_affinity} returns an affinity matrix
#'
#' @param i The index of which sbm result in the dat
#' @param j The number of individuals in each sbm result
#' @param dat The data which get from cleanSBM function
#' @return An affinity matrix
#' @examples
#' A <- create_affinity(1,339,sbmresults)


create_affinity <- function(i,j,dat){
  A <- matrix(0,j,j)
  cluster_results_levels <- dat[[i]] %>%
    dplyr::filter(type=="grid")

  #for (i in 0:(nlevels(as.factor(cluster_results_levels$level))-1)) {
  cluster_results <- cluster_results_levels %>%
    dplyr::filter(level==0) %>%
    dplyr::select(node,cluster)
  cluster_results$cluster <- as.factor(cluster_results$cluster)
  namesorder <- cluster_results$node

  #create affinity matrix
  affinitymatrix_list <- list()
  rownames=c()
  colnames=c()
  for (l in 1:nlevels(cluster_results$cluster)) {
    cluster_group <- subset(cluster_results,cluster==levels(cluster_results$cluster)[l])
    affinitymatrix <- matrix(1,nrow(cluster_group),nrow(cluster_group))
    rownames <- c(rownames,cluster_group$node)
    colnames <- c(colnames,cluster_group$node)
    affinitymatrix_list[[l]] <- affinitymatrix
  }
  #combine each cluster matrix into one matrix
  library(Matrix)
  affinitymatrix <- as.matrix(bdiag(affinitymatrix_list))
  rownames(affinitymatrix) <- rownames
  colnames(affinitymatrix) <- colnames
  #affinitymatrix <- affinitymatrix[namesorder,namesorder]
  #affinitymatrix <- as(affinitymatrix[namesorder,namesorder],"sparseMatrix")
  #make addition of the matrix
  #A <- A+1/(nlevels(as.factor(cluster_results_levels$level)))*affinitymatrix[namesorder,namesorder]
  A <- A+affinitymatrix[namesorder,namesorder]
  #}
  print(paste(j,"affinity matrix end",sep = ""))
  return(A <- A)
}

