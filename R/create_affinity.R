#' Create similarity matrix based on each sbm result
#'
#' \code{create_affinity} returns an affinity matrix
#'
#' @param i The index of which sbm result in the dat
#' @param j The number of individuals in each sbm result
#' @param dat The data which get from cleanSBM function
#' @return An affinity matrix
#' @examples
#' A <- create_affinity(1,339,sbmresults,TRUE,0.8)


create_affinity <- function(i,j,dat,resample,p){
  A <- matrix(0,j,j)
  cluster_results_levels <- dat[[i]] %>%
    dplyr::filter(type=="grid")

  #for (i in 0:(nlevels(as.factor(cluster_results_levels$level))-1)) {
  cluster_results <- cluster_results_levels %>%
    dplyr::filter(level==0) %>%
    dplyr::select(node,cluster)
  cluster_results$cluster <- as.factor(cluster_results$cluster)
  namesorder <- as.character(cluster_results$node)

  #create affinity matrix
  affinitymatrix_list <- list()
  identitymatrix_list <- list()
  rownames=c()
  colnames=c()
  if(resample == FALSE){
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
    rownames(affinitymatrix) <- as.character(rownames)
    colnames(affinitymatrix) <- as.character(colnames)
    #affinitymatrix <- affinitymatrix[namesorder,namesorder]
    #affinitymatrix <- as(affinitymatrix[namesorder,namesorder],"sparseMatrix")
    #make addition of the matrix
    #A <- A+1/(nlevels(as.factor(cluster_results_levels$level)))*affinitymatrix[namesorder,namesorder]
    A <- A+affinitymatrix[namesorder,namesorder]
    identitymatrix <- NULL
  } else {
     #sample from original data
     sample_node <- sample(namesorder,size = p*j,replace = FALSE)
     sample_node_rest <- namesorder[which(!(namesorder %in% sample_node))]
     cluster_results <- cluster_results[which(cluster_results$node %in% sample_node),]
     #create identity matrix
     identitymatrix_list[[1]] <- matrix(1,length(sample_node),length(sample_node))
     identitymatrix_list[[2]] <- matrix(0,length(sample_node_rest),length(sample_node_rest))
     identitymatrix <- as.matrix(bdiag(identitymatrix_list))
     rownames(identitymatrix) <- as.character(c(sample_node,sample_node_rest))
     colnames(identitymatrix) <- as.character(c(sample_node,sample_node_rest))
     identitymatrix <- identitymatrix[namesorder,namesorder]
     
     for (l in 1:nlevels(cluster_results$cluster)) {
       cluster_group <- subset(cluster_results,cluster==levels(cluster_results$cluster)[l])
       affinitymatrix <- matrix(1,nrow(cluster_group),nrow(cluster_group))
       rownames <- as.character(c(rownames,cluster_group$node))
       colnames <- as.character(c(colnames,cluster_group$node))
       affinitymatrix_list[[l]] <- affinitymatrix
     }
     affinitymatrix_list[[l+1]] <- identitymatrix_list[[2]]
     
     #combine each cluster matrix into one matrix
     library(Matrix)
     affinitymatrix <- as.matrix(bdiag(affinitymatrix_list))
     rownames(affinitymatrix) <- as.character(c(rownames,sample_node_rest))
     colnames(affinitymatrix) <- as.character(c(colnames,sample_node_rest))
     #affinitymatrix <- affinitymatrix[namesorder,namesorder]
     #affinitymatrix <- as(affinitymatrix[namesorder,namesorder],"sparseMatrix")
     #make addition of the matrix
     #A <- A+1/(nlevels(as.factor(cluster_results_levels$level)))*affinitymatrix[namesorder,namesorder]
     A <- A+affinitymatrix[namesorder,namesorder]
    }
  print(paste(i,"affinity matrix end",sep = ""))
  return(list(A <- A,identitymatrix <- identitymatrix))
}

