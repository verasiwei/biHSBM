#' Build the Membership matrix
#'
#' \code{membership} returns a dataframe
#'
#' @param dat It is the result from partition_tree
#' @param A adjencency matrix

#' @return A dataframe of membership matrix
#' @examples
#' member <- membership(cluster_result,A)


membership <- function(dat,A){

  cluster_id_func <- function(k,tree.path,A){
    node_index <- which(sapply(1:length(cluster_result$tree.path),function(j) substr(cluster_result$tree.path[j],1,2*k)==tree.path))

    node_number <- as.integer(gsub('[^[:digit:].]+', '', cluster_result$tree.path[node_index]))
    node_number <- node_number[which(!is.na(node_number))]

    cluster_id <- rownames(A)[as.integer(node_number)]
    return(cluster_id <- cluster_id)
  }

  tree_path_type <- unique(gsub('[[:digit:]]+', '', dat$tree.path))
  tree_path_type <- tree_path_type[c(-1)]
  level_total <- c()
  for (i in 1:length(tree_path_type)) {
    level_total <- c(level_total,str_count(tree_path_type[i],"/"))
  }
  level_total <- max(level_total)


  grid_list <- list()
  grid <- list()
  for (i in 1:level_total) {
    for (j in tree_path_type[which(str_count(tree_path_type,"/")==i)]) {
      grid[[j]] <- cluster_id_func(i,j,A)
    }
    grid_list[[i]] <- grid
    grid <- list()
  }

  cluster_dat <- data.frame("grid"=c(grid_list[[1]][[1]],grid_list[[1]][[2]]))
  for (i in 1:level_total) {
    level <- data.frame()
    for (j in 1:length(grid_list[[i]])) {
      dat <- data.frame("grid"=grid_list[[i]][[j]],"group"=rep(names(grid_list[[i]][j]),length(grid_list[[i]][[j]])))
      level <- rbind(level,dat)
    }
    cluster_dat <- join(cluster_dat,level,by="grid")
  }
  for (i in 2:(level_total+1)) {
    colnames(cluster_dat)[i] <- paste("level.",i-1,sep = "")
  }

  for (i in 2:(level_total+1)) {
    colnames(cluster_dat)[i] <- paste("level-",i-1,sep = "")
    cluster_dat[,i] <- gsub("L",0,cluster_dat[,i])
    cluster_dat[,i] <- gsub("R",1,cluster_dat[,i])
    cluster_dat[,i] <- str_remove_all(cluster_dat[,i],"/")
    cluster_dat[,i] <- str_sub(cluster_dat[,i],-1,-1)
    cluster_dat[is.na(cluster_dat)] <- "x"
    cluster_dat[,c(-1)] <- lapply(cluster_dat[,c(-1)],factor)
  }
  return(cluster_dat <- cluster_dat)
}


