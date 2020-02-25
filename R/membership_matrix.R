#' Build the Membership matrix
#'
#' \code{membership} returns a dataframe
#'
#' @param dat It is the result from build_tree
#' @param A adjencency matrix

#' @return A dataframe of membership matrix and a dataframe for plotting
#' @examples
#' member <- membership_dat(cluster_result)


membership_str <- function(dat,A){
  tree_path_type <- unique(gsub('[[:digit:]]+', '', dat$tree.path))
  tree_path_type <- tree_path_type[c(-1)]
  level_total <- c()
  for (i in 1:length(tree_path_type)) {
    level_total <- c(level_total,str_count(tree_path_type[i],"/"))
  }
  level_total <- max(level_total)

  cluster_id_func <- function(dat,k,tree.path,A){
    node_index <- which(sapply(1:length(dat$tree.path),function(j) substr(dat$tree.path[j],1,2*k)==tree.path))

    node_number <- as.integer(gsub('[^[:digit:].]+', '', dat$tree.path[node_index]))
    node_number <- node_number[which(!is.na(node_number))]

    cluster_id <- rownames(A)[as.integer(node_number)]
    return(cluster_id <- cluster_id)
  }
  grid_list <- list()
  grid <- list()
  for (i in 1:level_total) {
    for (j in tree_path_type[which(str_count(tree_path_type,"/")==i)]) {
      grid[[j]] <- cluster_id_func(dat,i,j,A)
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

  nnodes=0
  for (i in 2:ncol(cluster_dat)) {
    nnodes <- nnodes+nlevels(cluster_dat[,i])
  }

  n <- 1
  for (i in 2:ncol(cluster_dat)) {
    for (j in 1:nlevels(cluster_dat[,i])) {
      levels(cluster_dat[,i])[j] <- n
      n <- n+1
    }
  }

  weight <- c()
  text <- c(0)
  n <- 1
  for (i in 2:(ncol(cluster_dat)-1)) {
    for (j in 1:nlevels(cluster_dat[,i])) {
      elements <- cluster_dat[which(cluster_dat[,i]==levels(cluster_dat[,i])[j]),i+1]
      if(!is.na(elements[1])){
        weight <- c(weight,sum((cluster_dat[,i]==as.character(n)),na.rm = TRUE))
      }
      if(!is.na(elements[1])){
        text <- c(text,levels(cluster_dat[,i])[j])
      }
      n <- n+1
    }
  }
  n <- 1
  for (i in 2:(ncol(cluster_dat)-1)) {
    for (j in 1:nlevels(cluster_dat[,i])) {
      elements <- cluster_dat[which(cluster_dat[,i]==levels(cluster_dat[,i])[j]),i+1]
      if(is.na(elements[1])){
        weight <- c(weight,sum((cluster_dat[,i]==as.character(n)),na.rm = TRUE))
      }
      if(is.na(elements[1])){
        text <- c(text,levels(cluster_dat[,i])[j])
      }
      n <- n+1
    }
  }
  text <- c(text,levels(cluster_dat[,ncol(cluster_dat)]))
  for (j in 1:nlevels(cluster_dat[,(ncol(cluster_dat)-1)+1])) {
    weight <- c(weight,sum((cluster_dat[,(ncol(cluster_dat)-1)+1]==as.character(n)),na.rm = TRUE))
    n <- n+1
  }
  weight <- insert(weight,1,sum(weight[1]+weight[2]))


  endnode <- c()
  for (i in 3:ncol(cluster_dat)) {
    endnode <- c(endnode,as.character(unique(cluster_dat[which(is.na(cluster_dat[,i])),i-1])))
    endnode <- endnode[!is.na(endnode)]
  }

  parent <- c(0,0)
  i=1
  while (length(parent) < nnodes) {
    if(i %in% as.numeric(endnode)){
      i <- i+1
    }else{
      parent <- c(parent,rep(i,2))
      i <- i+1
    }
  }
  dat <- data.frame("parent"=parent,"node"=c(1:nnodes),"text"=rep(c(0,1),nnodes/2))
  return(list("cluster_dat"=cluster_dat,"plot_dat"=dat,"weight"=weight))
}


