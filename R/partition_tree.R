#' Using non-backtracking strategy for stopping rule and build the tree
#'
#' \code{build_tree} returns a dataframe
#'
#' @param f It is the adjencency matrix
#' @param xi.loc.labels lists of index of each cluster
#' @param ncl number of clusters
#' @param cl.labels a vector of index of individuals
#' @return Lists of cluster results
#' @examples
#' cluster_result <- build_tree(adj,xi.loc.labels=list(), ncl=0, cl.labels=1:339,n.min=25,D=NULL)

build_tree <- function(f,xi.loc.labels, ncl, cl.labels,n.min=25,D=NULL){
  nisol = which(rowSums(f) > 0)
  isol = which(rowSums(f) == 0)
  cl.labels.full <- cl.labels
  ## sanity check -- do not start if there are too many isolated nodes
  if((length(nisol)<=8)||(length(isol)>=5*length(nisol))||(length(nisol)<2*n.min)){
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    tree.path <- c("",as.character(cl.labels))
    mod.path <- c(0,rep(0,length(cl.labels)))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels[isol]))
    print('Too few connected nodes, not even started!')
    return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path,mod.path=mod.path))
  }
  #print(paste(length(isol),"isolated nodes",cl.labels[isol]))
  all.deg = rowSums(f)
  f = f[nisol, nisol] ### only focus on non-isolated nodes - isolated nodes with be attached to the smaller child
  cl.labels.full <- cl.labels
  all.deg = all.deg[nisol]
  cl.labels = cl.labels[nisol]
  n1 = dim(f)[1]

  K = 2
  d <- colSums(f)
  n <- nrow(f)
  I <- as(diag(rep(1,n)),"dgCMatrix")
  D <- as(diag(d),"dgCMatrix")
  r <- sqrt(sum(d^2)/sum(d)-1)
  B <- as(matrix(0,nrow=2*n,ncol=2*n),"dgCMatrix")
  B[(n+1):(2*n),1:n] <- -1*I
  B[1:n,(n+1):(2*n)] <- D-I
  B[(n+1):(2*n),(n+1):(2*n)] <- f
  ss <- Re(eigs(B,k=2,which="LM")$values)
  split.flag <- sum(abs(ss)>r)>1
  if(split.flag) {
    #f.adj <- f + 3*mean(all.deg)/nrow(f)
    #f.adj <- f + 0.1/nrow(f)
    g <- graph_from_adjacency_matrix(as.matrix(f), weighted=T, mode="undirected")
    g <- simplify(g)
    adj <- get.adjacency(g,type = "both",attr = "weight")
    deg.adj <- diag(strength(g))+diag(rep(0.1,nrow(f)))

    l_matrix <- solve(sqrtm(deg.adj))%*%adj%*%solve(sqrtm(deg.adj))
    adj <- NULL
    deg.adj <- NULL
    eval <- eigs_sym(l_matrix,2,which = "LM")
    embed.Y <- data.frame(eval$vectors)
    #embed.Y <- data.frame(embed.s$X[,1:2])
    row.names(embed.Y) <- V(g)$name
    clustering <- kmeans(embed.Y,centers=2,iter.max=30,nstart=10)$cluster
    clus = list(clustering=clustering)
    xi.f = clus$clustering
    xi.labels = lapply(1:2, function(x){which(xi.f == x)})
    smaller.cluster <- xi.labels[[which.min(sapply(xi.labels,length))]]
    f1 <- f[smaller.cluster,smaller.cluster]
    a1.labels <- cl.labels[smaller.cluster]
    if(length(dim(f1)) > 0) {
      n1 <- dim(f1)[1]
    } else if(length(f1) > 0) { ### case when only 1 node is in this cluster, make it 2, so the later on rank check code still works
      n1 <- 1
    } else {
      n1 = 0
    }

    if(n1 > 2*n.min) { ## only do further clustering on cluster larger  2*n.min
      res = build_tree(f1,xi.loc.labels, ncl, a1.labels,n.min,D=D-1)
      xi.loc.labels = res$xi.loc.labels
      ncl = res$ncl
      L.tree.path <- res$tree.path
      if(length(isol)>0){
        xi.loc.labels[[ncl]] = c(xi.loc.labels[[ncl]],cl.labels.full[isol]) ### attached the isolated nodes in this level with the clusters under the smaller split

        path.head <- L.tree.path[length(L.tree.path)]
        path.head <- gsub('[[:digit:]]+', '', path.head)
        iso.path <- paste(path.head,cl.labels.full[isol],sep="")
        L.tree.path <- c(L.tree.path,iso.path)
      }

      L.mod.path <- res$mod.path

    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a1.labels
      if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],xi.loc.labels[[ncl]])
      L.tree.path <- as.character(xi.loc.labels[[ncl]])
      L.mod.path <- rep(0,length(xi.loc.labels[[ncl]]))
      #print('Too small left cluster, Branch End')
    }
    f2 = f[-smaller.cluster, -smaller.cluster]

    a2.labels = cl.labels[-smaller.cluster]
    if(length(dim(f2)) > 0) {
      n1 <- nrow(f2)

    } else if(length(f2) > 0) {
      n1 <- 1
    } else {
      n1 <- 0
    }
    if(n1 > 2*n.min) {
      res = build_tree(f2,xi.loc.labels, ncl, a2.labels,n.min,D=D-1)
      xi.loc.labels = res$xi.loc.labels
      R.tree.path <- res$tree.path
      R.mod.path <- res$mod.path
      ncl = res$ncl

    } else {
      ncl = ncl + 1
      xi.loc.labels[[ncl]] = a2.labels
      R.tree.path <- as.character(a2.labels)
      R.mod.path <- rep(0,length(a2.labels))
      #print('Too small right cluster, Branch End')
    }
    L.tree.path <- paste("L",L.tree.path,sep="/")
    R.tree.path <- paste("R",R.tree.path,sep="/")
    tree.path <- c("",L.tree.path,R.tree.path)

  } else {
    ncl = ncl + 1
    xi.loc.labels[[ncl]] = cl.labels
    if(length(isol)>0) xi.loc.labels[[ncl]] = c(cl.labels.full[isol],cl.labels)
    tree.path <- c("",as.character(cl.labels))
    if(length(isol)>0) tree.path <- c(tree.path,as.character(cl.labels.full[isol]))
    #print('One cluster, Branch End, not even started!')
  }
  return(list(xi.loc.labels = xi.loc.labels, ncl = ncl,tree.path=tree.path))
}

