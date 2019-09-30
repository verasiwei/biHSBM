#' Create aggregation matrix and apply different methods
#'
#' \code{create_affinity} returns an affinity matrix
#'
#' @param N Aggregation times
#' @param M The number of total sbm results
#' @param n The number of individuals
#' @param seed A vector of the seed you want to set for each aggregation time
#' @param rep The number of repeats for each aggregation time
#' @param method The method, kmeans, hclust or eigensplit
#' @param bucket The bucket where you save the sbm affinity matrix
#' @param reg_par The parameter value of regularization
#' @param var The percentage explaine the variation
#' @param individual The individual id in the cluster that you want to see consistency
#' @return
#' @examples
#' aug_result <- augmentation_matrix(20,2000,339,c(1:10),10,c("kmeans","hclust","eigensplit"),"jak2",0.1,0.9)

augmentation_matrix <-function(N,M,n,seed,rep,method,bucket,reg_par,var,individual){

  augmentation_index <- function(N,M,seed,rep){
    index_list <- list()
    index_rest_list <- list()
    # set.seed(seed)
    # index_list[[1]] <- sample(1:M,N,replace = FALSE)
    # index_rest <- c(1:M)[!(c(1:M) %in% index_list[[1]])]
    # index_rest_list[[1]] <- sample(index_rest,N,replace = FALSE)
    for (r in 1:rep) {
      set.seed(seed[r])
      index_list[[r]] <- sample(1:M,N,replace = FALSE)
      index_rest <- c(1:M)[!(c(1:M) %in% index_list[[r]])]
      set.seed(seed[r])
      index_rest_list[[r]] <- sample(index_rest,N,replace = FALSE)
      #index_rest <- index_rest[!(index_rest %in% index_list[[i]])]
    }
    index_list_total <- list(index_list,index_rest_list)
    return(index_list_total <- index_list_total)
  }

    embed.kmeans_list <- list()
    embed.kmeans2_list <- list()
    embed.hc_list <- list()
    embed.hc2_list <- list()
    embed.eigen_list <- list()
    for (i in 1:2) {
      embed <- list()
      embed.kmeans <- list()
      embed.kmeans2 <- list()
      embed.hc <- list()
      embed.hc2 <- list()
      embed.eigen <- list()
      for (k in 1:rep) {
        A <- matrix(0,n,n)
        for (j in augmentation_index(N,M,seed,rep)[[i]][[k]]) {
          A <- A+s3readRDS(object = paste("sbmmatrix_",j,".rds",sep = ""),bucket = bucket)
        }
        print(paste("finish augment",i,"_",k,sep = ""))
        A <- A[which(rownames(A) %in% individual),which(colnames(A) %in% individual)]
        #A <- A+as.numeric(reg_par)/nrow(A)
        g <- graph_from_adjacency_matrix(A, weighted=T, mode="undirected")
        g <- simplify(g)
        adj <- get.adjacency(g,type = "both",attr = "weight")
        deg.adj <- diag(strength(g))+diag(rep(as.numeric(reg_par),nrow(A)))
        A <- NULL
        #diag(adj.adj) <- 0
        if(method %in% c("eigensplit","all")){
          adj.adj <- adj+as.numeric(reg_par)/nrow(adj)
          eval <- eigs_sym(adj.adj,2,which = "LM")
          embed.eigen[[k]] <- ifelse(eval$vectors[,2]<0,1,2)
          adj.adj <- NULL
        }

        #deg.adj <- deg+diag(rep(reg_par,nrow(deg)))
        #all.deg <- rowSums(adj)
        #adj.reg <- adj + 0.1*mean(all.deg)/nrow(adj)
        #g <- graph_from_adjacency_matrix(adj.reg, weighted=T, mode="undirected")
        #adj <- NULL
        #g <- simplify(g)
        l_matrix <- solve(sqrtm(deg.adj))%*%adj%*%solve(sqrtm(deg.adj))
        deg.adj <- NULL
        adj <- NULL
        #embed.s <- embed_laplacian_matrix(g, no=10, type='DAD',scaled = FALSE)
        eval <- eigs_sym(l_matrix,10,which = "LM")
        eval.adj <- eval$vectors%*%diag(eval$values)

        embed.slist <- eval$values[1]
        m=1
        while(embed.slist<as.numeric(var)*sum(eval$values)){
          m=m+1
          embed.slist <- embed.slist+eval$values[m]
        }
        embed.Y <- data.frame(eval.adj[,1:m])
        embed.Y2 <- data.frame(eval.adj[,1:2])
        row.names(embed.Y) <- V(g)$name
        row.names(embed.Y2) <- V(g)$name
        g <- NULL
        if(method %in% c("kmeans","all")){
          #embed.kmeans[[k]] <- kmeans(embed.Y,centers=2,iter.max=30,nstart=10)$cluster
          #embed.kmeans2[[k]] <- kmeans(embed.Y2,centers=2,iter.max=30,nstart=10)$cluster
          embed.kmeans[[k]] <- pam(embed.Y,k=2)$clustering
          embed.kmeans2[[k]] <- pam(embed.Y2,k=2)$clustering
        }
        if(method %in% c("hclust","all")) {
          embed.hc_result <- hclust(dist(embed.Y), method="ward.D2")
          embed.hc[[k]] <- cutree(embed.hc_result,k=2)
          embed.hc2_result <- hclust(dist(embed.Y2), method="ward.D2")
          embed.hc2[[k]] <- cutree(embed.hc2_result,k=2)
        }
        print(paste("finish clustering",i,"_",k,sep = ""))
      }
      if(method=="kmeans"){
        embed.kmeans_list[[i]] <- embed.kmeans
        embed.kmeans2_list[[i]] <- embed.kmeans2
      }
      else if(method=="hclust") {
        embed.hc_list[[i]] <- embed.hc
        embed.hc2_list[[i]] <- embed.hc2
      }
      else if(method=="eigensplit"){
        embed.eigen_list[[i]] <- embed.eigen
      }
      else if(method=="all"){
        embed.kmeans_list[[i]] <- embed.kmeans
        embed.kmeans2_list[[i]] <- embed.kmeans2
        embed.hc_list[[i]] <- embed.hc
        embed.hc2_list[[i]] <- embed.hc2
        embed.eigen_list[[i]] <- embed.eigen
      }
    }

    if(method=="kmeans"){
      embed_total <- list("embed.kmeans_list"=embed.kmeans_list,"embed.kmeans2_list"=embed.kmeans2_list)
      } else if(method=="hclust") {
        embed_total <- list("embed.hc_list"=embed.hc_list,"embed.hc2_list"=embed.hc2_list)
        } else if(method=="eigensplit"){
          embed_total <- embed.eigen_list
        } else if(method=="all"){
          embed_total <-list("embed.kmeans_list"=embed.kmeans_list,"embed.kmeans2_list"=embed.kmeans2_list,
                             "embed.hc_list"=embed.hc_list,"embed.hc2_list"=embed.hc2_list,
                             "embed.eigen_split"=embed.eigen_list)
        }

  return(embed_total)

}


