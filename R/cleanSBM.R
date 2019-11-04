#' Run SBM and clean SBM results
#'
#' \code{cleanSBM} returns lists of sbm results
#'
#' @param n The number of times repeating running sbm
#' @param m The number of cores required, notice n/m should be integer
#' @param dat The data which satisfy the requirements of the formate of runnning sbm
#' @return Lists of sbm results
#' @examples
#' sbmresults <- cleanSBM(1000,10,pheno_dat)


cleanSBM <- function(n, m, dat){
  phenome_dat_long_list <- lapply(1:(n/m),function(x) dat)
  mycluster <- makeCluster(m,type = "PSOCK")
  registerDoParallel(mycluster)
  getDoParWorkers()
  getDoParName()
  start <- Sys.time()
  node_to_cluster <- foreach(i=1:m) %dopar% {
    library(dplyr)
    devtools::load_all("/home/siwei/pheSBMR")
    lapply(phenome_dat_long_list, run_sbm)
  }
  stopCluster(mycluster)
  end <- Sys.time()
  time <- end-start

  for (i in 1:m) {
    miss <- which(sapply(node_to_cluster[[i]],function(x) nrow(x)==0))
    while(length(miss)!=0){
      
        for (j in 1:length(miss)) {
          node_to_cluster[[i]][[miss[j]]] <-run_sbm(phenome_dat_long)
        }
      
      miss <- which(sapply(node_to_cluster[[i]],function(x) nrow(x)==0))
    }
  }

  sbmresults <- list()
  k <- 1
  for (i in 1:m) {
    for (j in 1:(n/m)) {
      sbmresults[[k]] <- node_to_cluster[[i]][[j]]
      k <- k+1
    }
  }

  return(sbmresults)
}

