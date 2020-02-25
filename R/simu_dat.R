#' Generate Simulation Data
#'
#' \code{simu_dat} returns lists of several dataframes
#'
#' @param N The number of a nodes clusters.
#' @param M The number of a nodes
#' @param n The number of stable b nodes that plant the pattern of a nodes clusters
#' @param p1 The minimum probability of turning on b node in each a node cluster(randomly generating 0-1 vector based on this uniform distribution on the interval from p1 to p2)
#' @param p2 The maximum probability of turning on b node in each a node cluster
#' @param num_noise The number of noise nodes
#' @param noise_p The individual lambda for connections between the noise b-nodes and any a-node
#' @param planted_p_on The lambda for a planted pattern b-node that has been turned on in a pattern.
#' @param planted_p_off The lambda for a planted pattern b-node that has been turned off in pattern.
#' @param ind Logic TRUE or FALSE represent whether add totally independent a cluster
#' @param ind_n The number of stable b nodes you want to add for independent a cluster
#' @param ind_sub The sample size of the independent cluster you want to add
#' @return Lists of lists of dataframes
#' @examples
#' dat <- my_patterns_function(10,400,10,60,0.4,0.4,1000,0.01,0.95,0.05,FALSE,5,20)
#' 


#created patterns of noise function
my_patterns_function <- function(N,M,sd,n,p1,p2,num_noise_node,noise_p,planted_p_on,planted_p_off,ind,ind_n,ind_sub){
  #function to generate N random integers that sum to M in R(generate size in )
rand_vect <- function(N, M, sd, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

  my_patterns <- t(as.numeric(rbinom(n,1,runif(1,min = p1,max = p2))))
    for (i in 2:N) {
      set.seed(i)
      pattern_subject <- t(as.numeric(rbinom(n,1,runif(1,min = p1,max = p2))))
      while (pattern_subject==my_patterns[i-1,]) {
        pattern_subject <- t(as.numeric(rbinom(n,1,runif(1,min = p1,max = p2))))
      }
      my_patterns <- rbind(my_patterns,pattern_subject)
      }
    my_patterns <- data.frame(cbind(my_patterns,rand_vect(N,M,sd)))
    colnames(my_patterns)[ncol(my_patterns)] <- "size"
    colnames(my_patterns)[1:n] <- paste("b",1:n,sep="")
#  } else {
#    for (i in 2:N) {
#      set.seed(i)
#      pattern_subject <- as.numeric(rbinom(n,1,runif(1,min = p1,max = p2)))
#      while (pattern_subject==my_patterns[i-1,]) {
#        pattern_subject <- as.numeric(rbinom(n,1,runif(1,min = p1,max = p2)))
#      }
#      my_patterns <- rbind(my_patterns,pattern_subject)
#    }
#    add_patterns <- t(rep(1,ind_n))
#    my_patterns <- data.frame(as.matrix(bdiag(my_patterns,add_patterns)))
#    size <- c(rand_vect(N,M),ind_sub)
#    my_patterns <- data.frame(cbind(my_patterns,size))
#    #rownames(my_patterns)[1] <- "1"
#    #rownames(my_patterns) <- as.integer(rownames(my_patterns))
#    colnames(my_patterns)[ncol(my_patterns)] <- "size"
#    colnames(my_patterns)[1:(n+ind_n)] <- paste("b",1:(n+ind_n),sep="")
#  }
  
  #2nd step
  #planted pattern
  planted_model_params <- setup_planted_pattern_model(
  my_patterns, 
  num_noise_nodes = num_noise_node,
  noise_p = noise_p,
  planted_p_on = planted_p_on,
  planted_p_off = planted_p_off
)

draw_from_planted <- planted_model_params %$%
  draw_from_model(
    b_a, b_b, Lambda, 
    binary_connections = TRUE, 
    a_name = "Subjects", 
    b_name = "Phecodes" )

#3rd step 
##generate simulation data
###change the index of noise b nodes
noise_node <- c()
for (i in (ncol(my_patterns)):(num_noise_node+n)) {
  noise_node <- c(noise_node,rep(i,M))
}
draw_from_planted$b_group[(M*n+1):nrow(draw_from_planted)] <- noise_node

pheno_dat_long <- draw_from_planted[draw_from_planted$num_edges==1,1:2]
colnames(pheno_dat_long) <- c("grid","phecode")
#add independent cluster
if(ind == TRUE){
ind_long <- data.frame("grid" = c((M+1):(M+ind_sub)),
                         "phecode" = rep(n+num_noise_node+1,ind_sub))
for (i in (n+num_noise_node+2):(n+num_noise_node+ind_n)) {
  ind_long <- rbind(ind_long,data.frame("grid" = c((M+1):(M+ind_sub)),
                         "phecode" = rep(i,ind_sub)))
}
pheno_dat_long <- rbind(pheno_dat_long,ind_long)
}


return(list("draw_from_planted"=draw_from_planted,
            "planted_model_params"=planted_model_params,
            "my_patterns"=my_patterns,
            "pheno_dat_long"=pheno_dat_long))
}
