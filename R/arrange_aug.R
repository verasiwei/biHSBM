#' Arrange aggregation matrix and apply different methods
#'
#' \code{arrange_aug} returns a data frame
#'
#' @param N The number of aggregation times
#' @param rep The number of repeats for each aggregation time
#' @param dat1 The list of one method from augmentation_matrix function
#' @param dat2 The list of one method from augmentation_matrix function
#' @return
#' @examples
#' dat <- arrange_aug(50,10,aug_result1,aug_result2)

arrange_aug <- function(N,rep,dat1,dat2){
  for (i in 1:N) {
    for (j in 1:rep) {
      if(sum(dat1[[i]][[1]][[j]]==1)>sum(dat1[[i]][[1]][[j]]==2)){
        dat1[[i]][[1]][[j]]=ifelse(dat1[[i]][[1]][[j]]==1,2,1)
      }
      if(sum(dat2[[i]][[1]][[j]]==1)>sum(dat2[[i]][[1]][[j]]==2)){
        dat2[[i]][[1]][[j]]=ifelse(dat2[[i]][[1]][[j]]==1,2,1)
      }
      if(sum(dat1[[i]][[2]][[j]]==1)>sum(dat1[[i]][[2]][[j]]==2)){
        dat1[[i]][[2]][[j]]=ifelse(dat1[[i]][[2]][[j]]==1,2,1)
      }
      if(sum(dat2[[i]][[2]][[j]]==1)>sum(dat2[[i]][[2]][[j]]==2)){
        dat2[[i]][[2]][[j]]=ifelse(dat2[[i]][[2]][[j]]==1,2,1)
      }
    }
  }

  nmi_val <- list()
  nmi <- c()
  mean <- c()
  min <- c()
  max <- c()
  se <- c()
  for(i in 1:N){
    for (j in 1:rep) {
      nmi <- c(nmi,aricode::NMI(dat1[[i]][[1]][[j]],dat2[[i]][[2]][[j]]))
    }
    nmi_val[[i]] <- nmi
    min <- c(min,min(nmi))
    max <- c(max,max(nmi))
    mean <- c(mean,mean(nmi))
    se <- c(se,sqrt(var(nmi)/length(nmi)))
    nmi <- c()
  }
 return(dat <- data.frame(y=mean,min=min,max=max,se=se))
}
