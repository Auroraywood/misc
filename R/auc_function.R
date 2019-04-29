## function of compute auc
auc <- function(x){
  l_a <- which(x==1)
  l <- length(l_a)
  sum <- 0
  for(i in 1:l){
    tmp <- x[1:l_a[i]]
    sum <- sum + length(which(tmp==0))
  }
  a <- 1 - sum/(length(which(x==1))*length(which(x==0)))
  return(a)
}