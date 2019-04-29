# compute auc of single file
library(PRROC) 
# fg: vector of postive prob; bg: vector of negative prob;
# roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)  
# roc$auc
# pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
# pr$auc.davis.goadrich

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
dat <- read.csv("/home/zdx/Desktop/score/hsp90/docked_gnina_hsp90n.csv", header = T, stringsAsFactors = F)
auc(dat$label)
aupr(dat$label)
active <- read.table("/home/zdx/Desktop/scoring_smina/active.smina.tsv", header = T, stringsAsFactors = F)
active$label <- 1
decoy <- read.table("/home/zdx/Desktop/scoring_smina/decoy.smina.tsv", header = T, stringsAsFactors = F)
decoy$label <- 0
dat <- rbind(active, decoy)
colnames(dat)
dat <- dat[order(dat$Best_energy, decreasing = F), ]
head(dat)

probs <- dat$cnnscore
fg <- probs[dat$label == 1]
bg <- probs[dat$label == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
roc$auc
plot(roc)
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr$auc.davis.goadrich
plot(pr)
  
  
  
  
  
  
  
  
  
  
  
  






