library(openxlsx)
setwd("~/Mycode/Python")
files <- dir()
active <- read.table(files[1], stringsAsFactors = F, header = T)
decoy <- read.table(files[2], stringsAsFactors = F, header = T)
FDA <- read.table(files[3], stringsAsFactors = F, header = T)
FDA_real <- read.table()
active$label <- rep(1, length(active$ID))
decoy$label <- rep(0, length(decoy$ID))

act_dec <- rbind(active, decoy)
act_dec <- act_dec[order(act_dec$cnnscore, decreasing = T), ]
l_a <- which(act_dec$label==1) # location of active
l <- length(l_a)
sum <- 0 # sum of decoy seen
for(i in 1:l){
  tmp <- act_dec[1:l_a[i], ]
  sum <- sum + length(which(tmp$label==0))
}
auc = 1 - sum/(length(active$ID)*length(decoy$ID))
