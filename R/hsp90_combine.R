## 
## @Aim: Study how to impove auc by combine two score methods.
## @Author: zdx
## @dataset: hsp90

library(stringr)
## normalize socre
normalize <- function(x){ # score must be in 2nd column
  t <- x
  t <- (t - mean(t))/sd(t) %>% # compute zscore
    tanh()
  return(t)
}
## function of compute auc
auc_compute <- function(df){
  l_a <- which(df$label==1) # location of active
  l <- length(l_a)
  sum <- 0 # sum of decoy seen
  for(i in 1:l){
    tmp <- df[1:l_a[i], ]
    sum <- sum + length(which(tmp$label==0))
  }
  auc = 1 - sum/(length(which(df$label==1))*length(which(df$label==0)))
  return(auc)
}
## deepdock
hsp <- read.csv("C:/Users/F/Desktop/score/hsp/deepdock/hsp90_0424_all_dude.csv", header = T, stringsAsFactors = F)
colnames(hsp) <- c("ID", "label", "deepdock_score")
hsp <- hsp[order(hsp$deepdock_score, decreasing = T), ] # rank by score
extract_ID <- function(x){ # function of extract ligand name
  return(sub("_pose_00", "", x[1])) 
}
hsp$ID <- apply(hsp, 1, extract_ID)

## gnina
## active - decoy
active <- read.table("C:/Users/F/Desktop/score/hsp/gnina/active_score", stringsAsFactors = F, header = T)
decoy <- read.table("C:/Users/F/Desktop/score/hsp/gnina/decoy_score", stringsAsFactors = F, header = T)
active$label <- rep(1, length(active[,1]))
decoy$label <- rep(0, length(decoy[,1]))
act_dec <- rbind(active, decoy)
act_dec <- act_dec[order(act_dec[,2], decreasing = T), ]
t <- c('CHEMBL1213898',
       'CHEMBL1213899',
       'CHEMBL190182',
       'CHEMBL191074',
       'CHEMBL200102',
       'CHEMBL252164',
       'CHEMBL254875',
       'CHEMBL278315',
       'CHEMBL361078',
       'CHEMBL362307',
       'CHEMBL364968',
       'CHEMBL365159',
       'CHEMBL365617',
       'CHEMBL365700',
       'CHEMBL366215',
       'CHEMBL371380',
       'CHEMBL399530',
       'CHEMBL402712',
       'CHEMBL404630',
       'CHEMBL424811',
       'CHEMBL426446',
       'CHEMBL467399')
l <- length(t)
for(i in 1:l){  # remove repeat ones
  need_del <- t[i]
  row <- which(act_dec$ID==need_del)
  act_dec <- act_dec[-row,]
}
rm(active, decoy, i, l, need_del, row, t)
head(act_dec)
head(hsp)
## merge two score
merge_score <- merge(hsp[, c("ID", "deepdock_score")], act_dec, by.x = "ID")
# =========================== data has been prepared =================================

main <- function(){
  tmp <- merge_score
  tmp$com_score <- tmp[,2] * tmp[,3] #---------------------- one 0.8233
  tmp$com_score <- 1 - (1-tmp[,2]) * (1-tmp[,3]) # --------- two 0.7765
  tmp$com_score <- (normalize(tmp[,2]) + normalize(tmp[,3]))/2
  tmp <- tmp[order(tmp$com_score, decreasing = T), ]
  head(tmp)
  auc_compute(tmp)
}
