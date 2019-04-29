library(openxlsx)

work_dir <- "C:/Users/F/Desktop/"
active_path <- "active_score"
decoy_path <- "decoy_score"
FDA_path <- "FDA_Approved.smina.tsv"
Merck_path <- "Merck_Kinase.smina.tsv"

setwd(work_dir)
## active - decoy
active <- read.table(active_path, stringsAsFactors = F, header = T)
decoy <- read.table(decoy_path, stringsAsFactors = F, header = T)
active$label <- rep(1, length(active[,1]))
decoy$label <- rep(0, length(decoy[,1]))
act_dec <- rbind(active, decoy)
act_dec <- act_dec[order(act_dec[,2], decreasing = T), ]
#write.csv(act_dec, file = "ours_score.csv", row.names = F)
## FDAFDA_Approved
FDA_score <- read.table(FDA_path, stringsAsFactors = F, header = T)
colnames(FDA_score) <- c("ID","score")
FDA_score <- FDA_score[order(FDA_score$score, decreasing = T), ]
FDA_real <- read.xlsx("topoII.true_label.xlsx", sheet = 1)
FDA_real <- FDA_real[,c("Cmpd_ID", "Active")]
colnames(FDA_real) <- c("ID","label")
FDA_final <- plyr::join(FDA_score, FDA_real, by = "ID", type = "left")
#write.csv(FDA_final, file = "FDA_Approved_score.csv", row.names = F)
## Merck_Kinase
Merck_Kinase_score <- read.table(Merck_path, stringsAsFactors = F, header = T)
colnames(Merck_Kinase_score) <- c("ID", "score")
Merck_Kinase_score <- Merck_Kinase_score[order(Merck_Kinase_score$score, decreasing = T), ]
Merck_Kinase_real <- read.xlsx("topoII.true_label.xlsx", sheet = 2)
Merck_Kinase_real <- Merck_Kinase_real[,c("Cmpd_ID", "Active")]
colnames(Merck_Kinase_real) <- c("ID","label")
Merck_Kinase_final <- plyr::join(Merck_Kinase_score, Merck_Kinase_real, by = "ID", type = "left")
#write.csv(Merck_Kinase_final, file = "Merck_Kinase_score.csv", row.names = F)

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

decoy_active <- auc_compute(act_dec)

Merck_auc <- auc_compute(Merck_Kinase_final)

FDA_auc <- auc_compute(FDA_final)

cat("FDA_auc: ", FDA_auc, "\n", "decoy_active: ", decoy_active,"\n", "Merck_auc: ", Merck_auc)
