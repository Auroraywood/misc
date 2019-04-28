## combine active and decoy score files
## Usage: Rscript com
setwd("~/Desktop/tmp/hsp/gnina")
active <- read.table("active_score", stringsAsFactors = F, header = T)
decoy <- read.table("decoy_score", stringsAsFactors = F, header = T)
active$label <- rep(1, length(active[,1]))
decoy$label <- rep(0, length(decoy[,1]))
docked <- rbind(active, decoy)
docked <- docked[order(docked[,2], decreasing = T), ]
head(docked)
write.csv(docked, file = "docked_gnina_hsp90n.csv", row.names = F, quote = F)
