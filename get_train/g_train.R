# label receptor ligand
setwd("/home/zdx/cnndock/get_train")
active <- read.table("active_files.txt", header = F, stringsAsFactors = F)
decoy <- read.table("decoy_files.txt", header = F, stringsAsFactors = F)
receptor <- read.table("receptor.txt", header = F, stringsAsFactors = F)
receptor <- receptor[1,1]

train <- rbind(active, decoy)
l <- length(train$V1)
random_id <- sample(1:l, l) # shuffle
train <- train[random_id, ]

receptors <- rep(receptor, l)
result <- data.frame(label=train$V1, receptor=receptors, ligand=train$V2)

write.table(result, file = "trainfile.type", col.names = FALSE, row.names = FALSE, quote = F)

