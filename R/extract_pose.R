#
#USE: Rscript extract_pose.R <data root dir>
library(stringr)

args=commandArgs(T)
setwd(args[1])
file.remove("tmp")

#files <- dir()
#f_l <- length(files) # number of files

# if(!dir.exists("tmp")){
#   dir.create("tmp")
# }
# 
# for(i in 1:f_l){
#   i <- 1
#   f_n <- files[i]
#   dat <- readLines(f_n)
#   tmp1 <- as.character(str_match(dat, "ENDMDL"))
#   tmp2 <- which(tmp1=="ENDMDL")
#   f_t <- unlist(str_split(f_n, "\\.")) # split
#   path <- paste0("./tmp/", f_t[1], "new", ".", f_t[2])
#   dat_out <- dat[2:(tmp2[1]-1)]
#   write.table(dat_out, file = path, quote = F, row.names = F, col.names = F)
# }
