#
#USE: Rscript extract_pose.R <data root dir>
library(stringr)

args=commandArgs(T)
setwd(args[1])
files <- dir() # set working dir
if(dir.exists("tmp")){# if "tmp" folder has been created, "tmp" needed to be deleted from files
  files <- setdiff(files, "tmp")
}else{
  dir.create("tmp")
}

f_l <- length(files) # number of files

for(i in 1:f_l){
  tryCatch({
    cat(i, "\n")
    f_n <- files[i]
    dat <- readLines(f_n)
    tmp1 <- as.character(str_match(dat, "ENDMDL"))
    tmp2 <- which(tmp1=="ENDMDL")
    dat_out <- dat[2:(tmp2[1]-1)]
    f_t <- unlist(str_split(f_n, "\\.")) # split
    path <- paste0("./tmp/", f_t[1], "new", ".", f_t[2])
    write.table(dat_out, file = path, quote = F, row.names = F, col.names = F)
  }, error=function(e){
      cat("ERROR :",f_n, "maybe has problem.", "\n")
    }
  )
}
