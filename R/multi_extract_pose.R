# Usage: Rscript multi_extract_pose.R ~/Documents/asfvtk_virtual_screening
library(stringr)

args <- commandArgs(T)

# data_root <- "~/Documents/asfvtk_virtual_screening"
data_root <- args[1]
setwd(data_root)
folders <- dir() # datasets need to be processed
folders_path <- paste0(data_root,"/", folders, "/docking_pdbqt") # path of datasets need to be processed

count <- 1
for(set in folders_path){
  # set <- folders_path[1]
  setwd(set)
  folders <- dir()
  if(count==1){
    sub_sets <- paste0(set, "/", folders)
    count <- count + 1
  }else{
    tmp <- paste0(set, "/", folders)
    sub_sets <- c(sub_sets, tmp)
    count <- count + 1
  }
}

for(item in sub_sets){
  setwd(item)
  files <- dir() # set working dir
  if(dir.exists("tmp")){# if "tmp" folder has been created, "tmp" needed to be deleted from files
    system("rm -r tmp")
    files <- setdiff(files, "tmp")
    dir.create("tmp")
  }else{
    dir.create("tmp")
  }
  
  f_l <- length(files) # number of files
  
  for(i in 1:f_l){
    # i = 1
    tryCatch({
      cat(i, "\n")
      f_n <- files[i]
      dat <- readLines(f_n)
      tmp1 <- as.character(str_match(dat, "ENDMDL"))
      tmp2 <- which(tmp1=="ENDMDL")
      dat_out <- dat[2:(tmp2[1]-1)]
      #f_t <- unlist(str_split(f_n, "\\.")) # split
      #path <- paste0("./tmp/", f_t[1], ".", f_t[2])
      path <- paste0("./tmp/", f_n)
      write.table(dat_out, file = path, quote = F, row.names = F, col.names = F)
    }, error=function(e){
      cat("ERROR :",f_n, "maybe has problem.", "\n")
    }
    )
  }
}
cat("\nSuccess!\n\n")
