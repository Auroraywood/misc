library(stringr)
library(magrittr)
library(dplyr)
library(tidyverse)

dat <- readLines("file:///C:/Users/F/Desktop/ATP subsequence.txt") %>%
  paste(collapse = " ") %>% # implode a vector of strings into a string
  str_split(">", simplify = T)#split by ">"
dat1 <- grep('PROTEIN_KINASE_ATP', dat, value = T) %>% # retain vectors have "PROTEIN_KINASE_ATP"
  trimws(which = "both") # remove leading and/or trailing whitespace from character strings
dat2 <- data.frame(v1=dat1) # separate one column into multiple columns
dat3 <- separate(dat2, col = "v1", into = c("pro_id", "ATP_bind_pocket"),  sep = " ")
bind_sq <- dat3

extract_pro_id <- function(x){
  t1 <- x[1]
  t2 <- str_extract(t1, "\\|......\\|")
  t3 <- sub("^\\|", "", t2)
  t4 <- sub("\\|$", "", t3)
  return(t4)
}

bind_sq$pro_id <- apply(dat3, 1, extract_pro_id) # extract protein id
mkl_pro_id <- read.table("file:///E:/desktop/work/MKL/pairwiseMKL-master/DTI_data/Uniprot_IDs_226kin.txt", stringsAsFactors = F, header = F)
colnames(mkl_pro_id) <- "pro_id"
res <- plyr::join(mkl_pro_id, bind_sq, by = "pro_id", type = "left")
res$ATP_bind_pocket <- toupper(res$ATP_bind_pocket)

keep_letter <- function(x){ # function of retain capital letter 
  t1 <- x[2]
  t2 <- unlist(str_extract_all(t1, "[A-Z]"))
  t3 <- paste(t2, collapse = "")
  return(t3)
}
keep_letter(res[3,])

res$ATP_bind_pocket <- apply(res, 1, keep_letter)
test <- res[1:5, 2, drop=F]
write.csv(test, file = "test1.csv")
