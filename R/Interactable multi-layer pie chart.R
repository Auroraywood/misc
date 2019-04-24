#Interactable multi-layer pie chart
dat <- read.csv("file:///C:/Users/F/Desktop/2019-03-26-test.csv", stringsAsFactors = F)
m_receptor <- unique(dat$l1) #elemant in level 1

result <- NULL
for(i in 1:8){
  if(i==1){
    tmp_receptor <- m_receptor[i] #current elemant of level 1
    tmp_df <- dat[which(dat$l1==tmp_receptor), ] #dataframe of current elemant
    m <- data.frame(table(tmp_df$l2))
    t_l <- length(m[,1])
    m$l1 <- rep(tmp_receptor, t_l)
    result <- m
  }
  else{
    tmp_receptor <- m_receptor[i] #current elemant of level 1
    tmp_df <- dat[which(dat$l1==tmp_receptor), ] #dataframe of current elemant
    m <- data.frame(table(tmp_df$l2))
    t_l <- length(m[,1])
    m$l1 <- rep(tmp_receptor, t_l)
    result <- rbind(result, m)
  }
}
res <- data.frame(l1=result$l1, l2=result$Var1, freq=result$Freq)

past_fun <- function(line){
  left <- line[1]
  right <- line[2]
  return(paste0(left, "-", right))
}
sequences <- data.frame(V1=apply(res, 1, past_fun), V2=res$freq)
