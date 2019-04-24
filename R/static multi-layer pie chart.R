# Membrane receptor: ME
# Enzyme: EN
# Other: OT
# Transcription factor: TR
# Auxiliary transport protein: AU
# Ion channel: IO
# Epigenetic regulator: EP
# Other cytosolic protein: OC

dat <- read.csv("file:///C:/Users/F/Desktop/Tclin100_protein_family_annotation.csv", na.strings = "", stringsAsFactors = F)
dat <- read.csv("file:///C:/Users/F/Desktop/kinase_326_combine_283_287_protein_family_annotation.csv", na.strings = "", stringsAsFactors = F)
dat <- read.csv("file:///C:/Users/F/Desktop/DUDE_protein_family_annotation.csv", na.strings = "", stringsAsFactors = F)
dat[is.na(dat)] <- "other"
m_receptor <- unique(dat$l1) #elemant in level 1
l <- length(m_receptor)
result <- NULL
for(i in 1:l){
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

res$total <- with(res, ave(freq, l1, FUN = sum))
level1 <- res[,c(1, 4)]
library(tidyverse)
level1_unique <- level1[!duplicated(level1$l1), ]
#write.csv(res, file = "res.csv", row.names = F)

par(mai=c(1,1,1,2.2))#bottom, left, top, and right
library(colorspace)
percentage <- res$freq/sum(res$freq)
percentage <- (round(percentage, digits = 4))*100
lab <- paste0(res$l2, ":", res$freq,"(", percentage, "%",")")
pie_zdx <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                     init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                     col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.15) * P$x, c(1, 1.15) * P$y)
      text(1.19 * P$x, 1.19 * P$y, labels[i], xpd = TRUE, adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}
col1 <- c("#9FA2FF", "#F1F1F1", "#BAAE00")
pie_zdx(res$freq, border = NA, radius = 1, col = col1, labels = lab, cex = 0.5) #first circle
par(new=TRUE)
pie(1, border = NA, radius = 0.85, col = "white",labels = NA)
pie_zdx1 <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                     init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                     col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.7) * P$x, c(1, 1.7) * P$y)
      text(1.8 * P$x, 1.8 * P$y, labels[i], xpd = TRUE, adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}
percentage <- level1_unique$total/sum(level1_unique$total)
percentage <- (round(percentage, digits = 2))*100
lab1 <- paste0(level1_unique$l1, ":", level1_unique$total,"(", percentage, "%",")")
#par1 <- choose_palette() 
par(new=TRUE)
pie_zdx1(level1_unique$total, border = NA, radius = 0.8, col = rainbow(length(level1_unique$total)), labels = NA, cex = 0.5)
legend("topleft", title = "Level 1", legend = lab1, bty = "n",cex = 0.5, pch = 15, col = rainbow(length(level1_unique$total)))
par(new=TRUE)
pie(1, border = NA, radius = 0.65, col = "white",labels = NA)
text(locator(8), labels = lab1, col = "blue", cex = 0.5)

## loop creates pie slices
plot.new() 
par(omi = c(0.5,0.5,0.75,0.5), mai = c(0.1,0.1,0.1,0.1), las = 1)
