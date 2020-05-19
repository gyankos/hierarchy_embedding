library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(scales)
library(ggrepel)
library(directlabels)
t <- read.csv("final.csv",header=FALSE)
names(t) <- c("Branching","Height", "Algorithm", "Spearman", "Precision", "Recall", "Time")
iterations <- unique(t$Branching)
h <- unique(t$Height)
metrics <- c("Spearman", "Precision", "Recall", "Time")
levels(t$Algorithm) <- revalue(levels(t$Algorithm), c("EuclideanEmbedding"="EE", "Local Density (beta=0.75)"="LD", "Multiple Descent (beta=0.75 alpha=0.5)"="MD1", "Multiple Descent (beta=1 alpha=0.5)"="MD2", "Poincarré k=50"="P50", "Poincarré k=7"="P7", "Relevance Vector"="RV","ParHier Cluster k=3"="PHC","ParHier Sum k=3"="PHS","EEEL k=7"="E7","EEEL k=50"="E50"))
l <- t$Algorithm[!duplicated(t$Algorithm)]

doPlot <- function(m) {
  df <- as.data.frame(unique(t[,c(1,2)]))
  for (x in l) { 
    v <- t[t$Algorithm==x,m]
    length(v) = nrow(df)
    df[x] <- v
  }
  write.csv(df, paste(c(m,"plots.csv"),collapse="_"),row.names = FALSE,quote=FALSE)
  # return (ggplot(Melt1, aes(x = Height, y = get(m), group = Algorithm))  +
  #   geom_line(aes(color=Algorithm,linetype=Algorithm))  +
  #   geom_point(aes(shape=Algorithm)) +
  #   scale_shape_manual(values=1:nlevels(Melt1$Algorithm)) +
  #   ylab(label=m) +
  #   xlab(paste(c("Tree Hight (Branching=",b,")"),collapse="")) +
  #   geom_dl(aes(label = Algorithm), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
  #   geom_dl(aes(label = Algorithm), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8)))
}

for (m in metrics) {
  doPlot(m)
}

m = "Time"
df <- as.data.frame(unique(t[,c(1,2)]))
for (x in c("EE","P50","P7")) { 
  v <- t[t$Algorithm==x,m]
  length(v) = nrow(df)
  df[x] <- v
}
write.csv(df, paste(c(m,"plots.csv"),collapse="_"),row.names = FALSE,quote=FALSE)
