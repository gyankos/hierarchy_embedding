library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library("scales")
library(ggrepel)
setwd("/media/giacomo/Data/hierarchy_paper/projects/hierarchy_tests")
t <- read.csv("final_experiments.csv",header=FALSE)
names(t) <- c("Branching","Height", "Algorithm", "Spearman", "Precision@<=k", "Recall@>k", "MultiThreaded Time")
iterations <- unique(t$Branching)
h <- unique(t$Height)
metrics <- c("Spearman", "Precision@<=k", "Recall@>k", "MultiThreaded Time")
levels(t$Algorithm) <- revalue(levels(t$Algorithm), c("EuclideanEmbedding"="EE", "Local Density (beta=0.75)"="LD", "Multiple Descent (beta=0.75 alpha=0.5)"="MD1", "Multiple Descent (beta=1 alpha=0.5)"="MD2", "Poincarré k=50"="P50", "Poincarré k=7"="P7", "Relevance Vector"="RV","ParHier Cluster k=3"="PHC","ParHier Sum k=3"="PHS","EEEL k=7"="E7","EEEL k=50"="E50"))

for (b in iterations) {
  for (m in metrics) {
    Plot1 <- t[t$Branching == b,]
    Melt1 <- Plot1[,c("Height", "Algorithm", m)]
    df <- data.frame(Height=h)
    l <- 
      group_split(Melt1 %>%
                    group_by(Algorithm))
    for (x in l) { 
      df[paste(as.character(unique(x$Algorithm)), sep="")] <- as.vector(x[,c(m)])
    }
    write.csv(df, paste(c(b,m,"plots.csv"),collapse="_"),row.names = FALSE)
    # # pdf( paste(c(b,m,"plots.pdf"),collapse="_"))
    # ggplot(Melt1, aes(x = Height, y = get(m), group = Algorithm))  +
    #   geom_line(aes(color=Algorithm))  + 
    #   geom_point(aes(shape=Algorithm)) +
    #   scale_shape_manual(values=1:nlevels(Melt1$Algorithm)) +
    #   ylab(label=m) +
    #   xlab(paste(c("Tree Hight (Branching=",b,")"),collapse="")) +
    #   geom_dl(aes(label = Algorithm), method = list(dl.trans(x = x + 0.4), "last.points", cex = 0.8)) +
    #   geom_dl(aes(label = Algorithm), method = list(dl.trans(x = x - 0.4), "first.points", cex = 0.8))
    #   #geom_label_repel(aes(label = Algorithm),
    #   #                 nudge_x = 1,
    #   #                 na.rm = TRUE)
    #   geom_dl(aes(label = Algorithm), method = list(dl.combine("first.points", "last.points"), cex = .6))
    # # dev.off()
  }
}
