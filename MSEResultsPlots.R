require(ggplot2)
require(dplyr)
require(tidyr)


mseMatPartialAsymmSameMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/MSEMatPartialAsymmSameMRNA.txt", sep="")


# call it DEGES not TCC
colnames(mseMatPartialAsymmSameMRNA)[5] <- "DEGES"
colnames(mseMatPartialAsymmSameMRNA)[4] <- "TMM"

# PS = partial, same 
matToPlotPS <- mseMatPartialAsymmSameMRNA %>% gather("method", "MSE", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotPS <- matToPlotPS %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=MSE, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=MSE, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Fold Change MSE") +
  ggtitle("Fold change MSE by proportion of DE \n asymmetry, same mRNA/cell") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw()



mseMatSymmSameMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/MSEMatSymmSameMRNA.txt",
           sep="")

# call it DEGES not TCC
colnames(mseMatSymmSameMRNA)[5] <- "DEGES"
colnames(mseMatSymmSameMRNA)[4] <- "TMM"

# SS = symm, same
matToPlotSS <- mseMatSymmSameMRNA %>% gather("method", "MSE", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotSS <- matToPlotSS %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=MSE, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=MSE, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Fold Change MSE") +
  ggtitle("Fold change MSE by proportion of DE \n symmetry, same mRNA/cell") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw()




mseMatPartialAsymmDiffMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/MSEMatPartialAsymmDiffMRNA.txt", sep="")


# call it DEGES not TCC
colnames(mseMatPartialAsymmDiffMRNA)[5] <- "DEGES"
colnames(mseMatPartialAsymmDiffMRNA)[4] <- "TMM"

#PD = partial, different
matToPlotPD <- mseMatPartialAsymmDiffMRNA %>% gather("method", "MSE", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotPD <- matToPlotPD %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=MSE, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=MSE, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Fold Change MSE") +
  ggtitle("Fold change MSE by proportion of DE \n asymmetry, different mRNA/cell") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw()




mseMatSymmDiffMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/MSEMatSymmDiffMRNA.txt", sep="")

# call it DEGES not TCC
colnames(mseMatSymmDiffMRNA)[5] <- "DEGES"
colnames(mseMatSymmDiffMRNA)[4] <- "TMM"

#SD = symm, diff
matToPlotSD <- mseMatSymmDiffMRNA %>% gather("method", "MSE", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotSD <- matToPlotSD %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=MSE, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=MSE, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Fold Change MSE") +
  ggtitle("Fold change MSE by proportion of DE \n symmetry, different mRNA/cell") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw()


multiplot(matToPlotPS, matToPlotSS, matToPlotPD, matToPlotSD, cols=2)



# this function borrowed from Cookbook for R
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}