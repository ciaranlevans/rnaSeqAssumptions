require(ggplot2)
require(dplyr)
require(tidyr)


efdrMatPartialAsymmSameMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/efdrMatPartialAsymmSameMRNA_prm2.txt", sep="")


# call it DEGES not TCC
colnames(efdrMatPartialAsymmSameMRNA)[5] <- "DEGES"
colnames(efdrMatPartialAsymmSameMRNA)[4] <- "TMM"

# PS = partial, same 
matToPlotPS <- efdrMatPartialAsymmSameMRNA %>% gather("method", "efdr", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotPS <- matToPlotPS %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=efdr, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=efdr, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Empirical FDR") +
  ggtitle("eFDR by proportion of DE \n asymmetry, same mRNA/cell") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype=2)




efdrMatSymmSameMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/efdrMatSymmSameMRNA_prm2.txt",
           sep="")

# call it DEGES not TCC
colnames(efdrMatSymmSameMRNA)[5] <- "DEGES"
colnames(efdrMatSymmSameMRNA)[4] <- "TMM"

# SS = symm, same
matToPlotSS <- efdrMatSymmSameMRNA %>% gather("method", "efdr", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotSS <- matToPlotSS %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=efdr, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=efdr, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Empirical FDR") +
  ggtitle("eFDR by proportion of DE \n symmetry, same mRNA/cell") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype=2)



efdrMatPartialAsymmDiffMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/efdrMatPartialAsymmDiffMRNA_prm2.txt", sep="")


# call it DEGES not TCC
colnames(efdrMatPartialAsymmDiffMRNA)[5] <- "DEGES"
colnames(efdrMatPartialAsymmDiffMRNA)[4] <- "TMM"

#PD = partial, different
matToPlotPD <- efdrMatPartialAsymmDiffMRNA %>% gather("method", "efdr", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotPD <- matToPlotPD %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=efdr, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=efdr, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Empirical FDR") +
  ggtitle("eFDR by proportion of DE \n asymmetry, different mRNA/cell") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype=2)




efdrMatSymmDiffMRNA <- 
  read.csv("C:/Users/cevans/Desktop/Summer After Grad/efdrMatSymmDiffMRNA_prm2.txt", sep="")

# call it DEGES not TCC
colnames(efdrMatSymmDiffMRNA)[5] <- "DEGES"
colnames(efdrMatSymmDiffMRNA)[4] <- "TMM"

#SD = symm, diff
matToPlotSD <- efdrMatSymmDiffMRNA %>% gather("method", "efdr", 1:6) %>% 
  mutate(propDE = rep(seq(from=0.05, to=0.95, by=0.05), 6))

matToPlotSD <- matToPlotSD %>% group_by(method) %>% ggplot() + 
  geom_point(aes(x=propDE, y=efdr, color=method, shape=method)) +
  geom_line(aes(x=propDE, y=efdr, color=method)) +
  xlab("Proportion of Differential Expression") +
  ylab("Empirical FDR") +
  ggtitle("eFDR by proportion of DE \n symmetry, different mRNA/cell") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        plot.title=element_text(size=16),
        legend.title=element_text(size=14)) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype=2)


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