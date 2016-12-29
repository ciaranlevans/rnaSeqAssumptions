# The code from this simulation borrows substantially from
# the code used by Law et al. (2014)

# packages needed for normalization and testing
library(DESeq2) # for the DE hypothesis testing
library(psych) # for calculating geometric mean
library(PoissonSeq) # for PoissonSeq normalization
library(truncnorm) # for generating sequencing depths
library(edgeR) # for TMM normalization
library(TCC) # for DEGES normalization
require(dplyr) # for data manipulation
require(tidyr) # for data manipulation

# Use equal or unequal library sizes
equal <- TRUE

# Use inverse chi-square or log-normal dispersion
invChisq <- TRUE

# will differential expression be symmetric?
symmetric <- TRUE

# if asymmetric, how asymmetric?
partiallyAsymm <- FALSE
completelyAsymm <- FALSE

# Get distribution function of abundance proportions
# This distribution was generated from a real dataset
load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))

# cutoff for BH  adjust
cutoff <- 0.05

# Generate baseline proportions for desired number of genes
ngenes <- 10000
baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
baselineprop <- baselineprop/sum(baselineprop)

n1 <- 2
n2 <- 2

nlibs <- n1+n2

# Library size 
if(equal){
  expected.lib.size <- rep(11e6,nlibs)
} else {
  expected.lib.size <- 20e6 * c(1,0.1,1,0.1,1,0.1)
}



# number of runs
nsamp <- 50

# percentages of differential expression
percDE <- seq(from=0.05, to=0.95, by=0.05)


# matrix to hold the results (mean eFDR)
efdrMat <- matrix(nrow=19, ncol=6)

# matrix to hold the standard devs of the nsamp efdr measures
efdrSDMat <- matrix(nrow=19, ncol=6)

# matrix to hold the power results (mean power)
powerMat <- matrix(nrow=19, ncol=6)

# matrix to hold the std devs of the power 
powerSDMat <- matrix(nrow=19, ncol=6)

# matrix to hold the average MSE for each method
MSEMat <- matrix(nrow=19, ncol=6)

set.seed(21)


for(k in percDE){
  
  ndiffExpress <- ngenes*k
  
  # store FDR values
  efdrOracle <- c()
  efdrTotCount <- c()
  efdrDESeq2 <- c()
  efdrEdgeR <- c()
  efdrTCC <- c()
  efdrPSeq <- c()
  
  # store power values
  powerOracle <- c()
  powerTotCount <- c()
  powerDESeq2 <- c()
  powerEdgeR <- c()
  powerTCC <- c()
  powerPSeq <- c()
  
  # store MSE values
  mseOracle <- c()
  mseTotCount <- c()
  mseDESeq2 <- c()
  mseEdgeR <- c()
  mseTCC <- c()
  msePSeq <- c()
  
  for(j in 1:nsamp){
    i <- sample(1:ngenes, ndiffExpress)
    
    # baseline proportions for each gene; DE genes get multiplied by the fold change
    baselineprop1 <- baselineprop2 <- baselineprop
    
    # sequencing depths
    sizeRatio <- rtruncnorm(nlibs, a = 0.6, b = 1.4, mean=1, sd=0.2)
    
    # we can either have asymmetric differential expression, with unequal proportions
    # of up- and down-regulated genes, or symmetric differential expression
    if (symmetric){
      i1 <- i[1:(ndiffExpress/2)]
      i2 <- i[(ndiffExpress/2 + 1):ndiffExpress]
      
      baselineprop1[i1] <- baselineprop1[i1]*4
      baselineprop2[i2] <- baselineprop2[i2]*6
    } else if (completelyAsymm){
      baselineprop1[i] <- baselineprop1[i]*fc 
    } else {
      i1 <- i[1:floor(3*ndiffExpress/4)]
      i2 <- i[(floor(3*ndiffExpress/4) + 1):ndiffExpress]
      
      baselineprop1[i1] <- baselineprop1[i1]*fc
      baselineprop2[i2] <- baselineprop2[i2]*fc
    }
    
    # make these actual proportions
    baseSum1 <- sum(baselineprop1)
    baseSum2 <- sum(baselineprop2)
    baselineprop1 <- baselineprop1/sum(baselineprop1)
    baselineprop2 <- baselineprop2/sum(baselineprop2)
    
    # mean expression for each gene
    mu0.1 <- matrix(baselineprop1,ngenes,1) %*% 
      matrix(expected.lib.size[1:n1]*sizeRatio[1:n1],1,n1)
    mu0.2 <- matrix(baselineprop2,ngenes,1) %*% 
      matrix(expected.lib.size[(n1+1):(n1+n2)]*sizeRatio[(n1+1):(n1+n2)],1,n2)
    
    mu0 <- cbind(mu0.1,mu0.2)
    
    # status will keep track of which genes are differentially expressed
    status <- rep(0,ngenes)
    status[i] <- 1
    
    # Biological variation
    BCV0 <- 0.2+1/sqrt(mu0)
    if(invChisq){
      df.BCV <- 40
      BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
    } else {
      BCV <- BCV0*exp( rnorm(ngenes,mean=0,sd=0.25)/2 )
    }
    if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
    shape <- 1/BCV^2
    scale <- mu0/shape
    mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)
    
    # Technical variation
    counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)
    
    #Filter: following limma/voom code, remove genes with
    # rowsums < 10 in the read count matrix
    keep <- rowSums(counts)>=10
    nkeep <- sum(keep)
    counts2 <- counts[keep,]
    status2 <- status[keep]
    
    
    
    ## naive-DESeq2
    cond <- factor(rep(1:2, each=n1))
    dds <- DESeqDataSetFromMatrix(counts2, DataFrame(cond), ~ cond)
    dds <- estimateSizeFactors(dds)
    sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    
    res$padj <- p.adjust(res$pvalue, method="BH")
    efdrDESeq2[j] <- 
      sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)
    
    powerDESeq2[j] <- 
      sum(res$padj[status2==1] < cutoff)/sum(status2==1)
    
    
    # want non-DE genes with no zero counts
    countsToCompare <- as.data.frame(cbind(counts2[,1][status2==0], 
                                           counts2[,nlibs][status2==0]))
    colnames(countsToCompare) <- c("col1", "col2")
    countsToCompare <- filter(countsToCompare, col1!=0 & col2!=0)
    
    
    
    normCounts1 <- countsToCompare[,1]/(sizeFactors(dds)[1])
    normCounts2 <- countsToCompare[,2]/(sizeFactors(dds)[nlibs])
    mseDESeq2[j] <- mean(log(normCounts1/normCounts2)^2)
    
    
    
    # Oracle DESEq2
    normFactors <- sizeRatio/rep(c(baseSum1, baseSum2), each=n1)
    sizeFactors(dds) <- normFactors/geometric.mean(normFactors)
    
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    
    res$padj <- p.adjust(res$pvalue, method="BH")
    efdrOracle[j] <- 
      sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)
    
    powerOracle[j] <- 
      sum(res$padj[status2==1] < cutoff)/sum(status2==1)
    
    normCounts1 <- countsToCompare[,1]/(sizeFactors(dds)[1])
    normCounts2 <- countsToCompare[,2]/(sizeFactors(dds)[nlibs])
    mseOracle[j] <- mean(log(normCounts1/normCounts2)^2)
    
    
    
    
    # Total count normalization
    # rather than divide the counts by the respective total count, we
    # make size factor estimates to plug into DESeq2 by dividing each total
    # count by the geometric mean of the total counts
    totCounts <- colSums(counts2)
    sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    
    res$padj <- p.adjust(res$pvalue, method="BH")
    efdrTotCount[j] <- 
      sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)
    
    powerTotCount[j] <- 
      sum(res$padj[status2==1] < cutoff)/sum(status2==1)
    
    normCounts1 <- countsToCompare[,1]/(sizeFactors(dds)[1])
    normCounts2 <- countsToCompare[,2]/(sizeFactors(dds)[nlibs])
    mseTotCount[j] <- mean(log(normCounts1/normCounts2)^2)
    
    
    
    
    # edgeR normalization
    normFactors <- edgeR::calcNormFactors(counts2)
    normFactors <- normFactors*totCounts
    normFactors <- normFactors/geometric.mean(normFactors)
    sizeFactors(dds) <- normFactors
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    
    res$padj <- p.adjust(res$pvalue, method="BH")
    efdrEdgeR[j] <- 
      sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)
    
    powerEdgeR[j] <- 
      sum(res$padj[status2==1] < cutoff)/sum(status2==1)
    
    normCounts1 <- countsToCompare[,1]/(sizeFactors(dds)[1])
    normCounts2 <- countsToCompare[,2]/(sizeFactors(dds)[nlibs])
    mseEdgeR[j] <- mean(log(normCounts1/normCounts2)^2)
    
    
    
    ## TCC estimate
    ## TCC may fail when proportion of DE is too high
    ## so 0.95 and 1 are not simulated
    if(k < 0.95){
      tcc <- new("TCC", counts2, cond)
      tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                                  iteration = 3, FDR = 0.1, floorPDEG = 0.05)
      tccEsts <- tcc$norm.factors
      tccEsts <- tccEsts*colSums(counts2)/geometric.mean(colSums(counts2))
      
      sizeFactors(dds) <- tccEsts/geometric.mean(tccEsts)
      
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds)
      res <- results(dds)
      
      res$padj <- p.adjust(res$pvalue, method="BH")
      efdrTCC[j] <- 
        sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)
      
      powerTCC[j] <- 
        sum(res$padj[status2==1] < cutoff)/sum(status2==1)
      
      normCounts1 <- countsToCompare[,1]/(sizeFactors(dds)[1])
      normCounts2 <- countsToCompare[,2]/(sizeFactors(dds)[nlibs])
      mseTCC[j] <- mean(log(normCounts1/normCounts2)^2)
    }
    
    
    
    
    ## PoissonSeq estimate
    psests <- PS.Est.Depth(counts2)
    
    sizeFactors(dds) <- psests/geometric.mean(psests)
    
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    
    res$padj <- p.adjust(res$pvalue, method="BH")
    efdrPSeq[j] <- 
      sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)
    
    powerPSeq[j] <- 
      sum(res$padj[status2==1] < cutoff)/sum(status2==1)
    
    normCounts1 <- countsToCompare[,1]/(sizeFactors(dds)[1])
    normCounts2 <- countsToCompare[,2]/(sizeFactors(dds)[nlibs])
    msePSeq[j] <- mean(log(normCounts1/normCounts2)^2)
    
    
    print(paste(k,j)) 
  }
  efdrMat[(k*20),1] <- mean(efdrOracle, na.rm=T)
  efdrMat[(k*20),2] <- mean(efdrTotCount, na.rm=T)
  efdrMat[(k*20),3] <- mean(efdrDESeq2, na.rm=T)
  efdrMat[(k*20),4] <- mean(efdrEdgeR, na.rm=T)
  if(k < 0.95){
    efdrMat[(k*20),5] <- mean(efdrTCC, na.rm=T)
  } else{
    efdrMat[(k*20),5] <- NA
  }
  efdrMat[(k*20),6] <- mean(efdrPSeq, na.rm=T)
  
  efdrSDMat[(k*20),1] <- sd(efdrOracle, na.rm=T)
  efdrSDMat[(k*20),2] <- sd(efdrTotCount, na.rm=T)
  efdrSDMat[(k*20),3] <- sd(efdrDESeq2, na.rm=T)
  efdrSDMat[(k*20),4] <- sd(efdrEdgeR, na.rm=T)
  if(k < 0.95){
    efdrSDMat[(k*20),5] <- sd(efdrTCC, na.rm=T)
  } else{
    efdrSDMat[(k*20),5] <- NA
  }
  efdrSDMat[(k*20),6] <- sd(efdrPSeq, na.rm=T)
  
  powerMat[(k*20),1] <- mean(powerOracle, na.rm=T)
  powerMat[(k*20),2] <- mean(powerTotCount, na.rm=T)
  powerMat[(k*20),3] <- mean(powerDESeq2, na.rm=T)
  powerMat[(k*20),4] <- mean(powerEdgeR, na.rm=T)
  if(k < 0.95){
    powerMat[(k*20),5] <- mean(powerTCC, na.rm=T)
  } else{
    powerMat[(k*20),5] <- NA
  }
  powerMat[(k*20),6] <- mean(powerPSeq, na.rm=T)
  
  powerSDMat[(k*20),1] <- sd(powerOracle, na.rm=T)
  powerSDMat[(k*20),2] <- sd(powerTotCount, na.rm=T)
  powerSDMat[(k*20),3] <- sd(powerDESeq2, na.rm=T)
  powerSDMat[(k*20),4] <- sd(powerEdgeR, na.rm=T)
  if(k < 0.95){
    powerSDMat[(k*20),5] <- sd(powerTCC, na.rm=T)
  } else{
    powerSDMat[(k*20),5] <- NA
  }
  powerSDMat[(k*20),6] <- sd(powerPSeq, na.rm=T)
  
  
  MSEMat[(k*20),1] <- mean(mseOracle, na.rm=T)
  MSEMat[(k*20),2] <- mean(mseTotCount, na.rm=T)
  MSEMat[(k*20),3] <- mean(mseDESeq2, na.rm=T)
  MSEMat[(k*20),4] <- mean(mseEdgeR, na.rm=T)
  if(k < 0.95){
    MSEMat[(k*20),5] <- mean(mseTCC, na.rm=T)
  } else{
    MSEMat[(k*20),5] <- NA
  }
  MSEMat[(k*20),6] <- mean(msePSeq, na.rm=T)
  
  
  
}

efdrMat <- data.frame(efdrMat)
colnames(efdrMat) <- 
  c("Oracle", "Total Count", "DESeq", "edgeR", "TCC", "PoissonSeq")

write.table(efdrMat, "efdrMatSymmDiffMRNA_prm2.txt")


efdrSDMat <- data.frame(efdrSDMat)
colnames(efdrSDMat) <- 
  c("Oracle", "Total Count", "DESeq", "edgeR", "TCC", "PoissonSeq")

write.table(efdrSDMat, "efdrSDMatSymmDiffMRNA_prm2.txt")




powerMat <- data.frame(powerMat)
colnames(powerMat) <-
  c("Oracle", "Total Count", "DESeq", "edgeR", "TCC", "PoissonSeq")

write.table(powerMat, "powerMatSymmDiffMRNA_prm2.txt")


powerSDMat <- data.frame(powerSDMat)
colnames(powerSDMat) <-
  c("Oracle", "Total Count", "DESeq", "edgeR", "TCC", "PoissonSeq")

write.table(powerSDMat, "powerSDMatSymmDiffMRNA_prm2.txt")


MSEMat <- data.frame(MSEMat)
colnames(MSEMat) <- 
  c("Oracle", "Total Count", "DESeq", "edgeR", "TCC", "PoissonSeq")

write.table(MSEMat, "MSEMatSymmDiffMRNA_prm2.txt")