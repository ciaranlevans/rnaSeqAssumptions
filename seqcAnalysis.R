require(seqc)

taqdat <- taqman

taqdat <- taqdat[-which(is.na(taqdat[,1])),]

taqdat <- taqdat[,1:18]

numPresent <- 
  (taqdat$A1_detection=="P") + (taqdat$A2_detection=="P") + (taqdat$A3_detection=="P") +
  (taqdat$A4_detection=="P") + (taqdat$B1_detection=="P") + (taqdat$B2_detection=="P") +
  (taqdat$B3_detection=="P") + (taqdat$B4_detection=="P") 

taqdat <- taqdat[which(numPresent>=6),]

taqdat <- taqdat[,-c(4,6,8,10,12,14,16,18)]

taqdat <- taqdat[-which(duplicated(taqdat$EntrezID)==T),]


readmat <- ILM_aceview_gene_AGR[,1:132]

commonGenes <- which(readmat$EntrezID %in% taqdat$EntrezID)

commonGenes <- readmat$EntrezID[commonGenes]

taqdat <- taqdat[which(taqdat$EntrezID %in% commonGenes),]
readmat <- readmat[which(readmat$EntrezID %in% commonGenes),]

taqdat <- taqdat[order(taqdat$EntrezID),]
readmat <- readmat[order(readmat$EntrezID),]

counts <- readmat[,5:132]
counts <- as.matrix(counts)
cond <- factor(rep(1:2, each=64))

counts[which(is.na(counts))] = 0

presentCounts <- which(rowSums(counts) >= 5)

counts <- counts[presentCounts,]
taqdat <- taqdat[presentCounts,]



logfc = abs(log(rowMeans(taqdat[,3:6]), base=2) - log(rowMeans(taqdat[,7:10]), base=2))
nonDE <- logfc < 0.2

logfc1 = log(rowMeans(taqdat[,3:6]), base=2) - log(rowMeans(taqdat[,7:10]), base=2)
logfc2 = log(rowMeans(taqdat[,7:10]), base=2) - log(rowMeans(taqdat[,3:6]), base=2)

#symmDiffMRNA <- union(which(logfc1 > -3.5 & logfc1 < -1), which(logfc1 > 3.5))
#symmDiffMRNA <- union(symmDiffMRNA, which(logfc1 > -1 & logfc1 < 1))

#countsSD <- counts[symmDiffMRNA,]
#taqdatSD <- taqdat[symmDiffMRNA,]

#logfcSD <- logfc[symmDiffMRNA]
#nonDESD <- logfcSD < 0.2


require(DESeq2)

# DESeq2
dds <- DESeqDataSetFromMatrix(counts, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
#sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

sum(res$padj[nonDE == T] < 0.05)/sum(res$padj < 0.05)

restrict <- which(logfc < 0.2 | logfc > 1)
logfcRestrict <- logfc[restrict]
padjRestrict <- res$padj[restrict]

nonDERestrict <- 0 + (logfcRestrict < 0.2)

rocDESeq2 <- roc((0 + nonDERestrict), padjRestrict)


# Total Count
require(psych)

totCounts <- colSums(counts)
sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")
sum(res$padj[nonDE == T] < 0.05)/sum(res$padj < 0.05)

padjRestrict <- res$padj[restrict]
rocTot <- roc((0 + nonDERestrict), padjRestrict)


# edgeR normalization
normFactors <- edgeR::calcNormFactors(counts)
normFactors <- normFactors*totCounts
normFactors <- normFactors/geometric.mean(normFactors)
sizeFactors(dds) <- normFactors
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

padjRestrict <- res$padj[restrict]
rocEdgeR <- roc((0 + nonDERestrict), padjRestrict)


## TCC estimate
tcc <- new("TCC", counts, cond)
tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                            iteration = 3, FDR = 0.1, floorPDEG = 0.05)
tccEsts <- tcc$norm.factors
tccEsts <- tccEsts*colSums(counts)/geometric.mean(colSums(counts))

sizeFactors(dds) <- tccEsts/geometric.mean(tccEsts)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

padjRestrict <- res$padj[restrict]
rocTCC <- roc((0 + nonDERestrict), padjRestrict)


## PoissonSeq estimate
psests <- PS.Est.Depth(counts)

sizeFactors(dds) <- psests/geometric.mean(psests)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

padjRestrict <- res$padj[restrict]
rocPS <- roc((0 + nonDERestrict), padjRestrict)


plot(rocDESeq2, main="ROC Curves for DE Testing")
lines(rocTot, col="blue", lty=2)
lines(rocEdgeR, col="green", lty=3)
lines(rocTCC, col="orange", lty=4)
lines(rocPS, col="pink", lty=5)
text(0.6, 0.55, "PoissonSeq", cex=0.7)
#text(0.65, 0.75, "DESeq2", cex=0.7)
legend("bottomright", c("DEGES", "DESeq", "PoissonSeq", "TMM", "Total Count"),
       lty=c(4,1,5,3,2), col=c("orange", "black", "pink", "green", "blue"),
       cex=c(0.7, 0.7, 0.7, 0.7, 0.7),
       lwd=2)

# sampling the columns
efdrDESeq2_s2 <- c()
efdrTot_s2 <- c()
efdrEdgeR_s2 <- c()
efdrPS_s2 <- c()
efdrTCC_s2 <- c()
set.seed(3)

for(i in 1:100){
  classA <- sample(1:64, 2, replace=F)
  classB <- sample(65:128, 2, replace=F)
  countsSub <- counts[,c(classA, classB)]
  condSub <- factor(rep(1:2, each=2))
  
  presentCountsSub <- which(rowSums(countsSub) >= 5)
  
  countsSub <- countsSub[presentCountsSub,]
  taqdatSub <- taqdat[presentCountsSub,]
  
  
  
  logfcSub = 
    abs(log(rowMeans(taqdatSub[,3:6]), base=2) - log(rowMeans(taqdatSub[,7:10]), base=2))
  nonDESub <- logfcSub < 0.2
  
  
  
  dds <- DESeqDataSetFromMatrix(countsSub, DataFrame(condSub), ~ condSub)
  dds <- estimateSizeFactors(dds)
  #sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  
  efdrDESeq2_s2[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  # Total Count
  totCounts <- colSums(countsSub)
  sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrTot_s2[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  # edgeR normalization
  normFactors <- edgeR::calcNormFactors(countsSub)
  normFactors <- normFactors*totCounts
  normFactors <- normFactors/geometric.mean(normFactors)
  sizeFactors(dds) <- normFactors
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrEdgeR_s2[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  ## TCC estimate
  tcc <- new("TCC", countsSub, condSub)
  tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                              iteration = 3, FDR = 0.1, floorPDEG = 0.05)
  tccEsts <- tcc$norm.factors
  tccEsts <- tccEsts*colSums(countsSub)/geometric.mean(colSums(countsSub))
  
  sizeFactors(dds) <- tccEsts/geometric.mean(tccEsts)
  
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrTCC_s2[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  ## PoissonSeq estimate
  psests <- PS.Est.Depth(countsSub)
  
  sizeFactors(dds) <- psests/geometric.mean(psests)
  
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrPS_s2[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  print(i)
}



# sampling the columns
efdrDESeq2_s5 <- c()
efdrTot_s5 <- c()
efdrEdgeR_s5 <- c()
efdrPS_s5 <- c()
efdrTCC_s5 <- c()
set.seed(3)

for(i in 1:100){
  classA <- sample(1:64, 5, replace=F)
  classB <- sample(65:128, 5, replace=F)
  countsSub <- counts[,c(classA, classB)]
  condSub <- factor(rep(1:2, each=5))
  
  presentCountsSub <- which(rowSums(countsSub) >= 5)
  
  countsSub <- countsSub[presentCountsSub,]
  taqdatSub <- taqdat[presentCountsSub,]
  
  
  
  logfcSub = 
    abs(log(rowMeans(taqdatSub[,3:6]), base=2) - log(rowMeans(taqdatSub[,7:10]), base=2))
  nonDESub <- logfcSub < 0.2
  
  
  
  dds <- DESeqDataSetFromMatrix(countsSub, DataFrame(condSub), ~ condSub)
  dds <- estimateSizeFactors(dds)
  #sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  
  efdrDESeq2_s5[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  # Total Count
  totCounts <- colSums(countsSub)
  sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrTot_s5[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  # edgeR normalization
  normFactors <- edgeR::calcNormFactors(countsSub)
  normFactors <- normFactors*totCounts
  normFactors <- normFactors/geometric.mean(normFactors)
  sizeFactors(dds) <- normFactors
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrEdgeR_s5[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  ## TCC estimate
  tcc <- new("TCC", countsSub, condSub)
  tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                              iteration = 3, FDR = 0.1, floorPDEG = 0.05)
  tccEsts <- tcc$norm.factors
  tccEsts <- tccEsts*colSums(countsSub)/geometric.mean(colSums(countsSub))
  
  sizeFactors(dds) <- tccEsts/geometric.mean(tccEsts)
  
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrTCC_s5[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  
  ## PoissonSeq estimate
  psests <- PS.Est.Depth(countsSub)
  
  sizeFactors(dds) <- psests/geometric.mean(psests)
  
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  res$padj <- p.adjust(res$pvalue, method="BH")
  efdrPS_s5[i] <- sum(res$padj[nonDESub == T] < 0.05)/sum(res$padj < 0.05)
  
  print(i)
}


