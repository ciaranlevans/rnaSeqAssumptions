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

set.seed(3)
# sample 89 so that 75% of DE genes have LFC > 1 and other 25% have LFC < -1
asymm <- union(which(logfc1 > -1), sample(which(logfc1 < -1), 89, replace=F))

countsAsymm <- counts[asymm,]
taqdatAsymm <- taqdat[asymm,]

logfcAsymm <- logfc[asymm]
nonDEAsymm <- logfcAsymm < 0.2


require(DESeq2)

# DESeq2
dds <- DESeqDataSetFromMatrix(countsAsymm, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
#sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

#sum(res$padj[nonDE == T] < 0.05)/sum(res$padj < 0.05)

restrict <- which(logfcAsymm < 0.2 | logfcAsymm > 1)
logfcRestrict <- logfcAsymm[restrict]
padjRestrict <- res$padj[restrict]

nonDERestrict <- 0 + (logfcRestrict < 0.2)

rocDESeq2 <- roc((0 + nonDERestrict), padjRestrict)


# Total Count
require(psych)

totCounts <- colSums(countsAsymm)
sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")
#sum(res$padj[nonDE == T] < 0.05)/sum(res$padj < 0.05)

padjRestrict <- res$padj[restrict]
rocTot <- roc((0 + nonDERestrict), padjRestrict)


# edgeR normalization
normFactors <- edgeR::calcNormFactors(countsAsymm)
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
tcc <- new("TCC", countsAsymm, cond)
tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                            iteration = 3, FDR = 0.1, floorPDEG = 0.05)
tccEsts <- tcc$norm.factors
tccEsts <- tccEsts*colSums(countsAsymm)/geometric.mean(colSums(countsAsymm))

sizeFactors(dds) <- tccEsts/geometric.mean(tccEsts)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

padjRestrict <- res$padj[restrict]
rocTCC <- roc((0 + nonDERestrict), padjRestrict)


## PoissonSeq estimate
psests <- PS.Est.Depth(countsAsymm)

sizeFactors(dds) <- psests/geometric.mean(psests)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")

padjRestrict <- res$padj[restrict]
rocPS <- roc((0 + nonDERestrict), padjRestrict)


plot(rocDESeq2, main="ROC Curves for DE Testing, Asymmetry")
lines(rocTot, col="blue", lty=2)
lines(rocEdgeR, col="green", lty=3)
lines(rocTCC, col="orange", lty=4)
lines(rocPS, col="pink", lty=5)
text(0.87, 0.65, "Total Count", cex=0.7)
#text(0.65, 0.75, "DESeq2", cex=0.7)
legend("bottomright", c("DEGES", "DESeq", "PoissonSeq", "TMM", "Total Count"),
       lty=c(4,1,5,3,2), col=c("orange", "black", "pink", "green", "blue"),
       cex=c(0.7, 0.7, 0.7, 0.7, 0.7),
       lwd=2)