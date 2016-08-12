library(ggplot2)
library(plyr)
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

res <- read.table("results.txt", header=TRUE)
head(res)

saveRDS(res, "exampleGeneSets.RDS")

vulcano.plot <- function(df, fc.thresh = 0.5, alpha = 0.05){
    p <- ggplot() +
         geom_point(data = subset(df, padj > alpha & abs(log2FoldChange) < fc.thresh),
                    aes(log2FoldChange, -log10(padj)),
                    color = cbbPalette[1]) +        
         geom_point(data = subset(df, padj < alpha),
                    aes(log2FoldChange, -log10(padj)),
                    color = cbbPalette[7]) +
         geom_point(data = subset(df, abs(log2FoldChange) > fc.thresh),
                    aes(log2FoldChange, -log10(padj)),
                    color = cbbPalette[2]) +
         geom_point(data = subset(df, abs(log2FoldChange) > fc.thresh & padj < alpha),
                    aes(log2FoldChange, -log10(padj)),
                    color = cbbPalette[4]) +
         geom_vline(xintercept = c(-fc.thresh, fc.thresh), linetype = "dotted") +
         geom_hline(yintercept = c(-log10(alpha)), linetype = "dotted") +
         xlab(expression(paste("Foldchange [", log[2] ,"]"))) +
         ylab("FDR") +
         theme_bw() 
    return(p)
}

png("vulcanoplot.png", width = 12, height = 12, res = 300, unit = "cm")
vulcano.plot(df = res)
dev.off()

png("foldchange1.png", width = 12, height = 12, res = 300, unit = "cm")
ggplot(res, aes(x=log2FoldChange)) +
    geom_histogram(aes(y = (..count..)/sum(..count..)), fill = cbbPalette[2]) +
    theme_bw() + ylab("Relative Häufigkeit") +
    xlab(expression(paste("Foldchange [", log[2] ,"]")))
dev.off()
#

png("foldchange2.png", width = 12, height = 12, res = 300, unit = "cm")
ggplot(res, aes(x=exp(log2FoldChange))) +
    geom_histogram(aes(y = (..count..)/sum(..count..)), fill = cbbPalette[2]) +
    theme_bw() + ylab("Relative Häufigkeit") +
    xlab("Foldchange")
dev.off()

exp(res$log2FoldChange)


# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
#
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#
# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))




## ------------------------------------------------------------
## by J.Kruppa on Thursday, June 25, 2015 (14:40)
## edgeR
## https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

library(baySeq)
library(edgeR)
data(mobData)
data(mobAnnotation)

treatment <- c("MM", "MM", "WM", "WM", "WW", "WW")

d <- DGEList(counts = mobData, group = factor(treatment))

## ------------------------------------------------------------
## by J.Kruppa on Thursday, June 25, 2015 (15:05)
## remove to low counts

d.full <- d # keep the old one in case we mess up
head(d$counts)

apply(d$counts, 2, sum) # total gene counts per sample

keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

## reset the lib size
d$samples$lib.size <- colSums(d$counts)
d$samples

## http://www.rna-seqblog.com/which-method-should-you-use-for-normalization-of-rna-seq-data/
d <- calcNormFactors(d)


## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (10:43)
## dispersion

designMat <- model.matrix(~ 0 + d$samples$group)
colnames(designMat) <- levels(d$samples$group)

pdf("dispersionGLM.pdf", width = 6, height = 5)
estimateDispersion(d, designMat, method="power", figure_flag = TRUE)
dev.off()


d <- estimateDispersion(d, designMat, method = "power")


fit <- glmFit(d, designMat)

## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (10:44)
## Diff Analysis


glmFit <- function(object, design, contrast = "Tukey") {
    require(multcomp)
    contrastMat <- contrMat(rep(1, ncol(design)), type = contrast)
    colnames(contrastMat) <- levels(groups)
    ## gloabl fit
    fit <- edgeR::glmFit(object, design)
    ## by contrast
    fitList <- llply(1:nrow(contrastMat), function(i) {
                         glmLRT(fit, contrast = contrastMat[i,])
                     })
    return(fitList)
}


glmFitList <- glmFit(d, design = designMat, contrast = "Tukey") 

glmTest <- function(fitList, n){
    glmTopList <- llply(fitList, function(x) {
              data.frame(topTags(x, n))
          })
    names(glmTopList) <- laply(fitList, function(x) x$comparison)
    return(glmTopList)
}

glmTestList <- glmTest(glmFitList, n = 500)



vulcano <- function(df, fc.thresh = 0.5, alpha = 0.05){
    p <- ggplot() +
         geom_point(data = subset(df, FDR > alpha & abs(logFC) < fc.thresh),
                    aes(logFC, -log10(FDR)),
                    color = cbbPalette[1]) +        
         geom_point(data = subset(df, FDR < alpha),
                    aes(logFC, -log10(FDR)),
                    color = cbbPalette[7]) +
         geom_point(data = subset(df, abs(logFC) > fc.thresh),
                    aes(logFC, -log10(FDR)),
                    color = cbbPalette[2]) +
         geom_point(data = subset(df, abs(logFC) > fc.thresh & FDR < alpha),
                    aes(logFC, -log10(FDR)),
                    color = cbbPalette[4]) +
         geom_vline(xintercept = c(-fc.thresh, fc.thresh), linetype = "dotted") +
         geom_hline(yintercept = c(-log10(alpha)), linetype = "dotted") +
         xlab(expression(paste("Foldchange [", log[2] ,"]"))) +
         ylab("FDR") +
         theme_bw() 
    print(p)
}

pdf("vulcano.pdf", width = 5, height = 5)
l_ply(glmTestList, function(x) vulcano(x))
dev.off()




glmTopList[[1]][1]

laply(fitList, function(x) x$comparison)

glmTest(glmFitList, n = 10)

contrMat(rep(1, length(levels(groups))), type = "Tukey")

fit <- glmFit(d, design.mat)
# compare (group 1 - group 2) to 0:
# this is equivalent to comparing group 1 to group 2
lrt12 <- glmLRT(fit, contrast=c(1,-1,0))
lrt13 <- glmLRT(fit, contrast=c(1,0,-1))
lrt23 <- glmLRT(fit, contrast=c(0,1,-1))
topTags(lrt12, n=10)











d <- DGEList(counts = mobData, group = factor(treatment))
d


design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)

d <- estimateGLMCommonDisp(d, design.mat)

fit <- glmFit(d, designMat)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
##
fit <- glmFit(d2, design.mat)



## ------------------------------------------------------------
## by J.Kruppa on Thursday, June 25, 2015 (15:05)
## remove to low counts

d.full <- d # keep the old one in case we mess up
head(d$counts)

apply(d$counts, 2, sum) # total gene counts per sample

keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

## reset the lib size
d$samples$lib.size <- colSums(d$counts)
d$samples

## http://www.rna-seqblog.com/which-method-should-you-use-for-normalization-of-rna-seq-data/
d <- calcNormFactors(d)
d

## MDS plot
getMDSPlot <- function(object, method = "bcv", groups, figure_flag = FALSE){
    mdsStats <- plotMDS(object, method = method)
    plotDf <- data.frame(x = mdsStats$x,
                         y = mdsStats$y,
                         group = groups)
    ##
    limit <- max(abs(plotDf$x), abs(plotDf$y))
    p <- ggplot() +
             geom_text(data = plotDf, aes(x, y, label = row.names(plotDf), color = group),
                       size = 3) +
             scale_color_manual(name = "Group", values = cbbPalette[1:3]) +   
             xlab(paste(toupper(method), "distance 1")) + xlim(-limit, limit) + ylim(-limit, limit) +
             ylab(paste(toupper(method), "distance 2")) + theme_bw() +
             ggtitle("Multidimensional scaling plot")
    if(figure_flag){
        return(p)
    } else {
        return(mdsStats)
    }
}

mdsPlot <- getMDSPlot(d, method = "bcv", groups = groups, figure_flag = TRUE)
##
pdf("mds.pdf", width = 8, height = 7)
mdsPlot
dev.off()


designMat <- model.matrix(~ 0 + treatment)
colnames(designMat) <- levels(treatment)


estimateDispersion <- function(object, design, method = "power", figure_flag = FALSE){
    object <- estimateGLMCommonDisp(object, design)
    object <- estimateGLMTrendedDisp(object, design, method = method)
    object <- estimateGLMTagwiseDisp(object, design)
    plotDf <- data.frame(tagwise = object$tagwise.dispersion,
                         trendwise = object$trended.dispersion,
                         aveLogCPM = object$AveLogCPM)
    ##
    p <- ggplot() + geom_point(data = plotDf, aes(aveLogCPM, tagwise, color = "Gene")) +
            geom_line(data = plotDf, aes(aveLogCPM, trendwise, color = "Trend")) +
            geom_hline(yintercept = object$common.dispersion, color = cbbPalette[7]) +
            theme_bw() +
            scale_color_manual(name = "Method",
                               values = c("Gene" = cbbPalette[1], "Trend" = cbbPalette[6])) +   
            xlab("Average log CPM") + ylab("Biological coeffcient of variation") +
            ggtitle("Estimatated Dispersion")
    if(figure_flag){
        return(p)
    } else {
        return(object)
    }
}

pdf("dispersionGLM.pdf", width = 6, height = 5)
estimateDispersion(d, designMat, method="power", figure_flag = TRUE)
dev.off()


d <- estimateDispersion(d, designMat, method = "power")

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)

fit <- glmFit(d, design.mat)
lrt12 <- glmLRT(fit, contrast = c(1, -1, 0))
lrt13 <- glmLRT(fit, contrast = c(1, 0, -1))
lrt23 <- glmLRT(fit, contrast = c(0, 1, -1))
topTags(lrt12, n=10)



et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
et13 <- exactTest(d1, pair=c(1,3)) # compare groups 1 and 3
et23 <- exactTest(d1, pair=c(2,3)) # compare groups 2 and 3
topTags(et12, n=10)

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)



require(multcomp)

glmFit <- function(object, design, groups, contrast = "Tukey") {
    contrastMat <- contrMat(rep(1, length(levels(groups))), type = contrast)
    colnames(contrastMat) <- levels(groups)
    ## gloabl fit
    fit <- glmFit(object, design)
    ## by contrast
    fitList <- llply(1:nrow(contrastMat), function(i) {
                         glmLRT(fit, contrast = contrastMat[i,])
                     })
    return(fitList)
}


glmFitList <- glmFit(d, design = designMat, groups, contrast = "Tukey") 

contrMat(rep(1, length(levels(groups))), type = "Tukey")

fit <- glmFit(d2, design.mat)
# compare (group 1 - group 2) to 0:
# this is equivalent to comparing group 1 to group 2
lrt12 <- glmLRT(fit, contrast=c(1,-1,0))
lrt13 <- glmLRT(fit, contrast=c(1,0,-1))
lrt23 <- glmLRT(fit, contrast=c(0,1,-1))
topTags(lrt12, n=10)




### Old stuff-------------










## ------------------------------------------------------------
## by J.Kruppa on Thursday, June 25, 2015 (15:10)
## dispersion

## over all
d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)

plotBCV(d1)

## the whole thing with a GLM model
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
