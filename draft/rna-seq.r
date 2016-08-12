library(ggplot2)
library(plyr)
library(baySeq)
library(edgeR)
library(pipelineRNASeq)
library(MASS)

## The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (11:11)
## functions to store

## getMDSPlot <- function(object, method = "bcv", groups, imgFile = NULL){
##     mdsStats <- plotMDS(object, method = method)
##     plotDf <- data.frame(x = mdsStats$x,
##                          y = mdsStats$y,
##                          group = groups)
##     ##
##     limit <- max(abs(plotDf$x), abs(plotDf$y))
##     labels <- row.names(plotDf)
##     p <- ggplot() +
##              geom_text(data = plotDf, aes(x, y, label = labels, color = group),
##                        size = 3) +
##              scale_color_manual(name = "Group", values = cbbPalette[1:3]) +   
##              xlab(paste(toupper(method), "distance 1")) + xlim(-limit, limit) + ylim(-limit, limit) +
##              ylab(paste(toupper(method), "distance 2")) + theme_bw() +
##              ggtitle("Multidimensional scaling plot")
##     if(is.null(imgFile)){
##         return(mdsStats)        
##     } else {
##         pdf(imgFile, width = 6, height = 5)
##         print(p)
##         dev.off()
##     }
## }

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

glmTest <- function(fitList, n){
    glmTopList <- llply(fitList, function(x) {
              data.frame(topTags(x, n))
          })
    names(glmTopList) <- laply(fitList, function(x) x$comparison)
    return(glmTopList)
}


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

smear <- function(df, fc.thresh = 0.5, alpha = 0.05){
    p <- ggplot() + 
         geom_point(data = subset(df, FDR > alpha),
                   aes(logCPM, logFC)) + 
         geom_point(data = subset(df, FDR < alpha),
                   aes(logCPM, logFC), color = cbbPalette[7]) +
         theme_bw() + ylim(-max(abs(df$logFC)), max(abs(df$logFC))) +
         geom_hline(yintercept = c(-fc.thresh, fc.thresh), color = cbbPalette[6]) +
         xlab("Average log CPM") +
         ylab(expression(paste("Foldchange [", log[2] ,"]")))    
    print(p)
}

## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (11:13)
## START

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

getMDSPlot(d, groups = d$samples$group, imgFile = "mds_1.pdf")



## ------------------------------------------------------------
## ------------------------------------------------------------
## Differenzial ANALYSIS
## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (10:43)
## dispersion

designMat <- model.matrix(~ 0 + d$samples$group)
colnames(designMat) <- levels(d$samples$group)
## estimate the dispersion by GLM
d <- estimateDispersion(d, designMat, method = "power")

## check plot if everything is okay
pdf("dispersionGLM.pdf", width = 6, height = 5)
estimateDispersion(d, designMat, method="power", figure_flag = TRUE)
dev.off()

## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (10:44)
## Diff Analysis

## fit the glm model for each gene i.e. tag
glmFitList <- glmFit(d, design = designMat, contrast = "Tukey") 

## get the top N test results
glmTestList <- glmTest(glmFitList, n = 500)

## plot all in a vulcano plot
pdf("vulcano.pdf", width = 5, height = 5)
l_ply(glmTestList, function(x) vulcano(x))
dev.off()

## plot all in a smear plot
pdf("smear.pdf", width = 5, height = 5)
l_ply(glmTestList, function(x) smear(x))
dev.off()

## ------------------------------------------------------------
## by J.Kruppa on Friday, June 26, 2015 (11:12)
## END





## ------------------------------------------------------------
## by J.Kruppa on Monday, June 29, 2015 (09:17)
## New analysis with the goseq package




library(pipelineRNASeq)

library(org.Hs.eg.db)
library(ggdendro)
library(goseq)


Li2008data


                                        # Get the mapping from ENSEMBL 2 Entrez
en2eg <- as.list(org.Hs.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg <- as.list(org.Hs.egPATH)
# Define a function which gets all unique KEGG IDs
# associated with a set of Entrez IDs
grepKEGG <- function(id, mapkeys){unique(unlist(mapkeys[id], use.names=FALSE))}
# Apply this function to every entry in the mapping from
# ENSEMBL 2 Entrez to combine the two maps
kegg <- lapply(en2eg, grepKEGG, eg2kegg)
head(kegg)










## ------------------------------------------------------------
## by J.Kruppa on Tuesday, December 16, 2014 (09:52)
## Get Annotations
## listDatasets(ensembl)[grep("mus", listDatasets(ensembl)$description),]

## load ensemble data
id_type <- "ensembl_gene_id"
mart <- useMart("ensembl")
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
attributes <- listAttributes(ensembl)

## get the annotation
annotation <- getBM(attributes = c(id_type, "description", "chromosome_name",
                        "go_id", "name_1006", "namespace_1003"),
                    filters = c(id_type), values = row.names(countDf), mart = ensembl)

## prepare annotation for writing
annotation$description <- sub(";", ",", annotation$description)
index <- match(rownames(countDf), annotation$ensembl_gene_id)
annoFinalOut <- annotation[index,]
write.csv(annoFinalOut, file.path(procDir, "GO2Genes_Annotations.csv"),
          quote=FALSE, row.names=FALSE)



## counts
d$counts

sigResList <- llply(resList, function(x) subset(x, padj < 1e-4)$geneID)

l_ply(seq_along(comparisonList), function(i) {
png(file.path(imgDir, paste0("heatmap-", names(comparisonList)[i],".png")),
    width = 15, height = 15, unit = "cm", res = 300)
    ## define colors
    coln <- 10
    red <- rgb(seq(255, 0, length.out=coln), rep(0, coln), rep(0, coln), maxColorValue=255)
    blue <- rgb(rep(0, coln), rep(0, coln), seq(0, 255, length.out=coln), maxColorValue=255)
    heatcol <- c(red[1:(coln-1)], blue)
    ## build matrix
    normCounts <- counts(cdsList[[i]], normalized=TRUE)
    mat <- as.matrix(normCounts[sigResList[[i]], ])
    colnames(mat) <- comparisonList[[i]]
    ## plot heatmap
    heatmap.2(mat, trace="none", col=heatcol, main = names(cdsList)[[i]], margins=c(3,11), scale = "row")    
dev.off()
}, .progress = "text")




goResList <- llply(resList, function(x){ 
    genes <- rep(0, dim(counts(cdsList[[1]], normalized=TRUE))[1])
    genes[x$padj < 0.05] = 1
    names(genes) = rownames(countDf)    
    pwf <- nullp(genes, genome = "mm10", id = "ensGene")    
    pvals <- goseq(pwf, genome = "mm10", id = "ensGene",
                   gene2cat = data.frame(annoFinalOut$ensembl_gene_id, annoFinalOut$go_id))
    padjusted <- p.adjust(pvals$over_represented_pvalue, method="BH")
    goResDf <- data.frame(GOcategory = pvals$category,
                          p.raw = pvals$over_represented_pvalue,
                          p.adj= padjusted)
    ##goResDf <- goResDf[-1,]
    def <- match(goResDf$GOcategory, annoFinalOut$go_id)
    finalGoResDf <- data.frame(goResDf, GOdescription = annoFinalOut$name_1006[def],
                               GOdomain = annoFinalOut$namespace_1003[def])
    finalGoResDf[order(finalGoResDf$p.raw),]
    return(finalGoResDf)
}, .progress = "text")

writeRDS(goResList, procDir)



pdf("pca.pdf", width = 7, height = 6)
PC <- prcomp(t(d$counts), center=FALSE, scale=FALSE)
PC2 <- predict(PC)
plotDf <- data.frame(pc1 = PC2[,1], pc2 = PC2[,2],
                     group = d$samples$group)
ggplot(plotDf, aes(pc1, pc2, label = row.names(plotDf))) + geom_text(size=3) +
    xlab("Principal Component 1") + ylab("Principal Component 2") +
    theme_bw()
ggplot(plotDf, aes(pc1, pc2)) + geom_point(aes(color = group)) +
    xlab("Principal Component 1") + ylab("Principal Component 2") +
    theme_bw()
## ggplot(plotDf, aes(pc1, pc2)) + geom_point(aes(color = timepoints)) +
##     xlab("Principal Component 1") + ylab("Principal Component 2") +
##     theme_bw()
## ggplot(plotDf, aes(pc1, pc2)) + geom_point(aes(color = adjuvant)) +
##     xlab("Principal Component 1") + ylab("Principal Component 2") +
##     theme_bw()
dev.off()


## ------------------------------------------------------------
## by J.Kruppa on Tuesday, December 16, 2014 (12:51)
## Hclust
distMat <- dist(t(d$counts))
clustSingleMat <- hclust(distMat, method = "single")
clustCompleteMat <- hclust(distMat, method = "complete")

imgDir <- "."

png(file.path(imgDir, "dendrogram-single.png"),
        width = 15, height = 15, unit = "cm", res = 300)
ggdendrogram(clustSingleMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'single linkage' - Single samples")
dev.off()

png(file.path(imgDir, "dendrogram-complete.png"),
        width = 15, height = 15, unit = "cm", res = 300)
ggdendrogram(clustCompleteMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'complete linkage' - Single samples")
dev.off()

pdf(file.path(imgDir, "Dendrogramms.pdf"))
ggdendrogram(clustSingleMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'single' - Single samples")
clustSingleMat$labels <- finalSampleInfoDf$group
ggdendrogram(clustSingleMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'single' - Experimental group")
clustSingleMat$labels <- finalSampleInfoDf$timepoints
ggdendrogram(clustSingleMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'single' - Timepoints")
ggdendrogram(clustCompleteMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'complete' - Single samples")
clustCompleteMat$labels <- finalSampleInfoDf$group
ggdendrogram(clustCompleteMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'complete' - Experimental group")
clustCompleteMat$labels <- finalSampleInfoDf$timepoints
ggdendrogram(clustCompleteMat, rotate=FALSE) +
    ggtitle("Dendrogramm with method 'complete' - Timepoints")
for(i in seq_along(comparisonList)){
    distMat <- dist(t(countDf[, as.numeric(names(comparisonList[[i]]))]))
    clustSingleMat <- hclust(distMat, method = "complete")
    clustSingleMat$labels <- comparisonList[[i]]
    p <- ggdendrogram(clustSingleMat, rotate=FALSE) +
        ggtitle(paste("Dendrogramm with method 'complete' - ", names(comparisonList)[[i]]))
    print(p)
}
dev.off()


