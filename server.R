
library(shiny)
library(ggplot2)
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
##
resDf <- readRDS("exampleGeneSets.RDS")

vulcano.plot <- function(df, fc.thresh = 0.5, alpha = 0.05){
    p <- ggplot() +
        geom_point(data = subset(df, padj > alpha &
                                     abs(log2FoldChange) < fc.thresh),
                   aes(log2FoldChange, -log10(padj)),
                   color = cbbPalette[1]) +        
        geom_point(data = subset(df, padj < alpha),
                   aes(log2FoldChange, -log10(padj)),
                   color = cbbPalette[7]) +
        geom_point(data = subset(df, abs(log2FoldChange) > fc.thresh),
                   aes(log2FoldChange, -log10(padj)),
                           color = cbbPalette[2]) +
        geom_point(data = subset(df, abs(log2FoldChange) > fc.thresh &
                                     padj < alpha),
                   aes(log2FoldChange, -log10(padj)),
                   color = cbbPalette[4]) +
        geom_text(data = subset(df, abs(log2FoldChange) > fc.thresh &
                                     padj < alpha),
                  aes(log2FoldChange, -log10(padj),
                      label = subset(df, abs(log2FoldChange) > fc.thresh &
                                     padj < alpha)$Gene),
                  check_overlap = TRUE, hjust = 0) +
        geom_vline(xintercept = c(-fc.thresh, fc.thresh),
                   linetype = "dotted") +
        geom_hline(yintercept = c(-log10(alpha)), linetype = "dotted") +
        xlab(expression(paste("Foldchange [", log[2] ,"]"))) +
        ylab("FDR") +
        theme_bw() 
    return(p)
}

shinyServer(function(input, output) {
    output$plots <- renderPlot(
        vulcano.plot(df = resDf,
                     fc.thresh = input$fc.thresh,
                     alpha = input$alpha.error)
    )
  output$view <- renderTable({
      hitDf <- subset(resDf, abs(log2FoldChange) > input$fc.thresh &
                             padj < input$alpha.error)
      orderedDf <- hitDf[order(abs(hitDf$log2FoldChange), decreasing = TRUE),]
      row.names(orderedDf) <- NULL
      names(orderedDf) <- c("Gene symbol", "log foldchange", "p-value", "p-value [adjusted]")
      orderedDf
  })
})
