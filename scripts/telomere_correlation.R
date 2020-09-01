#------------------------------------------------------------------------------------------------------
# This script is based on the following walkthrough:
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
#
# For testing replace my_data below with:
# my_data <- read.csv("https://wiki.q-researchsoftware.com/images/b/b9/Ownership.csv",
#                     header = TRUE, fileEncoding = "latin1")
#------------------------------------------------------------------------------------------------------

library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(psych)

#------------------------------------------------------------------------
# flattenCorrMatrix
#   cormat : matrix of the correlation coefficients
#   pmat : matrix of the correlation p-values
#------------------------------------------------------------------------

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

#------------------------------------------------------------------------
# panel.cor
#------------------------------------------------------------------------

panel.cor <- function(x, y) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits = 2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

args = commandArgs(trailingOnly = TRUE)

input_fname <- toString(args[1])
output_prefix <- toString(args[2])

pearson_pdf = paste0(output_prefix, "pearson.pdf")
spearman_pdf = paste0(output_prefix, "spearman.pdf")

pearson_txt = paste0(output_prefix, "pearson.txt")
spearman_txt = paste0(output_prefix, "spearman.txt")
telomere_txt = paste0(output_prefix, "telomere.txt")

heatmap_margins = c(10, 10)

my_data <- read.table(file = input_fname, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)

res <- cor(my_data, use = "complete.obs")
write.table(round(res, 2), file = telomere_txt)

res2.pearson <- rcorr(as.matrix(my_data), type = c("pearson"))
write.table(flattenCorrMatrix(res2.pearson$r, res2.pearson$P), file = pearson_txt)

res2.spearman <- rcorr(as.matrix(my_data), type = c("spearman"))
write.table(flattenCorrMatrix(res2.spearman$r, res2.spearman$P), file = spearman_txt)

#------------------------------------------------------------------------
# Plot pearson
#------------------------------------------------------------------------
pdf(pearson_pdf)

corrplot(res2.pearson$r, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
chart.Correlation(my_data, histogram = TRUE, method = c("pearson"), pch = 19)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res2.pearson$r, col = col, symm = TRUE, margins = heatmap_margins)

upper.panel<-function(x, y) {
  points(x,y, pch = 19)
}

pairs(my_data, lower.panel = panel.cor, upper.panel = upper.panel)

pairs.panels(my_data, method = "pearson", hist.col = "#00AFBB", density = TRUE, ellipses = TRUE)

dev.off()

#------------------------------------------------------------------------
# Plot pearson
#------------------------------------------------------------------------
pdf(spearman_pdf)

corrplot(res2.spearman$r, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
chart.Correlation(my_data, histogram = TRUE, method = c("spearman"), pch = 19)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res2.spearman$r, col = col, symm = TRUE, margins = heatmap_margins)

pairs.panels(my_data, method = "spearman", hist.col = "#00AFBB", density = TRUE, ellipses = TRUE)

dev.off()
