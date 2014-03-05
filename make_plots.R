#!/usr/bin/Rscript

for (suffix in c(".percentages", ".log2")) {
  titles <- read.table(paste("plot_titles", suffix, ".txt", sep=""), sep="\t")$V1

  png(filename=paste("plots", suffix, ".png", sep=""), width=500*3, height=500*ceiling(length(titles)/3), res=200)

  par(mfrow=c(ceiling(length(titles)/3),3))

  for (i in seq_along(titles))
  {
    data <- read.table(paste("plot_data", suffix, ".txt", sep=""), header=F, skip=(i-1)*2, nrow=2)
    fig <- plot(as.numeric(data[2,]), as.numeric(data[1,]), pch=20, col=rgb(0,0,0,0.25), main=titles[i], ylab="Expected", xlab="Actual")
    reg <- lm(as.numeric(data[1,]) ~ as.numeric(data[2,]))
    abline(reg)
    mtext(bquote(y == .(coef(reg)[2])*x + .(coef(reg)[1])), adj=1, padj=0, cex=0.8)
  }

  dev.off()
}

# rarefaction
#plot(log10(a$V1), a$V2-1.493295, pch=20, col=rgb(0,0,0,0.2), xlab="log10(Reads Sampled)", ylab="RMS Error", main="DNAexp4")
#lines(lowess(a$V2-1.493295 ~ log10(a$V1)))
