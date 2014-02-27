png(filename="plots.png", width=2000, height=2000, res=200)

par(mfrow=c(4,3))

titles <- read.table("plot_titles.txt", sep="\t")$V1

for (i in seq_along(titles))
{
  data <- read.table("plot_data.txt", header=F, skip=(i-1)*2, nrow=2)
  fig <- plot(as.numeric(data[2,]), as.numeric(data[1,]), pch=20, col=rgb(0,0,0,0.25), main=titles[i], ylab="Expected", xlab="Actual")
  reg <- lm(as.numeric(data[1,]) ~ as.numeric(data[2,]))
  abline(reg)
  mtext(bquote(y == .(coef(reg)[2])*x + .(coef(reg)[1])), adj=1, padj=0, cex=0.8)
}

dev.off()
