# Rscript CMplot.R data.txt
library(CMplot)

# 从命令行参数中获取文件名
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]

# 读取数据文件
data <- read.table(file_name, header = TRUE)

# 绘制密度图
# pdf("density_plot.pdf")
CMplot(data, plot.type="d", col=c("darkgreen","yellow","red"), file="pdf", dpi=300, file.output=TRUE, verbose=TRUE)
# dev.off()

# 绘制 Q-Q 图
# 每个性状一个图
# pdf("multracks_Q-Q_plot.pdf")
CMplot(data, plot.type="q", col=c("dodgerblue1","olivedrab3","darkgoldenrod1"), threshold=1e6, conf.int.col="grey", box=FALSE, multracks=TRUE, file="pdf", dpi=300, file.output=TRUE, verbose=TRUE)
# dev.off()
# 几个性状放在一个图里
# pdf("multraits_Q-Q_plot.pdf")
CMplot(data, plot.type="q", col=c("dodgerblue1","olivedrab3","darkgoldenrod1"), threshold=1e6, conf.int.col="grey", box=FALSE, multraits=TRUE, file="pdf", dpi=300, file.output=TRUE, verbose=TRUE)
# dev.off()

# 绘制曼哈顿图
# 每个性状一个图
# pdf("multracks_Manhtn_plot.pdf")
CMplot(data, plot.type="m", multracks=T, threshold=c(1e-6,1e-4), threshold.col=c("red","blue"), threshold.lty=c(1,2), amplify=T, signal.cex=c(1,1), signal.pch=c(20,20), signal.col=c("red","orange"), chr.den.col=c("darkgreen","yellow","red"), file="pdf", dpi=300, file.output=TRUE, verbose=TRUE)
# dev.off()
# 几个性状的曼哈顿图放在一个图里
# pdf("multraits_Manhtn_plot.pdf")
CMplot(data, plot.type="m", multraits=T, threshold=c(1e-6,1e-4), threshold.col=c("red","blue"), threshold.lty=c(1,2), amplify=T, signal.cex=c(1,1), signal.pch=c(20,20), signal.col=c("red","orange"), chr.den.col=c("darkgreen","yellow","red"), file="pdf", dpi=300, file.output=TRUE, verbose=TRUE)
# dev.off()

# 绘制多圈环状图
# pdf("circle_Manhtn_plot.pdf")
CMplot(data, plot.type="c", r=0.4, col=c("grey30","grey60"), threshold=c(1e-6,1e-4), threshold.col=c("red","blue"), amplify=TRUE, signal.col=c("red","orange"), chr.labels=paste("Chr",c(1:20),sep=""), chr.den.col=c("darkgreen","yellow","red"), cir.chr.h=1.5, outward=FALSE, file="pdf", dpi=300, file.output=TRUE, verbose=TRUE)
# dev.off()
