# bedtools genomecov function histogram output produces a file with 5 columns
# see http://www.quinlanlab.org/pdf/bedtools.protocols.pdf
# pages 11.12.9 - 11.12.10
# 1) chromosome
# 2) depth
# 3) number of base pairs with depth = column2
# 4) size of the chromosome
# 5) fraction of base pairs with depth = column2

#   To take command line arguments
args <- commandArgs(TRUE)
#   This creates a vector of character strings for arguments
#   we will just take two arguments here, the stats directory
#   and the sample name
statsdir <- args[1]
samplename <- args[2]

#    Plot type is histogram
type <- 'h'

#    Line width
lwd <- 5

#    Color of plot
col <- 'blue'

# Labels for x and y axis
xlab <- 'Depth'
ylab <- 'Fraction at genome at depth'

inputfile <- paste(statsdir, samplename, 'coverage.hist.txt')

cov <- read.table(inputfile)

plot(dat.frame,
	type=type,
	lwd=lwd,
	col=col,
    xlab=xlab,
    ylab=ylab)



outputfile <- paste(statsdir, "/plots/", samplename, "_Forward_SeqqsPlots.pdf", sep="")
pdf(file=outputfile, width=8.8, 11)


pdf(file =outputfile, width=6, height=6)
cov <- read.table('Sample_T7013_T5009_SCN_2015-03-12_exon.txt', header=F)
#gcov <- cov[cov[,1] == 'exon',]
plot(cov[1:51,2], cov[1:51,5], type='h', col='blue', lwd=5, xlab = 'Depth', ylab = 'Fraction of genome at depth')
dev.off()