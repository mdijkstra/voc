#
## Initialize
#
rm(list=ls(all=T))
library(stringr); library(pheatmap)
setwd('/Users/mdijkstra/Dropbox/Documents/collaboration/lude/voc/r/analysis')
source('functions.R') # source functions

#
## Load data
#
source('load_data.R') # p, voc

#
## Preprocessing
#

# Remove columns with >20% 0's.
index.remove = apply(voc, 2, function(vec) .2 < length(which(0 == vec))/length(vec))
voc = voc[, !index.remove]
voc[which(0 == voc)] = NA

# Constants
n.individuals	= nrow(p)
n.phenotypes	= ncol(p)
n.vocs			= ncol(voc)

#
## Calculate correlations
#

# Calculate cor/p-values for _real_ data
real = correlation(p, voc)

# Randomize pvoc matrix
p.randomized	= p[sample(1:n.individuals), ]

# Calculate cor/p-values for _random_ data
random	= correlation(p.randomized, voc)

# Extract only relevant p-values (i.e. Pheno x VOC)
p.values.real		= real$cor.matrix.p
p.values.random		= random$cor.matrix.p


##################################
#
## Which p-values are significant after Bonferroni correction?
## Which have FDR < 0.05?
#
fdr = fdr.quick(p.values.real, p.values.random, n.elements = 2e3)

if (F) {
	pdf('fdr.pdf')
		n.elements.zoom = 1e2
		plot(log(fdr$p.unique[1:n.elements.zoom]), fdr$fdr[1:n.elements.zoom], t='b', lwd=2, axes=F, xlab="log( p-value )", ylab="FDR")
		axis(1); axis(2, las=2)
		alpha.bonf = 0.05 / n.phenotypes / n.vocs
		abline(v = log(alpha.bonf), lty=2, lwd = 2) # Bonferroni
		abline(h = 0.05, lty=2, lwd = 2) # FDR
	dev.off()
}

# Select p-values for which FDR < .05
p.values.significant = fdr$p.unique[which(fdr$fdr < .05)]

# Find indices of these p-values
index = NULL
for (p.value in p.values.significant) index = rbind(index, which(p == real$cor.matrix.p, arr.in = T))










#