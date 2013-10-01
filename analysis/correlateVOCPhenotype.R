#
## Initialize
#
rm(list=ls(all=F))
library(stringr); library(pheatmap)
setwd('/Users/mdijkstra/Dropbox/Documents/collaboration/lude/voc/analysis')
source('functions.R') # source functions

# Define unique stamp which we put in file name
non.zero.fraction = .5 # at least
select.gender = "Male" #Male, Female, All
ignore.zeros = T
unique.stamp = str_c(select.gender, "_", "Nonzerofraction", non.zero.fraction, if (ignore.zeros) "_IgnoreZeros" else "", "_", str_replace_all(format(Sys.time(), "%X"), ":", ""))

#
## Load data
#
source('load_data.R') # p, voc

#
## Preprocessing
#

# Gender selection:
index.male		= which(1 == p[, "GESLACHT"])
index.female	= which(2 == p[, "GESLACHT"])

if ("Male" == select.gender)	individuals.keep = index.male
if ("Female" == select.gender)	individuals.keep = index.female
if ("All" == select.gender)		individuals.keep = 1:nrow(p)

p	= p[individuals.keep, ]
voc	= voc[individuals.keep, ]

# Remove columns with >20% 0's. Also try >50%
index.keep = apply(voc, 2, function(vec) non.zero.fraction <= (1 - length(which(0 == vec))/length(vec)))
voc = voc[, index.keep]
if (ignore.zeros) voc[which(0 == voc)] = NA

# Constants
n.individuals	= nrow(p)
n.phenotypes	= ncol(p)
n.vocs			= ncol(voc)

#
## Calculate correlations
#

# Calculate cor/p-values for _real_ data
real = correlation.two.part(p, voc)

# Randomize pvoc matrix
p.randomized	= p[sample(1:n.individuals), ]

# Calculate cor/p-values for _random_ data
random	= correlation.two.part(p.randomized, voc)

# Extract only relevant p-values (i.e. Pheno x VOC)
p.values.real		= real$cor.matrix.p
p.values.random		= random$cor.matrix.p


##################################
#
## Which p-values are significant after Bonferroni correction?
## Which have FDR < 0.05?
#
fdr = fdr.quick(p.values.real, p.values.random, n.elements = 2e3)

# Select p-values for which FDR < .05
index.fdr = max(which(fdr$fdr <= .05))
hits = !is.infinite(index.fdr)
p.values.significant = if (!hits) NULL else fdr$p.unique[1:index.fdr] #which(fdr$fdr <= .05)]

if (T) {
	pdf(str_c('plots/', 'fdr_', unique.stamp, '.pdf'))
		n.elements.zoom = 1e2
		plot(log(fdr$p.unique[1:n.elements.zoom], 10), fdr$fdr[1:n.elements.zoom], t='b', pch=19, lwd=2, axes=F, xlab="log( p-value )", ylab="FDR", col="gray", main = str_c("#signif = ", if (hits) index.fdr else 0))
		axis(1); axis(2, las=2)
		alpha.bonf = 0.05 / n.phenotypes / n.vocs
		abline(v = log(alpha.bonf, 10), lty=2, lwd = 2) # Bonferroni
		abline(h = 0.05, lty=2, lwd = 2) # FDR
		
		if (hits)
		{
			par(new = T)
			plot(log(fdr$p.unique[1:n.elements.zoom], 10), fdr$fdr[1:n.elements.zoom], t = 'p', pch=c(rep(1, index.fdr), rep(NA, n.elements.zoom - index.fdr)), col="red", axes=F, xlab="", ylab="", lwd = 2)	
		}
	dev.off()
}

# Find indices of these p-values
summary = NULL
for (p.value in p.values.significant) {
	index = which(p.value == real$cor.matrix.p, arr.in = T)
	if (1 == nrow(index)) {
		summary = rbind(summary, c(pheno = rownames(index), voc = colnames(voc)[index[,2]], p.value = p.value))
	}
}

summary.pheno	= unique(summary[, "pheno"])
summary.voc		= unique(summary[, "voc"])
summary.matrix	= matrix(NA, nrow = length(summary.pheno), ncol = length(summary.voc), dimnames = list(summary.pheno, summary.voc))
if (hits) for (i in 1:nrow(summary)) {
	summary.matrix[summary[i, "pheno"], summary[i, "voc"]] = summary[i, "p.value"]
}
summary.matrix = apply(summary.matrix, 1, as.numeric)
if (is.numeric(summary.matrix)) {
	summary.matrix = as.matrix(summary.matrix, nrow = 1, ncol = 1)
	colnames(summary.matrix) = summary.pheno
}
rownames(summary.matrix) = summary.voc


if (hits & 1 < prod(dim(summary.matrix))) pheatmap(log(t(summary.matrix), 10), cluster_row = F, cluster_col = F, filename = str_c("plots/", "phenoVsVocsPvalueOnLog10Scale_", unique.stamp, ".pdf"))






#