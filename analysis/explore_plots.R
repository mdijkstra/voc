#
## Log-trans may enable t-test which is more powerful than non-parametric spearman
#

#png('voc845_distribution.png')
#	hist(voc[,19],br=100, main = "Histogram of Compound 845")
#dev.off()

#png('voc845_log-distribution.png')
#	hist(log(voc[,19]), main = "Histogram of Compound 845, log scale")
#dev.off()

#png('67vocs_distribution.png')
#	hist(voc,br=1000,xlim=c(0,100))
#dev.off()

#png('67vocs_log-distribution.png')
#	hist(log(voc))
#dev.off()

n.zeros = apply(voc.all, 1, function(voc.i) length(which(0 == voc.i)) / length(voc.i))
hist(n.zeros,br=100)
hist(log(n.zeros))

# But if we use all vocs, we don't see a nice 'normal' distribution
png('allVocs_log-distribution.png')
	voc.all. = as.numeric(as.matrix(voc.all))
	voc.all.[0 == voc.all.] = NA
	hist(log(voc.all.))
dev.off()

shapiro.test(sample(log(voc.all.), 5e3))
shapiro.test(sample(log(voc), 5e3))