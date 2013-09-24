
par(ask=T)
for (i in 1:nrow(index))
{
	x			= p[, index[i, 1]]
	y			= voc[, index[i, 2]]
	xy.cor		= round(cor.test(x, y, use = "complete.obs", method = "spearman")$estimate, 2)
	xy.cor.p	= cor.test(x, y, use = "complete.obs", method = "spearman")$p.value
	plot(x, y, xlab = rownames(index)[i], ylab = colnames(voc)[index[i, 2]], main = str_c("Cor: ", xy.cor,", p-val: ", xy.cor.p), axes=F)
	axis(1)
	axis(2, las=2)
	abline(lsfit(x,y))
}
par(ask=F)

voc.subset	= voc[, unique(index[,2])]
voc.cor		= correlation(voc.subset, voc.subset)


# Show that 1 == male and 2 == female
male = 1 == p[,"GESLACHT"]
female = 2 == p[,"GESLACHT"]
pdf('male.pdf')
	boxplot(p[male,"LENGTE"], p[female,"LENGTE"], las = 2, axes = F)
	axis(2, las=2)
dev.off()

# Investigate most signif. result!
pdf('signResult.pdf')
	boxplot(voc[male, 19], voc[female, 19], ylim=c(0,60), axes=F, xlab="gender", ylab=colnames(voc)[19])
	axis(2, las=2)
dev.off()

# Idem for next sign. result
pdf('signResult2.pdf')
	i = 2
	x			= p[, index[i, 1]]
	y			= voc[, index[i, 2]]
	xy.cor		= round(cor.test(x, y, use = "complete.obs", method = "spearman")$estimate, 2)
	xy.cor.p	= cor.test(x, y, use = "complete.obs", method = "spearman")$p.value

	plot(x, y, xlab = rownames(index)[i], ylab = colnames(voc)[index[i, 2]], main = str_c("Cor: ", xy.cor,", p-val: ", xy.cor.p), axes=F)
	axis(1)
	axis(2, las=2)
	abline(lsfit(x,y))
dev.off()