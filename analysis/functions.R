correlation = function(p, voc){
	cor.matrix.cor	= NULL # [phenotype.i, voc.j] correlation
	cor.matrix.p	= NULL # [phenotype.i, voc.j] p-value

	cor.test.PiVOCj = function(VOCj, Pi) cor.test(VOCj, Pi, use = "complete.obs", method = "spearman")
	n.phenotypes = ncol(p)
	for (i in 1:n.phenotypes)
	{	
		this.cor.list	= apply(voc, 2, cor.test.PiVOCj, p[, i])
		cor.values		= sapply(this.cor.list, function(v) v$estimate)
		p.values		= sapply(this.cor.list, function(v) v$p.value)
		cor.matrix.cor	= rbind(cor.matrix.cor, cor.values)
		cor.matrix.p	= rbind(cor.matrix.p, p.values)

		print(str_c(round(i / n.phenotypes * 100), "% completed..."))
	}
	
	rownames(cor.matrix.cor) = rownames(cor.matrix.p) = colnames(p)
	colnames(cor.matrix.cor) = colnames(cor.matrix.p) = colnames(voc)
	
	list(cor.matrix.cor = cor.matrix.cor, cor.matrix.p = cor.matrix.p)
}

correlation.two.part = function(p, voc){
	n.phenotypes = ncol(p)
	cor.matrix.cor	= NULL # [phenotype.i, voc.j] correlation
	cor.matrix.p	= NULL # [phenotype.i, voc.j] p-value

	# Determine which vocs are zero
	voc.is.zero = voc
	voc.is.zero = apply(t(voc), 1, function(x) 0 == x | is.na(x))

	cor.test.PiVOCj = function(VOCj, Pi) cor.test(VOCj, Pi, use = "complete.obs", method = "spearman")
	wilcox.test.PiVOCj = function(zero, Pi) {
		if (0 == length(which(zero)) | 0 == length(which(!zero)) | all(is.na(Pi[zero])) | all(is.na(Pi[!zero]))) {
			return(1)
		} else {
			pw = wilcox.test(Pi[zero], Pi[!zero])$p.value
			if (is.finite(pw)) return(pw) else return(1)
		}
	}
	
	for (i in 1:n.phenotypes)
	{
		#print(str_c("Phenotype ", i))
		this.cor.list		= apply(voc, 2, cor.test.PiVOCj, p[, i])
		cor.values			= sapply(this.cor.list, function(v) v$estimate)
		spearman.p.values	= sapply(this.cor.list, function(v) v$p.value)
		
		wilcox.p.values	= apply(voc.is.zero, 2, wilcox.test.PiVOCj, p[, i])
		
		
		cor.matrix.cor	= rbind(cor.matrix.cor, cor.values)
		cor.matrix.p	= rbind(cor.matrix.p, spearman.p.values * wilcox.p.values)

		print(str_c(round(i / n.phenotypes * 100), "% completed..."))
	}
	
	rownames(cor.matrix.cor) = rownames(cor.matrix.p) = colnames(p)
	colnames(cor.matrix.cor) = colnames(cor.matrix.p) = colnames(voc)
	
	list(cor.matrix.cor = cor.matrix.cor, cor.matrix.p = cor.matrix.p)
}

fdr.quick = function(p.values.real, p.values.random, n.elements = 2e3)
{
	a	= sort(as.vector(p.values.real[which(!is.na(p.values.real))]))
	b	= sort(as.vector(p.values.random[which(!is.na(p.values.random))]))
	
	p.unique	= NULL
	fdr			= NULL
	i = j = 1
	while ((i < length(a) | j < length(b)) & length(p.unique) < n.elements)
	{
		if (a[i] < b[j]) {
			while (a[i] == a[i + 1] & i < length(a)) i = i + 1
			p.unique	= c(p.unique, a[i])
			fdr			= c(fdr, min(1, (j - 1) / i))
			i = i + 1
		} else if (b[j] < a[i]) {
			while (b[j] == b[j + 1] & j < length(b)) j = j + 1
			p.unique	= c(p.unique, b[j])
			fdr			= c(fdr, min(1, j / (i - 1)))
			j = j + 1
		} else if (a[i] == b[j])
		{
			while (a[i] == a[i + 1] & i < length(a)) i = i + 1
			while (b[j] == b[j + 1] & j < length(b)) j = j + 1
			p.unique	= c(p.unique, a[i])
			fdr			= c(fdr, min(1, j / i))
			
			i = i + 1
			j = j + 1
		}
	}
	
	list(p.unique = p.unique, fdr = fdr)
}

fdr.manyIterations = function(p.values.real, p.values.random, n.random.iterations, n.elements = 2e3)
{
	a	= sort(as.vector(p.values.real[which(!is.na(p.values.real))]))
	b	= sort(as.vector(p.values.random[which(!is.na(p.values.random))]))
	
	p.values	= NULL
	fdr			= NULL
	i.max = min(length(a), n.elements)
	for (i in 1:i.max)
	{
		fdr[i] = length(which(b <= a[i])) / n.random.iterations / i
	}

	list(p = a[1:i.max], fdr = fdr)
}
