#
## Load data
#
#p.all		= read.csv('/Volumes/lldeep/data/PhenotypeData2017-08-06BloodCellCountParameters.txt', sep='\t', row.names=1)
if (exists(".p.all"))		p.all = .p.all				else .p.all =		p.all		= dget('/Volumes/lldeep/data/phenotypeDataBloodCellCountParametersAll.RData')
if (exists(".voc.all"))		voc.all = .voc.all			else .voc.all =		voc.all		= read.csv('/Volumes/lldeep/data/VOCData/VOCDataCombined.txt', sep='\t', row.names=1)
if (exists(".couple.all"))	couple.all = .couple.all	else .couple.all =	couple.all	= read.csv('/Volumes/lldeep/data/VOCData/GenotypeToVOCCoupling.txt', sep='\t', header=F, as.is=T)

#
## Preprocess VOC
#
# strip last 3 columns of VOC ("1000"    "1000.1"  "1000.2")
voc.all = voc.all[, -c(ncol(voc.all) - 2:0)]
colnames(voc.all) = as.numeric(str_sub(colnames(voc.all),2))

#
## Find out for which individuals we have both Pheno and VOC data
#
index.couple.genotype	= which(couple.all[, 1] %in% rownames(p.all))
index.couple.voc		= which(couple.all[, 2] %in% colnames(voc.all))
index.couple.overlap	= Reduce(intersect,  list(index.couple.genotype, index.couple.voc))

# Create couple table with only 'complete' couplings
couple = couple.all[index.couple.overlap, ]

# Create Pheno matrix with only individuals for which we also have VOC data
p = p.all[couple[, 1], ]

# Convert characters to numeric
p = as.matrix(p)
class(p) = "numeric"


# Create VOC matrix corresponding to Pheno matrix
voc.matching.subset = voc.all[, as.character(couple[, 2])]
voc = t(voc.matching.subset)

# Rename VOC_ID to Genotype_ID
rownames(voc) = as.vector(sapply(rownames(voc), function(gID) couple[which(couple[,2] == gID), 1]))

# Stop if p and voc have not equal row-order of Genotype_IDs
stopifnot(all(rownames(p) == rownames(voc)))