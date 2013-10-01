# Purpose: merge full data set with subset (LLDeep ID should also be present in couple.all)
p = read.csv('/Volumes/lldeep/data/PhenotypeData2017-08-06BloodCellCountParameters.txt', sep='\t', row.names=1)
p.all = read.csv('/Volumes/lldeep/data/phenotypes.txt', sep='\t', as.is=T, na.strings=c("#N/A","#NULL!")) # row headers are in first column, in this matrix!
couple.all	= read.csv('/Volumes/lldeep/data/VOCData/GenotypeToVOCCoupling.txt', sep='\t', header=F, as.is=T)

# BMI is not present in p.all
#colnames(p)[which(!colnames(p) %in% colnames(p.all))]
#[1] "BMI"

index.new.individuals = which(!is.na(p.all[,1]) & !p.all[,1] %in% rownames(p) & is.element(p.all[,1], couple.all[,1]))
p.new = p.all[index.new.individuals, ]

# Add BMI
p.new = cbind(p.new, BMI = p.new[, "GEWICHT"] / (p.new[, "LENGTE"]/100)^2)

# Create matrix in same style as target matrix p (rownames and column selection and order)
rownames(p.new) = p.new[,1]
p.new = p.new[, colnames(p)]


# Convert GENDER and NUCHTER

		#> p.all[5,c("LLDEEPID", "NUCHTER", "GESLACHT")]
		#     LLDEEPID NUCHTER GESLACHT
		#5 LLDeep_0005      Ja      Man

		# should be converted to
	
		#> p["LLDeep_0005",c("NUCHTER", "GESLACHT")]
		#            NUCHTER GESLACHT
		#LLDeep_0005       1        1
	
	# NUCHTER 1=Ja, 0=NEE
	# GESLACHT 1=Man, 2=Vrouw

p.new["Ja" == p.new[, "NUCHTER"], "NUCHTER"] = 1
p.new["Nee" == p.new[, "NUCHTER"], "NUCHTER"] = 0

p.new["Man" == p.new[, "GESLACHT"], "GESLACHT"] = 1
p.new["Vrouw" == p.new[, "GESLACHT"], "GESLACHT"] = 2

# Now we can merge!
p = rbind(p, p.new)

dput(p, file="/Volumes/lldeep/data/phenotypeDataBloodCellCountParametersAll.RData")
#write.csv(p, file="/Volumes/lldeep/data/phenotypeDataBloodCellCountParametersAll.txt", sep="\t", row.names=T, col.names=T)