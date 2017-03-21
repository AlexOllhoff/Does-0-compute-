# Demo GWAS in the BRIDG using the R NAM package
# Does[0]compute? March 21, 2017
# Alex Ollhoff

library(NAM)

# Load genotypic data
# ref allele = 2 = major allele
# 0 = minor allele, this allows minor alleles to have different effects if stratification is provided
gen = genos_forviz <- read.csv("~/Desktop/genos_GxE.csv", header = T)
gen = data.frame(gen)

# Load phenotypic and family data
# phe1 = pheno in col 3, phe2 = fam in col 2
y = read.csv("~/Desktop/phenos_GxE.csv", header = TRUE)[,-1]
y = data.frame(y)

# Set the number of families to select
FAMs <- 2

# Store family numbers
family_names<-c(1:88)

# Select BRIDG families without replacement
sample_fams <- family_names[sample(1:length(family_names), FAMs, replace = F)]
  
  # Filter for individual phenotypes present in the selected families 
  sampled_phenotypes<-y[(y[,3] %in% sample_fams),]
  
  # Filter for genotypes present in the sample
  sampled_genotypes<-gen[(gen[,1] %in% sampled_phenotypes$line_name),]
  
  # Remove sample names 
  gen_tomap <- sampled_genotypes[,-1]
  
  # QC for minor allele frequency
  # Larger threshold will retain more common variants, more stringent threshold will include more uncommon variants
  adjusted_genotypes <- snpQC(gen = gen_tomap, MAF = 0.05, impute = FALSE)
  rownames(adjusted_genotypes) <- sampled_phenotypes$line_name
  dim(adjusted_genotypes) # See how many markers have MAF > 0.05
  
  # Get the number of markers on each chromosome
  chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
  class(chr) # should be an integer
  
  # Run GWAS on your sample of BRIDG families
  # Include phenotypes, genotypes, stratification term for families, and SNPs/chromosome
  my_gwas = gwas2(sampled_phenotypes$BLUE,adjusted_genotypes,sampled_phenotypes$family,chr)
  
  # Plot your GWAS results in a manhattan plot
  # Marker order is not quite right here but it works for discussion
  # If there's no line you have no QTL
  plot( my_gwas, alpha = 0.05 / (ncol(adjusted_genotypes) * (1-0.05)))

  # Calculate the 0.05 false discovery rate threshold
  FDR = 0.05
  THR = -log10(FDR / (ncol(adjusted_genotypes) * (1-FDR)))
  THR
  
  # How many markers are above the FDR threshold?
  w = which((my_gwas$PolyTest$pval) > (THR))
  length(w)
  
###