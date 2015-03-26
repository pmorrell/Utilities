# R code for genotype quality control
# Written by Peter L. Morrell
# Includes edits by Chaochih Liu
# 2 March 2015 - Frontenac, MO 

# read in genotyping data file to be used
# includes row names for WBDC samples
WBDC_geno_count <- read.delim('~/Dropbox/Barley inversions/WBDC_genotype_count.txt', row.names=1)
# WBDC_geno_count <- read.delim('~/Dropbox/Documents/Work/manuscripts/barley inversions/WBDC_genotype_count.txt', row.names=1)
# read in list of samples used (or more importantly, to eliminate those that will not be used)
# this marker dataset from Fang et al. 2014 G3 - excludes accessions putatively subject to introgression
WBDC_284 <- read.delim('~/Dropbox/Barley Inversions/Table_S1.txt',comment.char='#') 
# WBDC_284 <- read.delim('~/Dropbox/Documents/Work/manuscripts/barley inversions/Table_S1.txt',comment.char='#') 


# define thresholds for missing data, heterozygous proportions, & mininimum minor allele frequency
# remove SNPs with missing data ≥ 25% 
# remove SNPs with observed heterozygosity > 10%
# remove SNPs with diversity = 0%; MAF = 0
miss_proportion <- 0.25
het_proportion <- 0.10
min_freq <- 0

# set project specific files to generic names used for code below
geno_dat <- WBDC_geno_count
included_samples <- WBDC_284[,1]
dim(WBDC_geno_count)

# from here on, code to entirely generalizable
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# only retain genotypes for samples in the list above
geno_dat <- geno_dat[row.names(geno_dat) %in% included_samples,]

# function to omit columns of data that are only NA values
all.na.omit.column <- function(dat.frame) {
mask <- apply(dat.frame,2, function(x) all(is.na(x)))
dat.frame <- dat.frame[!mask]
return(dat.frame)
}

# function to count missing data based on NA count
missing <- function(dat) {
miss <- sum(is.na(dat))
miss <- miss/(length(dat))
return(miss)
}

# function to identify proportion of SNPs that are heterozygous
# here they are all represented as genotype value of "1"
hets <- function(dat) {
het <- dat[dat == "1"]
het <- length(het)
hets_locus <- het/length(na.omit(dat))
return(hets_locus)
}

# below, sort genotype data frame by accessions and SNP name, omit SNPS with only missing data
# sort data frame by row names
geno_dat <- geno_dat[order(row.names(geno_dat)),]

# sort data frame by locus names
geno_dat <- geno_dat[,order(colnames(geno_dat))]

# omit SNPs that are missing only
geno_dat <- all.na.omit.column(geno_dat)

# all of the functions below are run using apply
# remove SNPs with missing data ≥ 25% 
# remove SNPs with observed heterozygosity > 10%
# remove SNPs with diversity = 0%; MAF = 0
miss_rm <- which(apply(geno_dat,2,missing) >= miss_proportion)
het_rm <- which(apply(geno_dat,2,hets) >= het_proportion)
mono_rm <- which(apply(geno_dat,2,sum, na.rm = TRUE) == min_freq)

geno_dat <- geno_dat[,c(-miss_rm,-het_rm,-mono_rm)]
dim(geno_dat)


# from here on, code for Structure file creation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# create duplicate data frame we can convert to Structure format
geno_struct <- geno_dat

 # convert genotypes counts to integers to be used in Structure
 # write the heterozygotes to '1'; shouldn't be many
 toStructure <- function(dat) {
	dat[dat == 1] <- 2
	dat[dat == 0] <- 1
	dat[is.na(dat)] <- -9
	return(dat)
}

geno_struct <- toStructure(geno_struct)
padding <- rep(1,dim(geno_struct)[1])
geno_struct <- cbind(padding,padding,geno_struct)

write.table(geno_struct, file="~/Desktop/WBDC_SNP.stu", quote=FALSE,sep='\t')
dim(geno_struct)

# from here on, shell code for cleaning up Structure file with name padding, etc.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#perl -i.bak -pln -e 's/X11_/11_/g' ~/Desktop/WBDC_SNP.stu
#perl -i.bak -pln -e 's/X12_/12_/g' ~/Desktop/WBDC_SNP.stu
#perl -i.bak -pln -e 's/padding\tpadding\t//g' ~/Desktop/WBDC_SNP.stu

