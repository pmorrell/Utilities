# remove all existing files
rm(list = ls())

library(LDheatmap)
library(genetics)

# read in SNP data set
SNPS <- read.table('~/R/SNP_ALLCAPs_3528_Cleaned_AC.txt', header=T)

# read in map file
map <- read.table("~/R/map.txt", header=T)
map_names <- colnames(map)

#spring_2row <- subset(SNPS, Growth_Habit == 'spring' && Row_Type == 2)

# pad names of SNPs in map so they match values in data column names
Xnames <- rep('X', dim(map)[1])
Marker_Name <- paste(Xnames,map$Marker,sep='')
map <- cbind(Marker_Name,map$Chr,map$Pos)
colnames(map) <- map_names
map <- as.data.frame(map)

# subset map and only use chromosome 7
chromo <- subset(map, Chr == 7)
chromo <- as.data.frame(chromo)

# retain genotype data from mapped markers 
chromo_SNPS <- SNPS[,colnames(SNPS) %in% as.character(chromo$Marker)]

# sort genotype data
#chromo_SNPS <- chromo_SNPS[,order(as.character(chromo$Marker))]

# using makeGenotypes from the genetics library to turn into genotype object
genos <- makeGenotypes(chromo_SNPS,sep='')

pdf(file='~/R/LG_7H_r2.pdf',width=10,height=10)
LDheatmap(genos, LDmeasure="r", distances='genetic', add.map=F, flip=T)
dev.off()
