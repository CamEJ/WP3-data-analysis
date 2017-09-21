## Making biom type file with abundunce, taxonomy and metadata in one. 
# then read it out into a csv so dont have to read in heavy mothur files each time 


library("phyloseq")
library("plyr")
library("vegan")
library("grid")
library("directlabels")
library("knitr")
library("clustsig")
library("ellipse")

# will use subsampled file with 70 samples and 28815 seq subsample depth.

sharedsubFile = read.table('stability.opti_mcc.0.03.subsample.shared')
sharedsubFile = t(sharedsubFile)
rownames(sharedsubFile) = sharedsubFile[,1]
colnames(sharedsubFile) = sharedsubFile[2,]
dim(sharedsubFile) # [1] 87682    90
sharedsubFile = sharedsubFile[,2:90]
sharedsubFile = sharedsubFile[4:87682,]
class(sharedsubFile) <- "numeric"
head(sharedsubFile)

taxFile = read.table('stability.cons.taxonomy', header=T, sep='\t')
head(taxFile)
fix = taxFile[,3]

fix = as.data.frame(fix)
write.table(fix, 'taxOnly.csv', row.names=FALSE,col.names=FALSE)
# open this in excel specifying ; as separators and add phylum, kingdom etc as col names
# select all and copy and paste into git/notepad and save as .txt file (HeadedConsTax.txt)
# read it back in

taxo = read.table('HeadedConsTax.txt', header=T, sep='\t')
# join original two cols from cons.tax file
library(tibble)
taxFile = taxFile[,1:2]
taxo1 = add_column(taxo, taxFile$OTU, .before = 1) # rename
taxo2 = add_column(taxo1, taxFile$Size, .before = 2) # rename

taxFile = taxo2
colnames(taxFile)[2] <- "Size"
colnames(taxFile)[1] <- "OTU"
head(taxFile)
# now it's ready, carry on with normal stuff..

rownames(taxFile) = taxFile[,1]
taxFile = taxFile[,2:8]
taxFile = as.matrix(taxFile)
head(taxFile)

## import metadata file
metaFile = read.table('MetaDataSubSampled.txt', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:6]
head(metaFile)

#metaFile$SamplingTime <- factor(metaFile$SamplingTime, levels = metaFile$SamplingTime)

## Create phyloseq object

#OTU = otu_table(sharedFile, taxa_are_rows = TRUE)
OTUsub = otu_table(sharedsubFile, taxa_are_rows = TRUE)
TAX = tax_table(taxFile)
META = sample_data(metaFile)
#physeq = phyloseq(OTU, TAX, META)
physeqSub = phyloseq(OTUsub, TAX, META)

## Get rid of any OTUs not present in any samples and get relative abundance

microSub <- prune_taxa(taxa_sums(physeqSub) > 0, physeqSub)
# remove any taxa whose row sums comes to 
# less than 2 total
trimSub <- prune_taxa(taxa_sums(physeqSub) > 1, physeqSub) 

microSubRel = transform_sample_counts(microSub, function(x) x / sum(x) )
trimSubRel = transform_sample_counts(trimSub, function(x) x / sum(x) )

## readout phyloseq object, untrimmed
dat <- psmelt(microSubRel)
write.csv(dat, file='untrimmedPhyloseqOTU_RelAbund.csv')

dat <- psmelt(trimSub)
write.csv(dat, file='NoSingPhyloseqOTU_Counts.csv')

dat <- psmelt(trimSubRel)
write.csv(dat, file='NoSingPhyloseqOTU_RelAbund.csv')

## that csv can then be read again and used to make rel abund plot etc 
## using ggplot
