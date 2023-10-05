#BiocManager::install("GenomicFiles")
#install.packages("BiocManager")
#install.packages("GenomicRanges")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("GenomicAlignments")
# install.packages("wavelets")
# install.packages("rtracklayer")

library(GenomicAlignments)
library(DESeq2)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicFiles)
library(wavelets)
library(GenomicAlignments)
library(BiocManager)
library(rtracklayer)
library(Rsamtools)
BiocManager::install("DESeq")
library(DESeq)
library(edgeR)





#####DESEQ2######



#load the 4 bam files into R 
bamfiles = BamFileList(dir(pattern=".bam"))


#to load GTF file 
gtf_file <- "Arabidopsis_thaliana.TAIR10.56.gtf"
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

# Extract the exons or other regions of interest
exons <- exonsBy(txdb, by="gene")

#load our metadata
metadata <- read.csv("metadata.csv")


#counting the reads
se = summarizeOverlaps(features=exons, 
                       reads=bamfiles, mode="Union", singleEnd=TRUE, 
                       ignore.strand=TRUE ) 



counts = assay(se) 

countgenenames = gsub("[.][1234567890]", "", row.names(counts)) 

rownames(counts)=countgenenames 

cds = DESeqDataSetFromMatrix(countData=counts,
colData=metadata,
design= ~ condition)


#to normalize the count data, by calculated size factors for each sample
cds = estimateSizeFactors(cds)

cds = estimateDispersions(cds)



dge <- DGEList(counts = counts, group = metadata$condition)
dge <- estimateDisp(dge)




































































































