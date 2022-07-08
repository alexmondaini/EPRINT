##
library(Rsubread)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
##

setwd("/my_dir/EPRINT")

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
gene.df = data.frame(genes(txdb))
#### convert to SAF format
gene.df = gene.df[, c(6,1,2,3,5)]
colnames(gene.df) = c("GeneID",	"Chr",	"Start",	"End",	"Strand")

validchrname = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",  
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                 "chr20", "chr21", "chr22", "chrX", "chrY")

gene.df = gene.df[which(gene.df$Chr %in% validchrname),]

### Generate counts per gene from input bam files
bam_files = paste0("data/bam/",list.files("data/bam/"))

### input bam files mapped to genes
gene.fc = featureCounts(as.character(bam_files), annot.ext = gene.df, isPairedEnd = T, strandSpecific = 1, 
                   countMultiMappingReads = F, requireBothEndsMapped = T, countChimericFragments = F)


genecount.df = cbind(gene.fc$annotation, gene.fc$counts)
colnames(genecount.df) = c(colnames(genecount.df)[1:6],paste0('eprint_',seq(4)),paste0('input_',seq(4)))

write.table(genecount.df, file="data/fus-eprint/fus.eprint_input_gene_counts.txt", sep="\t", row.names = F, col.names = T, quote = F)

### eprint bam files mapped to genes

filename.vec = bam_files$V1[bam_files$V3=="eprint"]
gene.fc = featureCounts(as.character(filename.vec), annot.ext = gene.df, isPairedEnd = T, strandSpecific = 1, 
                        countMultiMappingReads = F, requireBothEndsMapped = T, countChimericFragments = F)


genecount.df = cbind(gene.fc$annotation, gene.fc$counts)
colnames(genecount.df) = c(colnames(genecount.df)[1:6], bam_files$V2[bam_files$V3=="eprint"])

write.table(genecount.df, file="fus.eprint.gene.counts.txt", sep="\t", row.names = F, col.names = T, quote = F)


### eprint bam files mapped to peaks

filename.vec = bam_files$V1[bam_files$V3=="eprint"]

pkdf = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.eprint.merge.sort.bed", sep="\t", header=F)
colnames(pkdf) = c("chr", "start", "end", "name", "score", "strand")
pkdf$name = paste("peak", 1:dim(pkdf)[1], sep="_")
pkdf = pkdf[, c("name", "chr", "start", "end", "strand")]
colnames(pkdf) = c("GeneID", "Chr", "Start", "End", "Strand")


pk.fc = featureCounts(as.character(filename.vec), annot.ext = pkdf, isPairedEnd = T, strandSpecific = 1, 
                        countMultiMappingReads = F, requireBothEndsMapped = T, countChimericFragments = F)


pkcount.df = cbind(pk.fc$annotation, pk.fc$counts)
colnames(pkcount.df) = c(colnames(pkcount.df)[1:6], bam_files$V2[bam_files$V3=="eprint"])
colnames(pkcount.df)[1:6] = c("name", "chr", "start", "end", "strand", "length")
write.table(pkcount.df, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.txt", sep="\t", row.names = F, col.names = T, quote = F)


######################################################################################################
### map peaks to genes

pkdf = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.eprint.merge.sort.bed", sep="\t", header=F)
colnames(pkdf) = c("chr", "start", "end", "name", "score", "strand")
pkdf$name = paste("peak", 1:dim(pkdf)[1], sep="_")



pk.gr = makeGRangesFromDataFrame(pkdf,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="chr",
                                  start.field="start",
                                  end.field="end",
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)

pk.ann.gr = annotatePeak(pk.gr, 
                          tssRegion = c(0, 0),
                          TxDb = txdb,
                          level = "gene",
                          assignGenomicAnnotation = TRUE,
                          genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron",
                                                        "Downstream", "Intergenic", "Promoter"),
                          annoDb = NULL,
                          addFlankGeneInfo = FALSE,
                          flankDistance = 0,
                          sameStrand = TRUE,
                          ignoreOverlap = FALSE,
                          ignoreUpstream = FALSE,
                          ignoreDownstream = FALSE,
                          overlap = "all",
                          verbose = TRUE)

pk.ann.df = data.frame(pk.ann.gr)

write.table(pk.ann.df, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.gene.ann.txt", sep="\t", 
            row.names = F, col.names = T, quote = F)

######################################################################################################
### merge all information into one dataset

genecount.df = read.table("fus.input.gene.counts.txt", sep="\t", header=T)
pkcount.df = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.txt", sep="\t", header=T)
pk.ann.df = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.gene.ann.txt", sep="\t", header=T, quote="")

#colnames(genecount.df) = c("GeneID", "Chr", "Start", "End", "Strand", "Length", "sample5", "sample6", "sample7", "sample8")
rownames(pk.ann.df) = pk.ann.df$name

pk.ann.df[, c("sample1", "sample2", "sample3", "sample4")] = pkcount.df[match(pk.ann.df$name, pkcount.df$name), c("sample1", "sample2", "sample3", "sample4")]
pk.ann.df[, c("sample5", "sample6", "sample7", "sample8")] = genecount.df[match(pk.ann.df$geneId, genecount.df$GeneID), c("sample5", "sample6", "sample7", "sample8")]

library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
mapped_geneid <- mappedkeys(x)
mapped_sym <- unlist(as.list(x[mapped_geneid]))
dfid = data.frame(entrezgene_id = mapped_geneid, hgnc_symbol = mapped_sym)

pk.ann.df$symbol = dfid$hgnc_symbol[match(pk.ann.df$geneId, dfid$entrezgene_id)]


flag_annot<-function(annotation){
    flag = 0
    if( length(grep("UTR", annotation, ignore.case = TRUE)) ){ flag = 1 }
    if( length(grep("Intron", annotation, ignore.case = TRUE)) ){ flag = 1 }
    if( length(grep("Exon", annotation, ignore.case = TRUE)) ){ flag = 1 }
    return(flag)
}

pk.ann.df$intragenic = sapply(1:length(pk.ann.df$annotation), function(i){ flag_annot(pk.ann.df$annotation[i]) } )


write.table(pk.ann.df, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.gene.counts.txt", sep="\t", 
            row.names = F, col.names = T, quote = F)

###########################################################################################################
### get differential peaks

## get input sizefactors
genecount.df = read.table("data/fus-eprint/fus.input.gene.counts.txt", sep="\t", header=T)

genecount.mat = genecount.df[, c("sample5", "sample6", "sample7", "sample8")]
rownames(genecount.mat) = genecount.df$GeneID

tmpvec = apply(genecount.mat[, c("sample5", "sample6", "sample7", "sample8")], 1, sum)
genecount.mat = genecount.mat[which(tmpvec >= 10), ]

coldata = data.frame(treatment = c("siNEG", "siNEG", "siFUS", "siFUS"), sampletype = c("input", "input", "input", "input") )
rownames(coldata) = c("sample5", "sample6", "sample7", "sample8")
coldata$batch = c(1,2,1,2)

coldata$treatment = factor(coldata$treatment)
coldata$sampletype = factor(coldata$sampletype)
coldata$batch = factor(coldata$batch)



dds <- DESeqDataSetFromMatrix(countData = genecount.mat,
                              colData = coldata,
                              design = ~ batch + treatment)

dds = estimateSizeFactors(dds)
coldata.df = data.frame(colData(dds))
input.sizefactors = coldata.df$sizeFactor


vsd = vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("treatment", "sampletype"))
rm(vsd, sampleDistMatrix, sampleDists)

dds = DESeq(dds)
res = results(dds, contrast = c("treatment", "siFUS", "siNEG"))

res.df = data.frame(res)
res.df$symbol = pk.ann.df$symbol[match(rownames(res.df), pk.ann.df$geneId)]

write.table(res.df, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fusvsneg.input.gene.deseq.txt",
            sep="\t", row.names = F, col.names = T, quote=F)


rm(genecount.mat, coldata.df, coldata, dds, res, res.df)

### get pk size factors

pkcount.df = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.txt", sep="\t", header=T)
pkcount.mat = pkcount.df[, c("sample1", "sample2", "sample3", "sample4")]
rownames(pkcount.mat) = pkcount.df$name

tmpvec = apply(pkcount.mat[, c("sample1", "sample2", "sample3", "sample4")], 1, sum)
pkcount.mat = pkcount.mat[which(tmpvec >= 10), ]

coldata = data.frame(treatment = c("siNEG", "siNEG", "siFUS", "siFUS"), sampletype = c("eprint", "eprint", "eprint", "eprint") )
rownames(coldata) = c("sample1", "sample2", "sample3", "sample4")
coldata$batch = c(1,2,1,2)
coldata$sampleid = rownames(coldata)

coldata$treatment = factor(coldata$treatment)
coldata$sampletype = factor(coldata$sampletype)
coldata$batch = factor(coldata$batch)

dds <- DESeqDataSetFromMatrix(countData = pkcount.mat,
                              colData = coldata,
                              design = ~ batch + treatment)

dds = estimateSizeFactors(dds)
coldata.df = data.frame(colData(dds))
eprint.sizefactors = coldata.df$sizeFactor


vsd = vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))


library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("sampleid"))
rm(vsd, sampleDistMatrix, sampleDists)

rm(pkcount.mat, coldata.df, coldata, dds)

################################################
pk.ann.df = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.gene.counts.txt", sep="\t", header=T, quote="")

pk.ann.df = pk.ann.df[which(pk.ann.df$intragenic==1),]
rownames(pk.ann.df) = pk.ann.df$name

tmpvec = apply(pk.ann.df[, c("sample1", "sample2", "sample3", "sample4")], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec >= 10), ]

tmpvec = apply(pk.ann.df[, c("sample5", "sample6", "sample7", "sample8")], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec >= 10), ]

eprint.mat = pk.ann.df[, c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8")]
rownames(eprint.mat) = pk.ann.df$name

coldata = data.frame(treatment = c("siNEG", "siNEG", "siFUS", "siFUS", "siNEG", "siNEG", "siFUS", "siFUS"), 
                     sampletype = c("eprint", "eprint", "eprint", "eprint", "input", "input", "input", "input") )
rownames(coldata) =  c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8")

coldata$batch = c(1,2,1,2,1,2,1,2)

coldata$sampletype = factor(coldata$sampletype, levels = c("input", "eprint"))
coldata$treatment = factor(coldata$treatment, levels = c("siNEG", "siFUS"))
coldata$batch = factor(coldata$batch, levels = c("1", "2"))

dds <- DESeqDataSetFromMatrix(countData = eprint.mat,
                              colData = coldata,
                              design = ~ batch + treatment + sampletype + treatment:sampletype)

sizeFactors(dds) = c(eprint.sizefactors, input.sizefactors)
dds = estimateDispersions(dds)
#plotDispEsts(dds)
dds = nbinomLRT(dds, reduced = ~ batch + treatment + sampletype)

vsd = vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("treatment", "sampletype"))
rm(vsd, sampleDistMatrix, sampleDists)
##############################################################################################

res = results(dds, name="treatmentsiFUS.sampletypeeprint")

res.df = data.frame(res)
res.df$name = rownames(res.df)

write.table(res.df, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fus.merge.peak.counts.fusvsneg.deseq.lrt.txt",
            sep="\t", row.names = F, col.names = T, quote=F)

saveRDS(dds, "D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fus.merge.peak.counts.fusvsneg.deseq.lrt.dds")
rm(dds)
#############################################################################################################

