library(Rsubread)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
gene.df = data.frame(genes(txdb))

############################################################################################
### get fus eclip data


fus1.df = read.table(gzfile("D:/drabh/Dropbox/ePRINT/fus-eclip/fus-eclip-hepg2-hg19.bed.gz"), sep="\t", header=F)
fus2.df = read.table(gzfile("D:/drabh/Dropbox/ePRINT/fus-eclip/fus-eclip-k562-hg19.bed.gz"), sep="\t", header=F)
fus.df = rbind(fus1.df[, 1:6], fus2.df[,1:6])

#write.table(fus.df, file="D:/drabh/Dropbox/ePRINT/FUS.eclip.K562.encode.bed", sep="\t", col.names = F, row.names = F, quote=F)

colnames(fus.df) = c("chr", "start", "end", "name", "score", "strand")
fus.df$name = paste("peak", 1:dim(fus.df)[1], sep="_")
fus.df = fus.df[which(fus.df$chr %in% validchrname),]

fuspk.gr = makeGRangesFromDataFrame(fus.df,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end",
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)

fuspk.ann.gr = annotatePeak(fuspk.gr, 
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

fuspk.ann.df = data.frame(fuspk.ann.gr)
fuspk.ann.df$intragenic = sapply(1:length(fuspk.ann.df$annotation), function(i){ flag_annot(fuspk.ann.df$annotation[i]) } )
fuspk.ann.df$symbol = dfid$hgnc_symbol[match(fuspk.ann.df$geneId, dfid$entrezgene_id)]

fuspk.ann.df = fuspk.ann.df[which(fuspk.ann.df$intragenic==1),]
fuspk.ann.df = na.omit(fuspk.ann.df)
fus.gene.targets = unique(fuspk.ann.df$symbol)

write.table(fuspk.ann.df, "D:/drabh/Dropbox/ePRINT/fus-eclip/fus.eclip.encode.peaks.gene.targets.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(fus.gene.targets, "D:/drabh/Dropbox/ePRINT/fus-eclip/fus.eclip.encode.gene.targets.txt", row.names = F, col.names = F, quote = F)

################################################################################################

fuspk.ann.df = read.table("D:/drabh/Dropbox/ePRINT/fus-eclip/fus.eclip.encode.peaks.gene.targets.txt", sep="\t", header=T, quote="")
fuspk.ann.df = fuspk.ann.df[which(fuspk.ann.df$intragenic==1),]
fuspk.ann.df = na.omit(fuspk.ann.df)
fus.gene.targets = unique(fuspk.ann.df$symbol)


res.df = read.table(file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fus.merge.peak.counts.fusvsneg.deseq.lrt.txt",
                    sep="\t", header=T)
rownames(res.df) = res.df$name

dim(res.df)



pk.ann.df = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.gene.counts.txt", sep="\t", header=T, quote="")
rownames(pk.ann.df) = pk.ann.df$name
pk.ann.df = pk.ann.df[which(pk.ann.df$intragenic == 1),]

dim(pk.ann.df)

eprint.thresh = 50
expr.thresh = 100


tmpvec = apply(pk.ann.df[, c("sample1", "sample2", "sample3", "sample4")], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec > eprint.thresh), ]

tmpvec = apply(pk.ann.df[, c("sample5", "sample6", "sample7", "sample8")], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec > expr.thresh), ]


valid_peaks = intersect(rownames(res.df), rownames(pk.ann.df))

pk.ann.df = pk.ann.df[valid_peaks,]
res.df = res.df[valid_peaks, ]


fus.gene.targets = intersect(fus.gene.targets, pk.ann.df$symbol)

pvalmat = data.frame(matrix(1, nrow=1, ncol=5))

fc.thresh = log2(2.0)
p.thresh = 0.05

#######################################################################################
## down

#res.pk = rownames(head(res.df, top))
res.pk = rownames( res.df[which(res.df$padj < p.thresh & res.df$log2FoldChange < -1*fc.thresh),] )
length(res.pk)

res.gene.targets.dw = unique(pk.ann.df[res.pk, "symbol"])
N = length(unique(pk.ann.df$symbol))
m = length(fus.gene.targets)
n = N - m
(k = length(res.gene.targets.dw))
(q = length( intersect(fus.gene.targets, res.gene.targets.dw) ))
(p = phyper(q-1, m,n,k, lower.tail = F))
pvalmat[1,] = c(N, m, k, q, p)



#######################################################################################
## up

#res.pk = rownames(tail(res.df, top))
res.pk = rownames( res.df[which(res.df$padj < p.thresh & res.df$log2FoldChange > fc.thresh),] )

length(res.pk)

res.gene.targets.up = unique(pk.ann.df[res.pk, "symbol"])
N = length(unique(pk.ann.df$symbol))
m = length(fus.gene.targets)
n = N - m
(k = length(res.gene.targets.up))
(q = length( intersect(fus.gene.targets, res.gene.targets.up) ))
(p = phyper(q-1, m,n,k, lower.tail = F))

pvalmat[2,] = c(N, m, k, q, p)

#######################################################################################
## down only

res.gene.targets.filt.dw = setdiff(res.gene.targets.dw, res.gene.targets.up)
N = length(unique(pk.ann.df$symbol))
m = length(fus.gene.targets)
n = N - m
(k = length(res.gene.targets.filt.dw))
(q = length( intersect(fus.gene.targets, res.gene.targets.filt.dw) ))
(p = phyper(q-1, m,n,k, lower.tail = F))

pvalmat[3,] = c(N, m, k, q, p)

#######################################################################################
## up only

res.gene.targets.filt.up = setdiff(res.gene.targets.up, res.gene.targets.dw)
N = length(unique(pk.ann.df$symbol))
m = length(fus.gene.targets)
n = N - m
(k = length(res.gene.targets.filt.up))
(q = length( intersect(fus.gene.targets, res.gene.targets.filt.up) ))
(p = phyper(q-1, m,n,k, lower.tail = F))

pvalmat[4,] = c(N, m, k, q, p)

#######################################################################################
## all

res.pk = rownames( res.df[which(res.df$padj < p.thresh ), ])
res.pk = intersect(res.pk, rownames(pk.ann.df))


res.gene.targets = unique(pk.ann.df[res.pk, "symbol"])
N = length(unique(pk.ann.df$symbol))
m = length(fus.gene.targets)
n = N - m
(k = length(res.gene.targets))
(q = length( intersect(fus.gene.targets, res.gene.targets) ))
(p = phyper(q-1, m,n,k, lower.tail = F))

pvalmat[5,] = c(N, m, k, q, p)


####

colnames(pvalmat) = c("total_num_genes", "num_fus_targets", "num_pk_targets", "overlap_fus_pk", "hypergeo_pval")
pvalmat$cmp = c("down", "up", "down only", "up only", "all")

write.table(pvalmat, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fus.merge.peak.counts.overlap.pval.threshold.txt", sep="\t",
            row.names = F, col.names = T, quote = F)

####################
####################
##threshold independent analysis

binsize = 500
numsteps = 10

pvalmat = data.frame(matrix(1, nrow=1, ncol=7))

res.df$rankval = -1*log10(res.df$pvalue) * sign(res.df$log2FoldChange)
res.df = res.df[order(res.df$rankval),]

for(stepctr in 1:numsteps){
  
  top = binsize*stepctr
  i = stepctr + 3*(stepctr-1)
  print(c(top, i))
  #######################################################################################
  ## down
  
  res.pk = rownames(head(res.df, top))
  res.gene.targets.dw = unique(pk.ann.df[res.pk, "symbol"])
  N = length(unique(pk.ann.df$symbol))
  m = length(fus.gene.targets)
  n = N - m
  (k = length(res.gene.targets.dw))
  (q = length( intersect(fus.gene.targets, res.gene.targets.dw) ))
  (p = phyper(q-1, m,n,k, lower.tail = F))
  pvalmat[i,] = c(N, m, k, q, p, top, "down")
  
  
  
  #######################################################################################
  ## up
  
  res.pk = rownames(tail(res.df, top))
  length(res.pk)
  
  res.gene.targets.up = unique(pk.ann.df[res.pk, "symbol"])
  N = length(unique(pk.ann.df$symbol))
  m = length(fus.gene.targets)
  n = N - m
  (k = length(res.gene.targets.up))
  (q = length( intersect(fus.gene.targets, res.gene.targets.up) ))
  (p = phyper(q-1, m,n,k, lower.tail = F))
  
  pvalmat[i+1,] = c(N, m, k, q, p, top, "up")
  
  #######################################################################################
  ## down only
  
  res.gene.targets.filt.dw = setdiff(res.gene.targets.dw, res.gene.targets.up)
  N = length(unique(pk.ann.df$symbol))
  m = length(fus.gene.targets)
  n = N - m
  (k = length(res.gene.targets.filt.dw))
  (q = length( intersect(fus.gene.targets, res.gene.targets.filt.dw) ))
  (p = phyper(q-1, m,n,k, lower.tail = F))
  
  pvalmat[i+2,] = c(N, m, k, q, p, top, "downonly")
  
  #######################################################################################
  ## up only
  
  res.gene.targets.filt.up = setdiff(res.gene.targets.up, res.gene.targets.dw)
  N = length(unique(pk.ann.df$symbol))
  m = length(fus.gene.targets)
  n = N - m
  (k = length(res.gene.targets.filt.up))
  (q = length( intersect(fus.gene.targets, res.gene.targets.filt.up) ))
  (p = phyper(q-1, m,n,k, lower.tail = F))
  
  pvalmat[i+3,] = c(N, m, k, q, p, top, "uponly")
  
}

colnames(pvalmat) = c("total_num_genes", "num_fus_targets", "num_pk_targets", "overlap_fus_pk", "hypergeo_pval", "top", "geneset")

write.table(pvalmat, file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fus.merge.peak.counts.overlap.pval.ranked.txt", sep="\t",
            row.names = F, col.names = T, quote = F)



### match fus eclip peaks to eprint peaks

res.pk = rownames( res.df[which(res.df$padj < p.thresh & res.df$log2FoldChange < -1*fc.thresh),] )

fuspk.df = fuspk.ann.df[fuspk.ann.df$intragenic == 1, c("seqnames", "start", "end", "name", "score", "strand")]
colnames(fuspk.df)[1]  ="chr"

pkdf = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.eprint.merge.sort.bed", sep="\t", header=F)
colnames(pkdf) = c("chr", "start", "end", "name", "score", "strand")
pkdf$name = paste("peak", 1:dim(pkdf)[1], sep="_")
rownames(pkdf) = pkdf$name

pkdf = pkdf[res.pk, ]


pkcount.gr = makeGRangesFromDataFrame(pkcount.df,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field="chr",
                                      start.field="start",
                                      end.field="end",
                                      strand.field="strand",
                                      starts.in.df.are.0based=FALSE)

pk.gr = makeGRangesFromDataFrame(pkdf,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field="chr",
                                 start.field="start",
                                 end.field="end",
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)

fus.gr = makeGRangesFromDataFrame(fuspk.df,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="chr",
                                  start.field="start",
                                  end.field="end",
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)

target.gr.ls = list()
target.gr.ls[[1]] = fus.gr

pk.overlap.df = enrichPeakOverlap(pk.gr, target.gr.ls, txdb, nShuffle = 100)

################################################################################################