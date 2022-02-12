
library(transite)
library(BSgenome.Hsapiens.UCSC.hg19)
library(universalmotif)


setwd("D:/software/R/Rtemp")
min.pk.len = 50

res.df = read.table(file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/diffpeaks/fus.merge.peak.counts.fusvsneg.deseq.lrt.txt",
                    header=T, sep="\t")
rownames(res.df) = res.df$name
dim(res.df)

pk.ann.df = read.table("D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/fus.merge.peak.eprint.counts.gene.counts.txt", 
                       sep="\t", header=T, quote="")
rownames(pk.ann.df) = pk.ann.df$name
dim(pk.ann.df)

pk.ann.df = pk.ann.df[res.df$name, ]
pk.ann.df = cbind(pk.ann.df, res.df[pk.ann.df$name, ])

dim(pk.ann.df)

eprint.thresh = 50
expr.thresh = 100
p.thresh = 0.05
fc.thresh = log2(2)



tmpvec = apply(pk.ann.df[, c("sample1", "sample2", "sample3", "sample4")], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec > eprint.thresh), ]

tmpvec = apply(pk.ann.df[, c("sample5", "sample6", "sample7", "sample8")], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec > expr.thresh), ]


eprint.sampleid = c("sample1", "sample2", "sample3", "sample4")
input.sampleid = c("sample5", "sample6", "sample7", "sample8")

tmpvec = apply(pk.ann.df[, eprint.sampleid], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec > eprint.thresh), ]

tmpvec = apply(pk.ann.df[, input.sampleid], 1, sum)
pk.ann.df = pk.ann.df[which(tmpvec > expr.thresh), ]


rbp.m = get_motifs()
lenvec = sapply(1:length(rbp.m), function(i){ get_width(rbp.m[[i]])})

if(min.pk.len < max(lenvec)){ min.pk.len = max(lenvec)}

pkind = which( (pk.ann.df$end-pk.ann.df$start + 1) < min.pk.len)

tmpdf = sapply(pkind, function(i){ center = (pk.ann.df$start[i] + pk.ann.df$end[i])/2;
                                   start = as.integer( center - min.pk.len/2);
                                   end = as.integer( center + min.pk.len/2)
                                   return(c(start,end))
                                 })
tmpdf = t(tmpdf)
pk.ann.df[pkind, c("start", "end")]  = tmpdf
  
  
pk.gr = makeGRangesFromDataFrame(pk.ann.df[, c("seqnames", "start", "end", "name", "strand")],
                                 keep.extra.columns=TRUE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field="seqnames",
                                 start.field="start",
                                 end.field="end",
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)

pk.seq = getSeq(Hsapiens, pk.gr)
pk.seq.vec = gsub("T", "U", as.character(pk.seq)) 
pk.rna.seq = BStringSet(pk.seq.vec)


motif.hits.ls = list()
  
for(i in 1:length(rbp.m))
{
  print(i)
  mobj = rbp.m[[i]]
  
  name = get_rbps(mobj)
  id = get_id(mobj)
  
  motif.matrix =  t(as.matrix(get_motif_matrix(mobj)))
  motif = create_motif(motif.matrix, alphabet = "RNA", name = name, pseudocount = 1, strand = "+", type="PWM")

	mscan = scan_sequences(motif, pk.rna.seq, threshold.type = "pvalue", threshold = 0.0002, no.overlaps=F, RC=F)
	motif.hits.ls[[i]] = list(name = name, id = id, mhits = mscan) 
  rm(name, id, mscan)
}

saveRDS(motif.hits.ls, "D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/motifs/fusvsneg.counts.motif.hits.ls")

######################
## clean up motif names

#m.id = sapply(1:length(rbp.m), function(i){ paste0(get_id(rbp.m[[i]]), collapse=";")   })
#rbp.names = sapply(1:length(rbp.m), function(i){ paste0(get_rbps(rbp.m[[i]]), collapse=";")   })
#write.table(data.frame(motif=m.id, rbp=rbp.names), "D:/drabh/Dropbox/ePRINT/transite-rbp-motif-names.txt", sep="\t", row.names = F, col.names = T, quote = F)
##########################


statvec = -1*log10(pk.ann.df$pvalue) * sign(pk.ann.df$log2FoldChange)
names(statvec) = pk.ann.df$name
statvec = sort(statvec)

peakhits.ls = list()

for(i in 1:length(motif.hits.ls)){

  id = paste0(motif.hits.ls[[i]]$id, collapse = ":")
  peakhits.ls[[id]] = unique(motif.hits.ls[[i]]$mhits$sequence)
}

library(fgsea)

set.seed(12345)
fgsea.pos.res = fgseaMultilevel(peakhits.ls, statvec, minSize = 100, maxSize = 100000, eps=0, scoreType = "pos")

fgsea.pos.df = data.frame(fgsea.pos.res)
rownames(fgsea.pos.df) = fgsea.pos.df$pathway

pvalscore = -1*log10(fgsea.pos.df$padj)*sign(fgsea.pos.df$NES)
names(pvalscore) = fgsea.pos.df$pathway
pvalscore = rev(sort(pvalscore))

fgsea.pos.df = fgsea.pos.df[names(pvalscore), ]

rbp.df = read.table("D:/drabh/Dropbox/ePRINT/transite-rbp-motif-hgnc.txt", sep="\t", header=T)
fgsea.pos.df$rbp = rbp.df$hgnc[match(fgsea.pos.df$pathway, rbp.df$motif)]

write.table(fgsea.pos.df[, c(1:7,9)], file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/motifs/fusvsneg.counts.pos.knownmotif.enrich.txt", 
            sep="\t", row.names = F, col.names = T, quote=F)

###
set.seed(12345)
fgsea.neg.res = fgseaMultilevel(peakhits.ls, statvec, minSize = 100, maxSize = 100000, eps=0, scoreType = "neg")

fgsea.neg.df = data.frame(fgsea.neg.res)
rownames(fgsea.neg.df) = fgsea.neg.df$pathway

pvalscore = -1*log10(fgsea.neg.df$padj)*sign(fgsea.neg.df$NES)
names(pvalscore) = fgsea.neg.df$pathway
pvalscore = sort(pvalscore)

fgsea.neg.df = fgsea.neg.df[names(pvalscore), ]

rbp.df = read.table("D:/drabh/Dropbox/ePRINT/transite-rbp-motif-hgnc.txt", sep="\t", header=T)
fgsea.neg.df$rbp = rbp.df$hgnc[match(fgsea.neg.df$pathway, rbp.df$motif)]

write.table(fgsea.neg.df[, c(1:7,9)], file="D:/drabh/Dropbox/ePRINT/V0070-peaks/first_read_bundled/motifs/fusvsneg.counts.neg.knownmotif.enrich.txt", 
            sep="\t", row.names = F, col.names = T, quote=F)


###############################################################################################
### kmer enrichment

library(gtools)


match_kmers_per_sequence<-function(sequences, klen){

counts.df = oligonucleotideFrequency(RNAStringSet(sequences), klen)
counts = apply(counts.df,2,sum)
names(counts) = colnames(counts.df)

return(counts)

}

generate_kmer_counts<-function(sequences, klen){

  
  print(length(sequences))
  print(klen)
  
if(klen <= 7){
  
if(length(sequences) > 1000){ binsize = 1000 }
if(length(sequences) <= 1000){ binsize = length(sequences) }
  
}
  
if(klen > 8){
    
    if(length(sequences) > 100){ binsize = 100 }
    if(length(sequences) <= 100){ binsize = length(sequences) }
    
}

tmpdf = permutations(n=4, r=klen, v = c("A", "C", "G", "U"), repeats.allowed = T )

kmervec = sapply(1:dim(tmpdf)[1], function(i){ paste0(tmpdf[i,], collapse="")  })
  
binctr = 0
kmervec.counts = rep(0, length(kmervec))
names(kmervec.counts) = kmervec

while(binctr < length(sequences)){

   if((binctr %% 10000) == 0){ print(binctr)}
   if((binctr+binsize) > length(sequences)){

     seqvec = sequences[(binctr+1): length(sequences)]
   }else{

     seqvec = sequences[(binctr+1): (binctr + binsize)]
   }
  
   counts = match_kmers_per_sequence(seqvec, klen)   
   kmervec.counts[names(counts)] = kmervec.counts[names(counts)] + counts
   binctr = binctr + binsize                         
}
  
 return(kmervec.counts)

}


pk.up = pk.ann.df$name[which(pk.ann.df$padj < 0.1 & pk.ann.df$log2FoldChange > 0)]
pk.dw = pk.ann.df$name[which(pk.ann.df$padj < 0.1 & pk.ann.df$log2FoldChange < 0)]

#pk.ann.df$rankval = -1*log10(pk.ann.df$pvalue) * sign(pk.ann.df$log2FoldChange)
#pk.ann.df = pk.ann.df[order(pk.ann.df$rankval), ]
#pk.dw = rownames(head(pk.ann.df, 1000))

K = 6
kmervec.bg = generate_kmer_counts(pk.rna.seq, K)
kmervec.fg.dw = generate_kmer_counts(pk.rna.seq[pk.dw], K)
kmervec.fg.up = generate_kmer_counts(pk.rna.seq[pk.up], K)



#kmervec.bg = generate_kmers(pk.rna.seq, k=K)
#kmervec.fg = generate_kmers(pk.rna.seq[pk.dw], k=K)



kmer.df = compute_kmer_enrichment(kmervec.fg, kmervec.bg)
kmers = names(kmervec.fg)
kmers = gsub("T", "U", kmers)

kmer.df$kmers = kmers
rownames(kmer.df) = kmers

kmer.df = kmer.df[rev(order(kmer.df$enrichment)),]
sort(kmer.df$kmers[which(kmer.df$adj_p_value < 0.1 & kmer.df$enrichment > 1)])


#################################################################################################
